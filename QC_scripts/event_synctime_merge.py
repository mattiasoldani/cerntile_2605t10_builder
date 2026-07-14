#!/usr/bin/env python3
import sys
import getopt
from pathlib import Path
import numpy as np
#from _curses import A_RIGHT
#from more_itertools.more import difference

NAN = 'NaN'


def usage(exit_code=1):
    prog = Path(sys.argv[0]).name
    print(
        f'Usage: python3 {prog} --file1 FILE --file2 FILE --col1 COL --col2 COL '
        '... [-o OUTPUT]',
        file=sys.stderr,
    )
    print('', file=sys.stderr)
    print('  --file1        first tab-separated file (sorted by merge column)', file=sys.stderr)
    print('  --file2        second tab-separated file (sorted by merge column)', file=sys.stderr)
    print('  --col1         1-based merge column in file1', file=sys.stderr)
    print('  --col2         1-based merge column in file2', file=sys.stderr)
    print('  --scale1       scale to convert col1 to ns', file=sys.stderr)
    print('  --scale2       scale to convert col2 to ns', file=sys.stderr)
    # print('  --shift1       manual shift applied file1 (ns; default: automatic from 1st evt)', file=sys.stderr)
    # print('  --shift2       manual shift applied file2 (ns; default: automatic from 1st evt)', file=sys.stderr)
    print('  --jitter       maximum accepted time jitter (ns; default: 1000)', file=sys.stderr)
    print('  --iteration-steps   to how many ranges will the search be split, >4 (default: 8)', file=sys.stderr)
    print('  --iteration-overlap   to how many ranges will the search be split (0.5-1.0, default 0.6)', file=sys.stderr)
    print('  --iteration-debug   increased verbosity of shift finding algorithm', file=sys.stderr)
    print('  --print-skipped  prints unmatched events with NaNs (default: not used)', file=sys.stderr)
    print('  --initial-range  initial range of offset (to-from) in [s] (default: 1000)', file=sys.stderr)
    print('  --initial-offset initial offset in [s] (default: 0)', file=sys.stderr)
    print('  -o, --output   output file                                     (default: stdout)', file=sys.stderr)
    print('  -h, --help     show this help', file=sys.stderr)
    raise SystemExit(exit_code)


def strip_comment(line):
    return line.split('#', 1)[0].rstrip(' \r\n')


def open_output(filename):
    if not filename or filename == '-':
        return sys.stdout
    return open(filename, 'w', encoding='utf-8', newline='')


def read_next_row(state):
    for raw in state['fh']:
        state['line_no'] += 1
        line = strip_comment(raw)
        if not line.strip():
            continue

        fields = line.split('\t')

        if state['ncols'] is None:
            state['ncols'] = len(fields)
            if state['col_idx'] < 0 or state['col_idx'] >= state['ncols']:
                raise ValueError(
                    f'{state["label"]}: merge column {state["col_idx"] + 1} is out of range '
                    f'for a line with {state["ncols"]} columns in "{state["filename"]}", '
                    f'line {state["line_no"]}'
                )
        elif len(fields) != state['ncols']:
            raise ValueError(
                f'{state["label"]}: inconsistent number of columns in "{state["filename"]}", '
                f'line {state["line_no"]}: expected {state["ncols"]}, got {len(fields)}'
            )

        key_text = fields[state['col_idx']].strip()
        try:
            key = float(key_text)
        except ValueError as exc:
            raise ValueError(
                f'{state["label"]}: cannot parse merge key "{key_text}" '
                f'in "{state["filename"]}", line {state["line_no"]}'
            ) from exc

        return fields, key

    return None, None


def init_reader(filename, col_idx, label):
    state = {
        'filename': filename,
        'label': label,
        'fh': open(filename, 'r', encoding='utf-8', errors='replace'),
        'col_idx': col_idx,
        'line_no': 0,
        'ncols': None,
        'fields': None,
        'key': None,
    }
    advance(state)
    return state


def close_reader(state):
    state['fh'].close()


def advance(state):
    fields, key = read_next_row(state)
    state['fields'] = fields
    state['key'] = key


def pop_group(state):
    if state['fields'] is None:
        raise ValueError(f'{state["label"]}: internal error, no current row')

    key = state['key']
    rows = []

    while state['fields'] is not None and state['key'] == key:
        rows.append(state['fields'])
        advance(state)

    return key, rows


def write_rows(out, left_rows, right_rows, left_ncols, right_ncols):
    left_nan = [NAN] * left_ncols
    right_nan = [NAN] * right_ncols

    if left_rows and right_rows:
        for left in left_rows:
            for right in right_rows:
                out.write('\t'.join(left + right) + '\n')
    elif left_rows:
        for left in left_rows:
            out.write('\t'.join(left + right_nan) + '\n')
    else:
        for right in right_rows:
            out.write('\t'.join(left_nan + right) + '\n')

# def within_tolerance(last_ns1, last_ns2, current_ns1, current_ns2, max_drift=1000, jitter=100, clock_scale_factor=1.0):
#     elapsed1 = current_ns1 - last_ns1;
#     elapsed2 = current_ns2 - last_ns2;
#     difference = abs(elapsed1 - elapsed2 * clock_scale_factor)

import numpy as np

def load_scaled_column(path: str, column: int, scale: float) -> np.ndarray:
    """
    Read one whitespace-separated column from a text file,
    multiply it by scale, and return a NumPy float64 array.

    column is 1-based:
        column = 1  -> first column
        column = 2  -> second column
    """
    values = np.loadtxt(
        path,
        usecols=(column - 1,),   # convert 1-based to 0-based
        dtype=np.float64
    )
    values *= scale          # in-place scaling
    return values


def get_correlated_count(values1, values2, offset_from, offset_to):
    correlated = 0
    for value1 in values1:
        idx1 = values2.searchsorted(offset_from + value1, side="right")
        idx2 = values2.searchsorted(offset_to + value1, side="right")
        # print(f"DEBUG off_from={offset_from}, off_to={offset_to}, idx1={idx1},idx2={idx2}")
        if idx2 > idx1:
            correlated += 1
        # correlated += (idx2-idx1)
    return correlated


def recursive_search(values1, values2, offset_estimate, search_range_s, steps, stoprange_s, tolerance_ns,debug,overlap):
    if debug: print(f"#Debug recursive_search(values1, values2, offset_estimate={offset_estimate}, search_range_s={search_range_s}, steps, stoprange_s, tolerance_ns):")
    guess_offset = 0
    results = []  # list of (count, guess_offset)
    compute_offsets_values = [offset_estimate - search_range_s / 2 + i * search_range_s / (steps - 1) for i in range(steps)]
    step_s = overlap * search_range_s / (steps - 1)
    # print(f"#compute offset values={compute_offsets_values}, off={offset_estimate}, range={search_range_s}")
    for guess_offset in compute_offsets_values:

        offset_from_s = guess_offset - step_s - tolerance_ns * 0.000000001
        offset_to_s = guess_offset + step_s + tolerance_ns * 0.000000001
        correlated_rough = get_correlated_count(values1, values2, offset_from_s, offset_to_s)
        # print(f"#guess: {guess_offset:.6f} - {(guess_offset+search_range_s/(steps-1)):.6f}: correlated {correlated_rough}")
        if debug: print(f"{search_range_s:.9f} {(offset_to_s-offset_from_s)*1000000000:.0f} {offset_from_s:.9f} {offset_to_s:.9f} {correlated_rough}")
        results.append((correlated_rough, guess_offset))
    results.sort(key=lambda x: x[0],reverse=True)
    #results.sort(reverse=True)
    if debug: print(); print()

    for result in results:
        guess_offset = result[1]
        offset_from_s = guess_offset - step_s - tolerance_ns * 0.000000001
        offset_to_s = guess_offset + step_s + tolerance_ns * 0.000000001
        # offset_from_s = result[1] - tolerance_ns * 0.000000001
        # offset_to_s = result[1] + 2*search_range_s / (steps ) + tolerance_ns * 0.000000001
        # guess_offset = 0.5 * (offset_from_s + offset_to_s)
        if debug: print(f"#DEBUG guess_offset={guess_offset}")
        if search_range_s / (steps) > stoprange_s:
            guess_offset = recursive_search(values1, values2, guess_offset, 2 * search_range_s / (steps - 1), steps, stoprange_s, tolerance_ns, debug, overlap)
        else: 
            if debug: print(f"{guess_offset} #final offset")
        # TODO explore equal offset?
        break

    # print(f"#interm. result range={search_range_s}, steps={steps} => {results[0][0]} correlations in offset {results[0][1]} s")

    return guess_offset
# def is_better_path(list1,list2): #is list1 better than list2? List1 can never be shorter than list2
#     for i in range(len(list2)):
#         if list1[i][0]>list2[i][0]:
#             return True
#         elif list1[i][0]<list2[i][0]:
#             return False
#     return False
#
#
# def recursive_search(values1, values2, offset_estimate, search_range_s, steps, stoprange_s, tolerance_ns, debug, overlap, best_path):
#     if debug: print(f"#Debug recursive_search(values1, values2, offset_estimate={offset_estimate}, search_range_s={search_range_s}, steps, stoprange_s, tolerance_ns):")
#     guess_offset = 0
#     result_candidates = []  # list of (count, guess_offset)
#     compute_offsets_values = [offset_estimate - search_range_s / 2 + i * search_range_s / (steps - 1) for i in range(steps)]
#     step_s = overlap * search_range_s / (steps - 1)
#     # print(f"#compute offset values={compute_offsets_values}, off={offset_estimate}, range={search_range_s}")
#     for guess_offset in compute_offsets_values:
#         offset_from_s = guess_offset - step_s - tolerance_ns * 0.000000001
#         offset_to_s = guess_offset + step_s + tolerance_ns * 0.000000001
#         correlated_rough = get_correlated_count(values1, values2, offset_from_s, offset_to_s)
#         # print(f"#guess: {guess_offset:.6f} - {(guess_offset+search_range_s/(steps-1)):.6f}: correlated {correlated_rough}")
#         if debug: print(f"{search_range_s:.9f} {(offset_to_s-offset_from_s)*1000000000:.0f} {offset_from_s:.9f} {offset_to_s:.9f} {correlated_rough}")
#         result_candidates.append((correlated_rough, guess_offset))
#
#     result_candidates.sort(reverse=True)
#     if debug: print(); print()
#     new_path = best_path[1:]
#     best_count = result_candidates[0][0]
#     best_offset = result_candidates[0][1]
#     for result in result_candidates:
#         guess_count = result[0]
#         if guess_count < best_count:
#             break
#         guess_offset = result[1]
#         if len(best_path)>0:
#             if best_path[-1][0]>guess_count:
#                 print(f"better branch found. leaving this one")
#                 return [(guess_count,guess_offset)] #end of story, dead end
#         offset_from_s = guess_offset - step_s - tolerance_ns * 0.000000001
#         offset_to_s = guess_offset + step_s + tolerance_ns * 0.000000001
#         # offset_from_s = result[1] - tolerance_ns * 0.000000001
#         # offset_to_s = result[1] + 2*search_range_s / (steps ) + tolerance_ns * 0.000000001
#         # guess_offset = 0.5 * (offset_from_s + offset_to_s)
#         if debug: print(f"#DEBUG guess_offset={guess_offset}")
#         if search_range_s / (steps) > stoprange_s:
#             guess_path = recursive_search(values1, values2, guess_offset, 2 * search_range_s / (steps - 1),
#                                              steps, stoprange_s, tolerance_ns, debug, overlap, new_path)
#             if is_better_path(guess_path, new_path):
#                 print(f"new path is better {guess_path}")
#                 new_path = guess_path
#
#             else:
#                 print("not a better path")
#         else: 
#             if debug: print(f"{guess_offset} #final offset")
#             return [(guess_count,guess_offset)] 
#         # TODO explore equal offset?
#     return new_path
#
#     # print(f"#interm. result range={search_range_s}, steps={steps} => {result_candidates[0][0]} correlations in offset {result_candidates[0][1]} s")
#
#     return guess_offset


def main(argv):
    long_options = [
        'help',
        'file1=',
        'file2=',
        'col1=',
        'col2=',
        # 'offset=',
        'scale1=',
        'scale2=',
        'jitter=',
        'output=',
        'iteration-steps=',
        'iteration-overlap=',
        'iteration-debug',
        'print-skipped',
        'initial-range=',
        'initial-offset='
    ]
    try:
        opts, _ = getopt.getopt(argv, 'ho:', long_options)
    except getopt.GetoptError:
        usage(2)

    file1 = None
    file2 = None
    col1 = None
    col2 = None
    scale1 = 8.0
    scale2 = 1000000000.0
    # shift1 = None
    # shift2 = None
    jitter = 1000
    output = '-'
    iteration_overlap = 0.55 #min 0.5
    iteration_debug = 0  
    iteration_steps = 8
    initial_range = 500
    initial_offset = 0
    print_skipped = 0
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage(0)
        elif opt == '--file1':
            file1 = arg
        elif opt == '--file2':
            file2 = arg
        elif opt == '--col1':
            col1 = int(arg)
        elif opt == '--col2':
            col2 = int(arg)
        elif opt == '--scale1':
            scale1 = float(arg)
        elif opt == '--scale2':
            scale2 = float(arg)
        # elif opt == '--shift1':
        #     shift1 = float(arg)
        # elif opt == '--shift2':
        #     shift2 = float(arg)
        elif opt == '--jitter': 
            jitter = float(arg)
        elif opt in ('-o', '--output'):
            output = arg
        elif opt == '--iteration-steps':
            iteration_steps = int(arg)
        elif opt == '--iteration-debug':
            iteration_debug = 1
        elif opt == '--iteration-overlap':
            iteration_overlap = float(arg)
        elif opt == '--print-skipped':
            print_skipped = 1
        elif opt == '--initial-range':
            initial_range = float(arg)
        elif opt == '--initial-offset':
            initial_offset = float(arg)

    if file1 is None or file2 is None or col1 is None or col2 is None:
        usage(2)
    if col1 < 1 or col2 < 1:
        usage(2)

    left = None
    right = None
    out = None
    
    dropped_events1=0
    dropped_events2=0
    synchronized_events = 0

    try:
        values1 = load_scaled_column(file1, column=col1, scale=scale1) * 0.000000001
        values2 = load_scaled_column(file2, column=col2, scale=scale2) * 0.000000001
        offset = 0  # initial
        #it is necessary the array1 has 
        if len(values1) < len(values2):
            offset = recursive_search(values1=values1, values2=values2, offset_estimate=initial_offset,
                         search_range_s=initial_range, steps=iteration_steps, stoprange_s=0.000000001,
                         tolerance_ns=8, debug=iteration_debug, overlap=iteration_overlap)
        else:
            offset = -recursive_search(values1=values2, values2=values1, offset_estimate=-initial_offset,
                         search_range_s=initial_range, steps=iteration_steps, stoprange_s=0.000000001,
                         tolerance_ns=8, debug=iteration_debug, overlap=iteration_overlap)
        print(f"#offset {offset}")
            
        left = init_reader(file1, col1 - 1, 'file1')
        right = init_reader(file2, col2 - 1, 'file2')
        if left['ncols'] is None:
            print(f'event_time_merge.py: error: file1 has no data rows: "{file1}"', file=sys.stderr)
            return 2
        if right['ncols'] is None:
            print(f'event_time_merge.py: error: file2 has no data rows: "{file2}"', file=sys.stderr)
            return 2
        
        out = open_output(output)

        # first_ns1 = left['key'] * scale1
        # first_ns2 = right['key'] * scale2
         
        # if shift1 == None: shift1 = -left['key'] * scale1
        # if shift2 == None: shift2 = -left['key'] * scale2
        # estimated_scale_factor = clock_scale_factor
        
        # last_ns1 = left['key'] * scale1 - 1  # adding -1 for the first event to avoid division by 0
        # last_ns2 = right['key'] * scale2 - 1

        while left['fields'] is not None or right['fields'] is not None:
            # taking care of rest of the events in the file1
            if left['fields'] is None:
                _, right_rows = pop_group(right)
                dropped_events2 += 1
                if print_skipped: write_rows(out, [], right_rows, left['ncols'], right['ncols'])
                continue
            # taking care of rest of the events in the file2
            if right['fields'] is None:
                _, left_rows = pop_group(left)
                dropped_events1 += 1
                if print_skipped: write_rows(out, left_rows, [], left['ncols'], right['ncols'])
                continue
            
            # here we have both events filled in the left and right    
            ns1 = left['key'] * scale1 + offset*1000000000
            ns2 = right['key'] * scale2 
            
            # elapsed1 = ns1 - last_ns1;
            # elapsed2 = ns2 - last_ns2;
            # elapsed1 = ns1 - last_ns1;
            # elapsed2 = ns2 - last_ns2;
            difference = ns1 - ns2
            # print(f"difference={difference}")
            if abs(difference) > ( jitter):
                if difference < 0:
                    _, left_rows = pop_group(left)
                    dropped_events1 += 1
                    if print_skipped: write_rows(out, left_rows, [], left['ncols'], right['ncols'])
                else:
                    _, right_rows = pop_group(right)
                    dropped_events2 += 1
                    if print_skipped: write_rows(out, [], right_rows, left['ncols'], right['ncols'])
                continue
            else:
                
                synchronized_events += 1
                # last_ns1 = left['key'] * scale1 
                # last_ns2 = right['key'] * scale2
                # if last_ns2 != first_ns2: 
                #     clock_scale_factor = (first_ns1 - last_ns1) / (first_ns2 - last_ns2)
                # print(f'#Debug: new scale {clock_scale_factor}', file=sys.stderr)
                _, left_rows = pop_group(left)
                _, right_rows = pop_group(right)
                write_rows(out, left_rows, right_rows, left['ncols'], right['ncols'])
        print(f'#Summary offset={offset}\tdrop1={dropped_events1}\tdrop2={dropped_events2}\tsynchronized={synchronized_events}', file=sys.stderr)
        print(f'#Summary offset={offset}\tdrop1={dropped_events1}\tdrop2={dropped_events2}\tsynchronized={synchronized_events}')
    except ValueError as exc:
        print(f'event_time_merge.py: error: {exc}', file=sys.stderr)
        return 2
    finally:
        if left is not None:
            close_reader(left)
        if right is not None:
            close_reader(right)
        if out is not None and out is not sys.stdout:
            out.close()

    return 0 


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
