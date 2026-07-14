#!/usr/bin/env python3
import sys
import getopt
from pathlib import Path
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
    print('  --max-drift    maximum allowed clock drift (ns per minute; default=1000)', file=sys.stderr)
    print('  --jitter       maximum accepted time jitter (ns; default: 100)', file=sys.stderr)
    print('  --event-offset number of events to drop from file2 before start.',file=sys.stderr) 
    print('                 positive: drop from file2, negative=drop from file1.', file=sys.stderr)
    print('  --drop         amount of events to drop from the beginning before', file=sys.stderr)
    print('                 synchronization even starts', file=sys.stderr)
    print('  --clock-scale-factor    factor for clock2 to compensate the internal clock references', file=sys.stderr)
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
        'max-drift=',
        'jitter=',
        'drop=',
        'event-offset=',
        'clock-scale-factor=',
        'output=',
    ]
    try:
        opts, _ = getopt.getopt(argv, 'ho:', long_options)
    except getopt.GetoptError:
        usage(2)

    file1 = None
    file2 = None
    col1 = None
    col2 = None
    event_offset = 0
    scale1 = 8.0
    scale2 = 1.0
    # shift1 = None
    # shift2 = None
    drop = 0
    jitter = 100
    max_drift = 1000
    clock_scale_factor = 1.0
    output = '-'

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
        elif opt == '--drop':
            drop = int(arg)
        elif opt == '--event-offset':
            event_offset = int(arg)
        # elif opt == '--shift1':
        #     shift1 = float(arg)
        # elif opt == '--shift2':
        #     shift2 = float(arg)
        elif opt == '--jitter': 
            jitter = float(arg)
        elif opt == '--max-drift':
            event_offset = int(arg)
        elif opt == '--clock-scale-factor': 
            clock_scale_factor = float(arg)
        elif opt in ('-o', '--output'):
            output = arg

    if file1 is None or file2 is None or col1 is None or col2 is None:
        usage(2)
    if col1 < 1 or col2 < 1:
        usage(2)

    left = None
    right = None
    out = None
    
    dropped_events1=0
    dropped_events2=0
    synchronized_events=0

    try:
        left = init_reader(file1, col1 - 1,  'file1')
        right = init_reader(file2, col2 - 1, 'file2')
        
        if left['ncols'] is None:
            print(f'event_time_merge.py: error: file1 has no data rows: "{file1}"', file=sys.stderr)
            return 2
        if right['ncols'] is None:
            print(f'event_time_merge.py: error: file2 has no data rows: "{file2}"', file=sys.stderr)
            return 2

        out = open_output(output)

        for _ in range(event_offset): pop_group(right); dropped_events2 += 1
        for _ in range(-event_offset): pop_group(left); dropped_events1 += 1
        for _ in range(drop): pop_group(right);  pop_group(left); dropped_events2 += 1; dropped_events1 += 1
        
        first_ns1 = left['key'] * scale1
        first_ns2 = right['key'] * scale2
         
        # if shift1 == None: shift1 = -left['key'] * scale1
        # if shift2 == None: shift2 = -left['key'] * scale2
        # estimated_scale_factor = clock_scale_factor
        
        last_ns1 = left['key'] * scale1 - 1  # adding -1 for the first event to avoid division by 0
        last_ns2 = right['key'] * scale2 - 1

        while left['fields'] is not None or right['fields'] is not None:
            # taking care of rest of the events in the file1
            if left['fields'] is None:
                _, right_rows = pop_group(right)
                dropped_events2 += 1
                write_rows(out, [], right_rows, left['ncols'], right['ncols'])
                continue
            # taking care of rest of the events in the file2
            if right['fields'] is None:
                _, left_rows = pop_group(left)
                dropped_events1 += 1
                write_rows(out, left_rows, [], left['ncols'], right['ncols'])
                continue
            
            # here we have both events filled in the left and right    
            ns1 = left['key'] * scale1  
            ns2 = right['key'] * scale2 
            
            elapsed1 = ns1 - last_ns1;
            elapsed2 = ns2 - last_ns2;
            elapsed1 = ns1 - last_ns1;
            elapsed2 = ns2 - last_ns2;
            difference = elapsed1 - elapsed2 * clock_scale_factor
            if abs(difference) > (elapsed1 * .000000001 * max_drift + jitter):
                if difference < 0:
                    _, left_rows = pop_group(left)
                    dropped_events1 += 1
                    #write_rows(out, left_rows, [], left['ncols'], right['ncols'])
                else:
                    _, right_rows = pop_group(right)
                    dropped_events2 += 1
                    #write_rows(out, [], right_rows, left['ncols'], right['ncols'])
                continue
            else:
                synchronized_events += 1
                last_ns1 = left['key'] * scale1 
                last_ns2 = right['key'] * scale2
                if last_ns2 != first_ns2: 
                    clock_scale_factor = (first_ns1 - last_ns1) / (first_ns2 - last_ns2)
                # print(f'#Debug: new scale {clock_scale_factor}', file=sys.stderr)
                _, left_rows = pop_group(left)
                _, right_rows = pop_group(right)
                write_rows(out, left_rows, right_rows, left['ncols'], right['ncols'])
        print(f'#Summary evt-offset={event_offset}\tdrop1= {dropped_events1}\tdrop2= {dropped_events2}\ttotal= {synchronized_events}\tclock_scale= {clock_scale_factor}', file=sys.stderr)
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
