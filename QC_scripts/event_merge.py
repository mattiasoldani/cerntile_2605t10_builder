#!/usr/bin/env python3
import sys
import getopt
from pathlib import Path


NAN = 'NaN'


def usage(exit_code=1):
    prog = Path(sys.argv[0]).name
    print(
        f'Usage: python3 {prog} --file1 FILE --file2 FILE --col1 COL --col2 COL '
        '[--offset OFFSET] [-o OUTPUT]',
        file=sys.stderr,
    )
    print('', file=sys.stderr)
    print('  --file1           first tab-separated file (sorted by merge column)', file=sys.stderr)
    print('  --file2           second tab-separated file (sorted by merge column)', file=sys.stderr)
    print('  --col1            1-based merge column in file1', file=sys.stderr)
    print('  --col2            1-based merge column in file2', file=sys.stderr)
    print('  --offset          numeric offset added to file2 merge key          (default: 0)', file=sys.stderr)
    print('  -o, --output      output file                                     (default: stdout)', file=sys.stderr)
    print('  -h, --help        show this help', file=sys.stderr)
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
            key = float(key_text) + state['offset']
        except ValueError as exc:
            raise ValueError(
                f'{state["label"]}: cannot parse merge key "{key_text}" '
                f'in "{state["filename"]}", line {state["line_no"]}'
            ) from exc

        return fields, key

    return None, None


def init_reader(filename, col_idx, offset, label):
    state = {
        'filename': filename,
        'label': label,
        'fh': open(filename, 'r', encoding='utf-8', errors='replace'),
        'col_idx': col_idx,
        'offset': offset,
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


def main(argv):
    long_options = [
        'help',
        'file1=',
        'file2=',
        'col1=',
        'col2=',
        'offset=',
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
    offset = 0.0
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
        elif opt == '--offset':
            offset = float(arg)
        elif opt in ('-o', '--output'):
            output = arg

    if file1 is None or file2 is None or col1 is None or col2 is None:
        usage(2)
    if col1 < 1 or col2 < 1:
        usage(2)

    left = None
    right = None
    out = None

    try:
        left = init_reader(file1, col1 - 1, 0.0, 'file1')
        right = init_reader(file2, col2 - 1, offset, 'file2')

        if left['ncols'] is None:
            print(f'tsv-merge.py: error: file1 has no data rows: "{file1}"', file=sys.stderr)
            return 2
        if right['ncols'] is None:
            print(f'tsv-merge.py: error: file2 has no data rows: "{file2}"', file=sys.stderr)
            return 2

        out = open_output(output)

        while left['fields'] is not None or right['fields'] is not None:
            if left['fields'] is None:
                _, right_rows = pop_group(right)
                write_rows(out, [], right_rows, left['ncols'], right['ncols'])
                continue

            if right['fields'] is None:
                _, left_rows = pop_group(left)
                write_rows(out, left_rows, [], left['ncols'], right['ncols'])
                continue

            if left['key'] == right['key']:
                _, left_rows = pop_group(left)
                _, right_rows = pop_group(right)
                write_rows(out, left_rows, right_rows, left['ncols'], right['ncols'])
            elif left['key'] < right['key']:
                _, left_rows = pop_group(left)
                write_rows(out, left_rows, [], left['ncols'], right['ncols'])
            else:
                _, right_rows = pop_group(right)
                write_rows(out, [], right_rows, left['ncols'], right['ncols'])

    except ValueError as exc:
        print(f'tsv-merge.py: error: {exc}', file=sys.stderr)
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
