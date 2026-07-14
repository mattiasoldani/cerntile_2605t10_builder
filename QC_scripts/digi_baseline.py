#!/usr/bin/env python3
import sys
import getopt
import struct
from pathlib import Path

def usage(exit_code=1):
    prog = Path(sys.argv[0]).name
    print(f'Usage: cat wave*.txt | python3 {prog} -f FROM -t TO [-c CH]...', file=sys.stderr)
    print('', file=sys.stderr)
    print('  -f, --from   first sample index (inclusive)', file=sys.stderr)
    print('  -t, --to     last sample index (exclusive)', file=sys.stderr)
    print('  -c, --ch     channel to include; may be used multiple times', file=sys.stderr)
    print('  -h, --help   show this help', file=sys.stderr)
    raise SystemExit(exit_code)


def float32(x):
    return struct.unpack('f', struct.pack('f', float(x)))[0]


def read_value(lines, key):
    line = next(lines).strip()
    name, value = line.split(':', 1)
    if name != key:
        raise ValueError(f'expected "{key}", got "{name}"')
    return value.strip()


def read_event_baseline(lines, sample_from, sample_to):
    for line in lines:
        line = line.strip()
        if line:
            break
    else:
        return None

    if not line.startswith('Record Length:'):
        raise ValueError(f'expected "Record Length", got "{line}"')
    record_length = int(line.split(':', 1)[1])

    read_value(lines, 'BoardID')
    channel = int(read_value(lines, 'Channel'))
    read_value(lines, 'Event Number')
    read_value(lines, 'Pattern')
    read_value(lines, 'Trigger Time Stamp')
    read_value(lines, 'DC offset (DAC)')
    read_value(lines, 'RO Time (significand in s)')
    read_value(lines, 'RO Time (decimals in ns)')

    if sample_to > record_length:
        raise ValueError(
            f'window [{sample_from}, {sample_to}) exceeds record length {record_length} in channel {channel}'
        )

    for _ in range(sample_from):
        next(lines)

    selected_sum = 0
    selected_count = sample_to - sample_from
    for _ in range(selected_count):
        selected_sum += int(next(lines))

    for _ in range(sample_to, record_length):
        next(lines)

    return channel, selected_sum, selected_count


def main(argv):
    long_options = ['help', 'from=', 'to=', 'ch=']
    try:
        opts, _ = getopt.getopt(argv, 'hf:t:c:', long_options)
    except getopt.GetoptError:
        usage(2)

    sample_from = None
    sample_to = None
    wanted_channels = set()

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage(0)
        elif opt in ('-f', '--from'):
            sample_from = int(arg)
        elif opt in ('-t', '--to'):
            sample_to = int(arg)
        elif opt in ('-c', '--ch'):
            wanted_channels.add(int(arg))

    if sample_from is None or sample_to is None or sample_from < 0 or sample_to <= sample_from:
        usage(2)

    sum_by_channel = {}
    count_by_channel = {}

    try:
        lines = iter(sys.stdin)
        while True:
            event = read_event_baseline(lines, sample_from, sample_to)
            if event is None:
                break

            channel, selected_sum, selected_count = event
            if wanted_channels and channel not in wanted_channels:
                continue

            sum_by_channel[channel] = sum_by_channel.get(channel, 0) + selected_sum
            count_by_channel[channel] = count_by_channel.get(channel, 0) + selected_count

    except StopIteration:
        print('baseline.py: error: unexpected end of input', file=sys.stderr)
        return 2
    except ValueError as exc:
        print(f'baseline.py: error: {exc}', file=sys.stderr)
        return 2

    if not count_by_channel:
        print('baseline.py: error: no matching events found', file=sys.stderr)
        return 1

    print('#ch\tbaseline')
    for channel in sorted(count_by_channel):
        baseline = sum_by_channel[channel] / count_by_channel[channel]
        print(f'{channel}\t{float32(baseline)}')

    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
