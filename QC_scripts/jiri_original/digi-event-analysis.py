#!/usr/bin/env python3
import sys
import getopt
from pathlib import Path

import numpy as np

try:
    from scipy.signal import butter, lfilter, filtfilt
except Exception:
    butter = None
    lfilter = None
    filtfilt = None


NAN = float('nan')


def usage(exit_code=1):
    prog = Path(sys.argv[0]).name
    print(
        f'Usage: cat wave*.txt | python3 {prog} '
        '[-b BASELINE] [-m METHOD] [-o ORDER] [-c CUT] [-w WINDOW] '
        '[--peak-from FROM] [--peak-to TO] [--timestamps-only]',
        file=sys.stderr,
    )
    print('', file=sys.stderr)
    print('  -b, --baseline         baseline file from baseline.py          (default: none -> 0)', file=sys.stderr)
    print('  -m, --method           none | butter | filtfilt | ma | binom   (default: none)', file=sys.stderr)
    print('  -o, --order            butterworth order                       (default: 2)', file=sys.stderr)
    print('  -c, --cut              lowpass cutoff, fraction of Nyquist     (default: 0.12)', file=sys.stderr)
    print('  -w, --window           FIR window length for ma/binom          (default: 11)', file=sys.stderr)
    print('      --peak-from        first sample index for peak search      (default: 0)', file=sys.stderr)
    print('      --peak-to          last sample index for peak search       (default: record length)', file=sys.stderr)
    print('      --baseline-from    evnt-by-event baseline sample position  (default:-1)', file=sys.stderr)
    print('      --baseline-to      evnt-by-event baseline sample position  (default:-1)', file=sys.stderr)
    print('      --timestamps-only  do not parse samples, fill analysis with NaN', file=sys.stderr)
    print('      --max-readout-delay  physical time from trigger to readout (default: 7)', file=sys.stderr)
    print('  -h, --help             show this help', file=sys.stderr)
    raise SystemExit(exit_code)


def read_baselines(filename):
    baselines = {}
    if not filename:
        return baselines

    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) < 2:
                continue
            baselines[int(fields[0])] = float(fields[1])
    return baselines


def read_value(lines, key):
    line = next(lines).rstrip('\n')
    name, value = line.split(':', 1)
    if name != key:
        raise ValueError(f'expected "{key}", got "{name}"')
    return value.strip()


def read_event(lines, with_samples):
    for line in lines:
        line = line.rstrip('\n')
        if line:
            break
    else:
        return None

    if not line.startswith('Record Length:'):
        raise ValueError(f'expected "Record Length:", got "{line}"')

    record_length = int(line.split(':', 1)[1])
    board_id = read_value(lines, 'BoardID')
    channel = int(read_value(lines, 'Channel'))
    event_number = int(read_value(lines, 'Event Number'))
    pattern = read_value(lines, 'Pattern')
    trigger_timestamp = int(read_value(lines, 'Trigger Time Stamp'))
    dc_offset_dac = read_value(lines, 'DC offset (DAC)')
    ro_time_s = int(read_value(lines, 'RO Time (significand in s)'))
    ro_time_ns = int(read_value(lines, 'RO Time (decimals in ns)'))

    samples = None
    if with_samples:
        sample_text = ''.join(next(lines) for _ in range(record_length))
        samples = np.fromstring(sample_text, sep='\n', dtype=np.int32)
        if samples.size != record_length:
            raise ValueError(
                f'expected {record_length} samples in channel {channel}, got {samples.size}'
            )
    else:
        for _ in range(record_length):
            next(lines)

    return {
        'record_length': record_length,
        'board_id': board_id,
        'channel': channel,
        'event_number': event_number,
        'pattern': pattern,
        'trigger_timestamp': trigger_timestamp,
        'dc_offset_dac': dc_offset_dac,
        'ro_time_s': ro_time_s,
        'ro_time_ns': ro_time_ns,
        'samples': samples,
    }

def individual_baseline(samples,baseline_from,baseline_to):
    if baseline_to < baseline_from:
        return 0.0
    return float(np.mean(samples[baseline_from:baseline_to]))
    

def fir_filter(samples, kernel):
    left = (len(kernel) - 1) // 2
    right = len(kernel) - 1 - left
    padded = np.pad(samples, (left, right), mode='edge')
    return np.convolve(padded, kernel, mode='valid')


def moving_average(samples, window):
    kernel = np.ones(window, dtype=np.float64) / window
    return fir_filter(samples, kernel)


def binomial_filter(samples, window):
    kernel = np.array([1.0], dtype=np.float64)
    for _ in range(window - 1):
        kernel = np.convolve(kernel, np.array([1.0, 1.0], dtype=np.float64))
    kernel /= kernel.sum()
    return fir_filter(samples, kernel)


def butter_filter(samples, zero_phase, coeffs):
    if butter is None or lfilter is None or filtfilt is None:
        raise ValueError('butter/filtfilt needs scipy')
    b, a = coeffs
    if zero_phase:
        return filtfilt(b, a, samples)
    return lfilter(b, a, samples)


def apply_filter(samples, method, window, coeffs):
    if method == 'none':
        return samples
    if method == 'ma':
        return moving_average(samples, window)
    if method == 'binom':
        return binomial_filter(samples, window)
    if method == 'butter':
        return butter_filter(samples, False, coeffs)
    if method == 'filtfilt':
        return butter_filter(samples, True, coeffs)
    raise ValueError(f'unknown method "{method}"')


def interpolate_zero(x0, y0, x1, y1):
    if y1 == y0:
        return float(x0)
    return float(x0) - float(y0) * (float(x1) - float(x0)) / (float(y1) - float(y0))


def find_left_crossing(samples, peak_pos):
    if samples[peak_pos] <= 0:
        return NAN

    for i in range(peak_pos, 0, -1):
        y0 = samples[i - 1]
        y1 = samples[i]
        if y1 == 0:
            return float(i)
        if y0 <= 0 < y1:
            return interpolate_zero(i - 1, y0, i, y1)

    if samples[0] == 0:
        return 0.0
    return NAN


def find_right_crossing(samples, peak_pos):
    if samples[peak_pos] <= 0:
        return NAN

    n = len(samples)
    for i in range(peak_pos, n - 1):
        y0 = samples[i]
        y1 = samples[i + 1]
        if y0 == 0:
            return float(i)
        if y0 > 0 >= y1:
            return interpolate_zero(i, y0, i + 1, y1)

    if samples[-1] == 0:
        return float(n - 1)
    return NAN


def analyze_positive_pulse(samples, peak_from, peak_to):
    n = len(samples)
    lo = 0 if peak_from is None else peak_from
    hi = n if peak_to is None else peak_to

    if lo < 0 or hi <= lo or hi > n:
        raise ValueError(f'peak window [{lo}, {hi}) exceeds waveform length {n}')

    peak_offset = int(np.argmax(samples[lo:hi]))
    peak_pos = lo + peak_offset
    peak_value = float(samples[peak_pos])
    left_cross = find_left_crossing(samples, peak_pos)
    right_cross = find_right_crossing(samples, peak_pos)

    return {
        'peak_pos': float(peak_pos),
        'peak_value': peak_value,
        'left_cross': left_cross,
        'right_cross': right_cross,
    }


def fmt(value):
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, str):
        return value
    if value is None:
        return 'nan'
    value = float(value)
    if np.isnan(value):
        return 'nan'
    return f'{value:.10g}'


def write_header():
    columns = [
        'board_id',
        'channel',
        'event_number',
        'event_increment',
        'record_length',
        'trigger_timestamp',
        'timestamp_difference',
        'ro_time_s',
        'ro_time_ns',
        'ro_time',
        'peak_pos',
        'peak_value',
        'left_cross',
        'right_cross',
    ]
    sys.stdout.write('#' + '\t'.join(columns) + '\n')


def write_row(event, metrics):
    ro_time = float(event['ro_time_s']) + 1e-9 * float(event['ro_time_ns'])
    values = [
        event['board_id'],
        event['channel'],
        event['event_number'],
        event['event_increment'],
        event['record_length'],
        event['trigger_timestamp'],
        event['timestamp_difference'],
        event['ro_time_s'],
        event['ro_time_ns'],
        ro_time,
        metrics['peak_pos'],
        metrics['peak_value'],
        metrics['left_cross'],
        metrics['right_cross'],
    ]
    sys.stdout.write('\t'.join(fmt(v) for v in values) + '\n')


def main(argv):
    long_options = [
        'help',
        'baseline=',
        'method=',
        'order=',
        'cut=',
        'window=',
        'peak-from=',
        'peak-to=',
        'baseline-from=',
        'baseline-to=',
        'timestamps-only',
        'max-readout-delay='
    ]
    try:
        opts, _ = getopt.getopt(argv, 'hb:m:o:c:w:', long_options)
    except getopt.GetoptError:
        usage(2)

    baseline_file = None
    method = 'none'
    order = 2
    cut = 0.12
    window = 11
    peak_from = None
    peak_to = None
    baseline_from = None
    baseline_to = None
    timestamps_only = False
    max_readout_delay = 7.0  # more than 17 is bad, less than 2 s is also bad.

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage(0)
        elif opt in ('-b', '--baseline'):
            baseline_file = arg
        elif opt in ('-m', '--method'):
            method = arg
        elif opt in ('-o', '--order'):
            order = int(arg)
        elif opt in ('-c', '--cut'):
            cut = float(arg)
        elif opt in ('-w', '--window'):
            window = int(arg)
        elif opt == '--peak-from':
            peak_from = int(arg)
        elif opt == '--peak-to':
            peak_to = int(arg)
        elif opt == '--baseline-from':
            baseline_from = int(arg)
        elif opt == '--baseline-to':
            baseline_to = int(arg)
        elif opt == '--timestamps-only':
            timestamps_only = True
        elif opt == '--max-readout-delay':
            max_readout_delay = float(arg)

    if method not in ('none', 'butter', 'filtfilt', 'ma', 'binom'):
        print(f'event-analysis.py: error: unknown method "{method}"', file=sys.stderr)
        return 2
    if order < 1 or not (0.0 < cut < 1.0) or window < 1:
        usage(2)
    if peak_from is not None and peak_from < 0:
        usage(2)
    if peak_to is not None and peak_to < 0:
        usage(2)
    if peak_from is not None and peak_to is not None and peak_to <= peak_from:
        usage(2)

    if baseline_file is not None:
        baselines = read_baselines(baseline_file)

    coeffs = None
    if method in ('butter', 'filtfilt'):
        if butter is None or lfilter is None or filtfilt is None:
            raise ValueError('butter/filtfilt needs scipy')
        coeffs = butter(order, cut, btype='low')

    write_header() 
    print(f"#parameters: --baseline-from {baseline_from} --baseline-to {baseline_to} --peak-from {peak_from} --peak-to {peak_to} --method {method}")
    previous_trigger_ts_by_channel = {}
    previous_event_number_by_channel = {}
    first_uts_by_channel = {} #unix timestamp in seconds, float
    first_ts_by_channel = {} #clock counter
    trigger_overflow_addition = {i:0 for i in range(8)}
    last_readout_time = None
    try:
        lines = iter(sys.stdin)
        while True:
            event = read_event(lines, with_samples=not timestamps_only)
            if event is None:
                break

            #block to calculate increments from previous events
            channel = event['channel']
            event_number = event['event_number']
            ro_utime = event['ro_time_s'] + 0.000000001 * event['ro_time_ns']
            if not (channel in first_ts_by_channel):
                first_uts_by_channel[channel] = ro_utime 
                first_ts_by_channel[channel] = event['trigger_timestamp']

            event['trigger_timestamp'] += trigger_overflow_addition[channel]
            trigger_timestamp = event['trigger_timestamp']
            while (0.000000008 * (event['trigger_timestamp'] -first_ts_by_channel[channel]) + max_readout_delay) < (ro_utime-first_uts_by_channel[channel]):
                #print(f' trigger_timestamp={trigger_timestamp}  ,ro_utime-first_uts_by_channel[channel]={ro_utime-first_uts_by_channel[channel]}')
                print(f'#time overflow detected, moving 17.16 s forward for event {event_number}')
                trigger_overflow_addition[channel] += 2147483648
                event['trigger_timestamp'] += 2147483648
                trigger_timestamp += 2147483648
            # print(f' trigger_timestamp={trigger_timestamp}  ,ro_utime-first_uts_by_channel[channel]={ro_utime-first_uts_by_channel[channel]}')
            if channel in previous_trigger_ts_by_channel:
                event['timestamp_difference'] = trigger_timestamp - previous_trigger_ts_by_channel[channel]
            else:
                event['timestamp_difference'] = NAN
            previous_trigger_ts_by_channel[channel] = trigger_timestamp
            # if channel in previous_trigger_ts_by_channel:
            #     difference = trigger_timestamp - previous_trigger_ts_by_channel[channel]
            #     if difference < 0 :
            #         trigger_overflow_addition[channel] += 2147483648
            #         event['trigger_timestamp'] += 2147483648
            #         trigger_timestamp += 2147483648
            #         difference += 2147483648
            #     event['timestamp_difference'] = trigger_timestamp - previous_trigger_ts_by_channel[channel]
            # else:
            #     event['timestamp_difference'] = NAN
            # previous_trigger_ts_by_channel[channel] = trigger_timestamp
            

            if channel in previous_event_number_by_channel:
                event['event_increment'] = event_number - previous_event_number_by_channel[channel]
            else:
                event['event_increment'] = NAN
            previous_event_number_by_channel[channel] = event_number

            if timestamps_only:
                metrics = {
                    'peak_pos': NAN,
                    'peak_value': NAN,
                    'left_cross': NAN,
                    'right_cross': NAN,
                }
            else:
                if (baseline_from is not None) and (baseline_to is not None):
                    baseline = individual_baseline(event['samples'], baseline_from, baseline_to)
                else: 
                    baseline = baselines.get(event['channel'], 0.0)
                                            
                centered = event['samples'].astype(np.float64) - baseline
                filtered = apply_filter(centered, method, window, coeffs)
                metrics = analyze_positive_pulse(filtered, peak_from, peak_to)

            write_row(event, metrics)

    except StopIteration:
        print('event-analysis.py: error: unexpected end of input', file=sys.stderr)
        return 2
    except ValueError as exc:
        print(f'event-analysis.py: error: {exc}', file=sys.stderr)
        return 2

    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
