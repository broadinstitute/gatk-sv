#!/usr/bin/env python

"""
Collect SR file metrics. Writes stats to stdout.

Metrics:
  sr_left_<sample>  : number of left records
  sr_right_<sample> : number of right records

If the provided sample list has more than one id, <sample> will be "_merged".

"""

import gzip
import argparse
import sys
import svtest.utils.TestUtils as tu
import svtest.utils.IOUtils as iou

LEFT_KEY = "sr_left"
RIGHT_KEY = "sr_right"


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtest sr-file',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sr_file', type=str)
    parser.add_argument('sample_list', type=str)

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    # Read file
    with gzip.open(args.sr_file, mode='rb') as fsr:
        metrics = get_metrics(fsr, args.sample_list)

    # Write metrics
    write_metrics(metrics)


def write_metrics(metrics):
    for key in metrics:
        sys.stdout.write("%s\t%s\n" % (key, str(metrics[key])))


def get_metrics(sr_file, sample_list):
    samples = iou.read_samples_list(sample_list)
    samples_set = set(samples)
    side_metrics = [0, 0]
    for line in sr_file:
        tokens = line.decode().strip().split('\t')
        test_record(tokens, samples_set)
        side = tokens[2]
        if side == 'left':
            side_metrics[0] += 1
        elif side == 'right':
            side_metrics[1] += 1
        else:
            raise ValueError("Unrecognized orientation: %s" % side)

    if len(samples) == 1:
        metric_suffix = "_" + samples[0]
    else:
        metric_suffix = "_merged"

    return {
        LEFT_KEY + metric_suffix: side_metrics[0],
        RIGHT_KEY + metric_suffix: side_metrics[1]
    }


def test_record(columns, sample_ids):
    tu.test_iterable_size(columns, 5)
    tu.test_is_int(columns, 1)
    valid_strands = set(['right', 'left'])
    tu.test_column_in_iterable(columns, 2, valid_strands)
    tu.test_is_int(columns, 3)
    tu.test_column_in_iterable(columns, 4, sample_ids)


if __name__ == '__main__':
    main()
