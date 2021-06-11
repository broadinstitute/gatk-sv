#!/usr/bin/env python

"""
Collect PE file metrics. Writes stats to stdout.

Metrics:
  pe_LL_<sample>     : number of records with LL orientation
  pe_RR_<sample>     : number of records with RR orientation
  pe_RL_<sample>     : number of records with RL orientation
  pe_LR_<sample>     : number of records with LR orientation

If the provided sample list has more than one id, <sample> will be "_merged".

"""

import gzip
import argparse
import sys
import svtest.utils.TestUtils as tu
import svtest.utils.IOUtils as iou

PLUS_PLUS_KEY = "pe_LL"
MINUS_MINUS_KEY = "pe_RR"
PLUS_MINUS_KEY = "pe_LR"
MINUS_PLUS_KEY = "pe_RL"


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtest pe-file',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('pe_file', type=str)
    parser.add_argument('sample_list', type=str)

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    # Read file
    with gzip.open(args.pe_file, mode='rb') as fpe:
        metrics = get_metrics(fpe, args.sample_list)

    # Write metrics
    write_metrics(metrics)


def write_metrics(metrics):
    for key in metrics:
        sys.stdout.write("%s\t%s\n" % (key, str(metrics[key])))


def get_metrics(file, sample_list):
    samples = iou.read_samples_list(sample_list)
    samples_set = set(samples)

    data = [0, 0, 0, 0]  # ++, --, +-, -+
    for line in file:
        tokens = line.decode().strip().split('\t')
        test_record(tokens, samples_set)
        first = tokens[2]
        second = tokens[5]
        val = first + second
        if val == '++':
            data[0] += 1
        elif val == '--':
            data[1] += 1
        elif val == '+-':
            data[2] += 1
        elif val == '-+':
            data[3] += 1
        else:
            raise ValueError("Unrecognized orientation: %s / %s" %
                             (first, second))

    if len(samples) == 1:
        metric_suffix = "_" + samples[0]
    else:
        metric_suffix = "_merged"

    return {
        PLUS_PLUS_KEY + metric_suffix: data[0],
        MINUS_MINUS_KEY + metric_suffix: data[1],
        PLUS_MINUS_KEY + metric_suffix: data[2],
        MINUS_PLUS_KEY + metric_suffix: data[3]
    }


def test_record(columns, samples):
    tu.test_iterable_size(columns, 7)
    tu.test_is_int(columns, 1)
    valid_strands = set(['+', '-'])
    tu.test_column_in_iterable(columns, 2, valid_strands)
    tu.test_is_int(columns, 4)
    tu.test_column_in_iterable(columns, 5, valid_strands)
    tu.test_column_in_iterable(columns, 6, samples)


if __name__ == '__main__':
    main()
