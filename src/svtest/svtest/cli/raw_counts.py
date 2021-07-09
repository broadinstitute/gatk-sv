#!/usr/bin/env python

"""
Collect raw read count metrics. Writes stats to stdout.

Metrics:
  counts_qN             : Nth percentile
  counts_num_intervals  : number of records
  counts_intervals_size : sum of interval sizes

"""

import gzip
import argparse
import sys
import numpy as np
import svtest.utils.TestUtils as tu

Q25_KEY = "rd_q25"
Q50_KEY = "rd_q50"
Q75_KEY = "rd_q75"
MEAN_KEY = "rd_mean"
NUM_ZERO = "rd_num_zero"
NUM_INTERVALS = "rd_num_intervals"
INTERVALS_SIZE = "rd_intervals_size"

EXPECTED_COLUMNS = 4
HEADER_CHAR = '@'


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtest raw-counts',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('counts_file', type=str)
    parser.add_argument('sample_id', type=str)

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    # Read file
    with gzip.open(args.counts_file, mode='rb') as f:
        metrics = get_metrics(f, args.sample_id)

    # Write metrics
    write_metrics(metrics)


def write_metrics(metrics):
    for key in metrics:
        sys.stdout.write("%s\t%s\n" % (key, str(metrics[key])))


def get_metrics(file, sample_id):
    counts = []
    intervals_size = 0
    last_line_was_header = False
    for bin in file:
        line = bin.decode()
        if line.startswith(HEADER_CHAR):
            last_line_was_header = True
            continue
        if last_line_was_header:
            last_line_was_header = False
            continue  # skip columns header line
        tokens = line.strip().split('\t')
        tu.test_iterable_size(tokens, EXPECTED_COLUMNS)
        start = int(tokens[1])
        end = int(tokens[2])
        intervals_size += end - start + 1
        count = int(tokens[3])
        counts.append(int(count))
    counts_arr = np.asarray(counts)
    quantiles = np.quantile(counts_arr, [0.25, 0.50, 0.75])
    mean = np.mean(counts_arr)
    num_zero = counts_arr[counts_arr == 0].size
    return {
        Q25_KEY + "_" + sample_id: quantiles[0],
        Q50_KEY + "_" + sample_id: quantiles[1],
        Q75_KEY + "_" + sample_id: quantiles[2],
        MEAN_KEY + "_" + sample_id: mean,
        NUM_ZERO + "_" + sample_id: num_zero,
        NUM_INTERVALS + "_" + sample_id: len(counts_arr),
        INTERVALS_SIZE + "_" + sample_id: intervals_size
    }


if __name__ == '__main__':
    main()
