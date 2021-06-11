#!/usr/bin/env python

"""
Collect bincov matrix metrics. Writes stats to stdout.

Metrics:
  bincov_matrix_intervals : number of intervals (rows)
  bincov_matrix_intervals_all_samples_zero : number of intervals that were 0 across all samples
  bincov_matrix_intervals_at_least_one_zero : number of intervals with 0 in at least one sample
  bincov_matrix_q25 : 25th percentile count
  bincov_matrix_q50: 50th percentile count
  bincov_matrix_q75: 75th percentile count
  bincov_matrix_mean_<SAMPLE>: mean count per sample

"""

import argparse
import gzip
import sys
import numpy as np
import svtest.utils.TestUtils as tu
import svtest.utils.IOUtils as iou

INTERVALS_KEY = "bincov_matrix_intervals"
ALL_ZERO_KEY = "bincov_matrix_intervals_all_samples_zero"
ONE_ZERO_KEY = "bincov_matrix_intervals_at_least_one_zero"
Q25_KEY = "bincov_matrix_q25"
Q50_KEY = "bincov_matrix_q50"
Q75_KEY = "bincov_matrix_q75"
SAMPLE_MEAN_KEY = "bincov_matrix_mean"


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtest bincov-matrix',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bincov_matrix', type=str)
    parser.add_argument('sample_list', type=str)
    parser.add_argument('--low-mem-mode', action='store_true',
                        help='Only validate and calculate number of intervals')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    # Read file
    with gzip.open(args.bincov_matrix, mode='rb') as f:
        metrics = get_metrics(f, args.sample_list, args.low_mem_mode)

    # Write metrics
    write_metrics(metrics)


def write_metrics(metrics):
    for key in metrics:
        sys.stdout.write("%s\t%s\n" % (key, str(metrics[key])))


def get_metrics(matrix_file, sample_list, low_mem_mode):
    samples = iou.read_samples_list(sample_list)
    samples_set = set(samples)

    header = matrix_file.readline().decode().strip().split('\t')
    header_samples_set = set(header[3:])
    tu.test_sets_equal(header_samples_set, samples_set,
                       item_str="sample", name_a="header", name_b="samples list")

    data = []
    interval_size = None
    num_records = 0
    for line in matrix_file:
        num_records += 1
        tokens = line.decode().strip().split('\t')
        tu.test_is_int(tokens, 1)
        tu.test_is_int(tokens, 2)
        if interval_size is None:
            interval_size = int(tokens[2]) - int(tokens[1])
        else:
            if interval_size != int(tokens[2]) - int(tokens[1]):
                raise ValueError("Interval not of size {:d}: {:s}:{:d}-{:d}".format(interval_size, tokens[0],
                                                                                    int(tokens[1]), int(tokens[2])))
        counts = tokens[3:]
        test_record(counts, len(samples_set))
        if not low_mem_mode:
            data.append([int(x) for x in counts])

    if not low_mem_mode:
        arr = np.asarray(data)
        quantiles = np.quantile(arr, [0.25, 0.50, 0.75])
        max_over_samples = arr.max(axis=1)
        num_zero_in_all = len(max_over_samples[max_over_samples == 0])
        min_over_samples = arr.min(axis=1)
        num_zero_in_one = len(min_over_samples[min_over_samples == 0])
        metrics = {
            Q25_KEY: quantiles[0],
            Q50_KEY: quantiles[1],
            Q75_KEY: quantiles[2],
            INTERVALS_KEY: num_records,
            ALL_ZERO_KEY: num_zero_in_all,
            ONE_ZERO_KEY: num_zero_in_one
        }
        column_means = arr.mean(axis=0)
        col = 0
        for sample in header[3:]:
            metrics[SAMPLE_MEAN_KEY + "_" + sample] = column_means[col]
            col += 1
    else:
        metrics = {
            INTERVALS_KEY: num_records
        }

    return metrics


def test_record(columns, n_header_cols):
    tu.test_iterable_size(columns, n_header_cols)
    for i in range(1, len(columns)):
        tu.test_is_int(columns, i)


if __name__ == '__main__':
    main()
