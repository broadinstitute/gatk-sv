#!/usr/bin/env python

"""
Collect medcov metrics. Writes stats to stdout.

Metrics:
  medcov_mean     : mean median coverage
  medcov_mean_err : mean absolute deviation from baseline profile (if --baseline specified)
  medcov_max_err  : max absolute deviation from baseline profile (if --baseline specified)

"""

import argparse
import sys
import svtest.utils.TestUtils as tu
import svtest.utils.IOUtils as iou

MEAN_KEY = "medcov_mean"
MEAN_ERROR_KEY = "medcov_mean_abs_err"
MAX_ERROR_KEY = "medcov_max_abs_err"


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtest medcov',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('test_file', type=str)
    parser.add_argument('sample_list', type=str)
    parser.add_argument('--baseline-file', type=str, default=None)

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    samples = iou.read_samples_list(args.sample_list)

    # Read file
    with open(args.test_file, mode='r') as ftest:
        if args.baseline_file is None:
            metrics = get_metrics(ftest, None, samples)
        else:
            with open(args.baseline_file, mode='r') as fbase:
                metrics = get_metrics(ftest, fbase, samples)

    # Write metrics
    write_metrics(metrics)


def write_metrics(metrics):
    for key in metrics:
        sys.stdout.write("%s\t%s\n" % (key, str(metrics[key])))


def get_metrics(test_file, baseline_file, samples):
    test_header, test_data = get_medcov_file_data(test_file)
    tu.test_iterable_sizes_equal(
        test_header, test_data, name_a="test file header", name_b="test file data row")
    tu.test_sets_equal(test_header, samples, item_str="sample id",
                       name_a="test file header", name_b="sample list")
    metrics = {
        MEAN_KEY: float(sum(test_data)) / len(test_header)
    }
    if baseline_file is not None:
        metrics = get_baseline_metrics(
            metrics, baseline_file, test_data, samples)
    return metrics


def get_baseline_metrics(metrics, baseline_file, test_data, samples):
    baseline_header, baseline_data = get_medcov_file_data(baseline_file)
    tu.test_iterable_sizes_equal(baseline_header, baseline_data,
                                 name_a="baseline file header", name_b="baseline file data row")
    tu.test_sets_equal(baseline_header, samples, item_str="sample id",
                       name_a="test file header", name_b="sample list")
    n = len(baseline_header)
    error_list = [abs(test_data[i] - baseline_data[i]) for i in range(n)]
    metrics[MEAN_ERROR_KEY] = float(sum(error_list)) / n
    metrics[MAX_ERROR_KEY] = max(error_list)
    return metrics


def get_medcov_file_data(file):
    header = file.readline().strip().split('\t')
    data = [int(x) for x in file.readline().strip().split('\t')]
    return header, data


if __name__ == '__main__':
    main()
