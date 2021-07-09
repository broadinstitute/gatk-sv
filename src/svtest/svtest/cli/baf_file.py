#!/usr/bin/env python

"""
Collect BAF file metrics. Writes stats to stdout.

Metrics:
  baf_qN_<sample>    : Nth percentile BAF
  baf_count_<sample> : total BAF record count

If the provided sample list has more than one id, <sample> will be "_merged".

"""

import gzip
import argparse
import sys
import numpy as np
import svtest.utils.TestUtils as tu
import svtest.utils.IOUtils as iou

Q25_KEY = "baf_q25"
Q50_KEY = "baf_q50"
Q75_KEY = "baf_q75"
COUNT_KEY = "baf_count"

EXPECTED_COLUMNS = 4


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtest raw-baf',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('baf_file', type=str)
    parser.add_argument('sample_list', type=str)

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    # Read file
    with gzip.open(args.baf_file, mode='rb') as fbaf:
        metrics = get_metrics(fbaf, args.sample_list)

    # Write metrics
    write_metrics(metrics)


def write_metrics(metrics):
    for key in metrics:
        sys.stdout.write("%s\t%s\n" % (key, str(metrics[key])))


def get_metrics(baf_file, sample_list):
    samples = iou.read_samples_list(sample_list)
    samples_set = set(samples)

    data = []
    for line in baf_file:
        tokens = line.decode().strip().split('\t')
        test_record(tokens, samples_set)
        baf = float(tokens[2])
        data.append(baf)
    arr = np.asarray(data)
    quantiles = np.quantile(arr, [0.25, 0.50, 0.75])

    if len(samples) == 1:
        metric_suffix = "_" + samples[0]
    else:
        metric_suffix = "_merged"

    return {
        Q25_KEY + metric_suffix: quantiles[0],
        Q50_KEY + metric_suffix: quantiles[1],
        Q75_KEY + metric_suffix: quantiles[2],
        COUNT_KEY + metric_suffix: len(arr)
    }


def test_record(columns, samples):
    tu.test_iterable_size(columns, EXPECTED_COLUMNS)
    tu.test_is_float(columns, 2)
    tu.test_column_in_iterable(columns, 3, samples)


if __name__ == '__main__':
    main()
