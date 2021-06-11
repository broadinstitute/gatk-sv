#!/usr/bin/env python

"""
Collect genotyping cutoff metrics. Writes stats to stdout.

Metric format: gt_cutoffs_<COPY_STATE>_<MEAN/SD/CUTOFF>

"""

import argparse
import sys
import pandas as pd

KEY_PREFIX = "gt_cutoffs_"

MEAN_COL = 'mean'
SD_COL = 'sd'
CUTOFFS_COL = 'cutoffs'


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtest gt-cutoffs',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cutoffs', type=str)
    parser.add_argument('type', type=str,
                        help="Genotyping cutoffs type (e.g. depth_depth, depth_pesr, pesr_pesr, pesr_depth")

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    # Get metrics
    df = pd.read_csv(args.cutoffs, sep='\t')
    metrics = get_metrics(df, args.type)

    # Write metrics
    write_metrics(metrics)


def write_metrics(metrics):
    for key in metrics:
        sys.stdout.write("%s\t%s\n" % (key, str(metrics[key])))


def get_metrics(df, type):
    metrics = {}
    for i in range(df.shape[0]):
        row = df.iloc[i]
        name = KEY_PREFIX + type + "_" + str(int(row.copy_state))
        metrics[name + "_mean"] = row[MEAN_COL]
        metrics[name + "_sd"] = row[SD_COL]
        metrics[name + "_cutoff"] = row[CUTOFFS_COL]
    return metrics


if __name__ == '__main__':
    main()
