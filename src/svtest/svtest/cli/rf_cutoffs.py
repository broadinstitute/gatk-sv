#!/usr/bin/env python

"""
Collect random forest cutoff metrics. Writes stats to stdout.

Metric format: rf_cutoff_<TEST>_<SVTYPE>_<ALGTYPE>_min<MIN>_max<MAX>

"""

import argparse
import sys
import pandas as pd

KEY_PREFIX = "rf_cutoff_"


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtest rf-cutoffs',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cutoffs', type=str)

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    # Get metrics
    df = pd.read_csv(args.cutoffs, sep='\t')
    metrics = get_metrics(df)

    # Write metrics
    write_metrics(metrics)


def write_metrics(metrics):
    for key in metrics:
        sys.stdout.write("%s\t%s\n" % (key, str(metrics[key])))


def get_metrics(df):
    metrics = {}
    for i in range(df.shape[0]):
        row = df.iloc[i]
        name = KEY_PREFIX + \
            "_".join([row.test, row.svtype, row.metric, row.algtype])
        if not pd.isna(row.min_svsize):
            name += "_min" + str(int(row.min_svsize))
        if not pd.isna(row.max_svsize):
            name += "_max" + str(int(row.max_svsize))
        metrics[name] = row.cutoff
    return metrics


if __name__ == '__main__':
    main()
