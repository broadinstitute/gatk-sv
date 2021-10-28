#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import pandas as pd


def identify_mosaics(metrics):
    metrics = metrics.loc[metrics.Median_Separation.astype(
        str) != 'coverage_failure'].copy()
    metrics.Median_Separation = metrics.Median_Separation.astype(float)

    q1 = metrics.Median_Separation.quantile(0.25)
    q3 = metrics.Median_Separation.quantile(0.75)
    IQR = q3 - q1

    cutoff = q1 - (1.5 * IQR)

    mosaics = metrics.loc[metrics.Median_Separation < cutoff]

    return mosaics['CNVID SampleIDs Median_Separation'.split()]


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('metrics')
    parser.add_argument('fout')
    args = parser.parse_args()

    metrics = pd.read_table(args.metrics)

    mosaics = identify_mosaics(metrics)
    mosaics.to_csv(args.fout, sep='\t', index=False)


if __name__ == '__main__':
    main()
