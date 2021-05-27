#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2018 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import numpy as np
import pandas as pd


def get_column_count(filename):
    with open(filename, "r") as f:
        header = f.readline()
        return len(header.split("\t"))


def get_rows_to_skip(filename):
    rows = []
    with open(filename, "r") as f:
        c = 0
        for line in f:
            if "." in line:
                rows.append(c)
            c += 1
    return rows


def get_table(filename, keep_rows_with_missing_variants):
    column_count = get_column_count(filename)
    dtypes = {0: np.str, 1: np.int, 2: np.int, 3: np.str}
    dtypes.update(dict((i + 4, np.uint8) for i in range(column_count)))

    skip_rows = []
    if not keep_rows_with_missing_variants:
        skip_rows = get_rows_to_skip(filename)

    return pd.read_table(filename, dtype=dtypes, skiprows=skip_rows)


def get_stacked(table, value_col_label="genotype"):
    table = table.set_index(["chr", "start", "end", "cnvID"])
    table.columns = table.columns.astype('category')
    table = (table.stack()
             .rename_axis(["chr", "start", "end", "cnvID", "sample"])
             .rename(value_col_label)
             .reset_index())
    return table


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('genotypes')
    parser.add_argument('GQ')
    parser.add_argument('fout')
    parser.add_argument('--keep_rows_with_missing_variants', action='store_true')
    args = parser.parse_args()

    # Load and melt genotypes
    gt = get_table(args.genotypes, args.keep_rows_with_missing_variants).drop_duplicates()
    gt = get_stacked(gt)

    # Round genotype copy states
    #  gt['genotype'] = gt['genotype'].round().astype(int)

    # Pivot back out so input genotype matrix is all integers
    #  pivot = gt.pivot_table(index='chr start end cnvID'.split(),
                           #  columns='sample', values='genotype').reset_index()
    #  samples = [c for c in pivot.columns if c not in 'chr start end cnvID'.split()]
    #  pivot[samples] = pivot[samples].astype(int)
    #  pivot.to_csv(args.genotypes, index=False, sep='\t')

    # Load and melt GQ
    gq = get_table(args.GQ, args.keep_rows_with_missing_variants)
    gq = get_stacked(gq, "GQ")

    # Merge genotypes with GQ
    gt = pd.merge(gt, gq, on='chr start end cnvID sample'.split(), how='left')
    gt.to_csv(args.fout, header=False, index=False, sep='\t')


if __name__ == '__main__':
    main()
