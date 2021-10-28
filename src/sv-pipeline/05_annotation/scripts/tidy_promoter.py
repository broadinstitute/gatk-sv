#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import pandas as pd


def explode(df, idx_cols, list_col):
    """
    Explode a column of lists into one row per list element
    """
    df = df.set_index(idx_cols)
    df[list_col] = df[list_col].str.split(',')

    df = df[list_col].apply(pd.Series).reset_index()
    df = pd.melt(df, id_vars=idx_cols, value_name=list_col)

    df = df.set_index(idx_cols).drop('variable', axis=1).dropna()
    df = df.reset_index().sort_values(idx_cols)

    return df


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('anno_bed')
    parser.add_argument('fout')
    args = parser.parse_args()

    bed = pd.read_table(args.anno_bed)

    noncoding = bed['name NONCODING'.split()]
    noncoding = explode(noncoding, ['name'], 'NONCODING')

    noncoding['effect'] = 'PROMOTER'
    noncoding = noncoding.rename(columns=dict(NONCODING='gene'))

    cols = 'name effect gene'.split()
    noncoding[cols].to_csv(args.fout, sep='\t', index=False)


if __name__ == '__main__':
    main()
