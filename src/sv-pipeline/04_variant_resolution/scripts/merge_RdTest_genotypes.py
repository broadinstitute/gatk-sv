#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2018 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('genotypes')
    parser.add_argument('GQ')
    parser.add_argument('fout')
    args = parser.parse_args()

    # Load and melt genotypes
    gt = pd.read_table(args.genotypes).drop_duplicates()
    gt = pd.melt(gt, id_vars='chr start end cnvID'.split(),
                 var_name='sample', value_name='genotype')

    # Round genotype copy states
    #  gt['genotype'] = gt['genotype'].round().astype(int)

    # Pivot back out so input genotype matrix is all integers
    #  pivot = gt.pivot_table(index='chr start end cnvID'.split(),
                           #  columns='sample', values='genotype').reset_index()
    #  samples = [c for c in pivot.columns if c not in 'chr start end cnvID'.split()]
    #  pivot[samples] = pivot[samples].astype(int)
    #  pivot.to_csv(args.genotypes, index=False, sep='\t')

    # Load and melt GQ
    gq = pd.read_table(args.GQ)
    gq = pd.melt(gq, id_vars='chr start end cnvID'.split(),
                 var_name='sample', value_name='GQ')

    # Merge genotypes with GQ
    gt = pd.merge(gt, gq, on='chr start end cnvID sample'.split(), how='left')
    gt.to_csv(args.fout, header=False, index=False, sep='\t')


if __name__ == '__main__':
    main()
