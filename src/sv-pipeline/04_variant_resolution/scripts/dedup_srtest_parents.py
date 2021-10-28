#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import pandas as pd
import numpy as np
import scipy.stats as ss


def choose_parent_counts(stats):
    # If the children don't agree, pick the strongest
    is_child = (stats['sample'].str.endswith('p1') |
                stats['sample'].str.endswith('s1'))
    children = stats.loc[is_child]
    best_idx = children.groupby('coord')['log_pval'].idxmax()
    best = children.loc[best_idx, 'coord pos'.split()].set_index('coord').pos

    is_parent = (stats['sample'].str.endswith('fa') |
                 stats['sample'].str.endswith('mo'))
    match_bestA = (stats.coord == 'posA') & (stats.pos == best.posA)
    match_bestB = (stats.coord == 'posB') & (stats.pos == best.posB)

    parents = stats.loc[is_parent & (match_bestA | match_bestB)]
    parents = parents.drop_duplicates()

    summed = parents.groupby('sample')['called_median bg_median'.split()].sum()
    pval = ss.poisson.cdf(summed.bg_median, summed.called_median)
    summed['log_pval'] = np.abs(np.log10(pval))

    # Add metadata
    summed['coord'] = 'sum'
    summed['pos'] = 0
    summed['name'] = parents.iloc[0]['name']
    summed = summed.reset_index()

    parents = pd.concat([parents, summed])
    stats = pd.concat([children, parents])

    return stats


def dedup_parents(stats):
    stats['quad'] = stats['sample'].str.split('.').str[0]
    cleaned = stats.groupby('name quad'.split()).apply(choose_parent_counts)
    cleaned = cleaned.reset_index(drop=True)
    cleaned = cleaned[stats.columns].drop('quad', axis=1).head()
    return cleaned


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('statsfile')
    args = parser.parse_args()

    stats = pd.read_table(args.statsfile)
    stats = dedup_parents(stats)
    stats.to_csv(args.fout, index=False, sep='\t', na_rep='NA')


if __name__ == '__main__':
    main()
