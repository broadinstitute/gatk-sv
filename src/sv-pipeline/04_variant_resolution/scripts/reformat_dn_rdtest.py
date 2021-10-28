#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Convert RdTest de novo format to standard format
"""

import argparse
import numpy as np
import pandas as pd


def reformat_dn_rdtest(dn_metrics):
    id_vars = 'chr Start End CNVID Type Family AffectedMember'.split()
    metrics = pd.melt(dn_metrics, id_vars=id_vars)

    metrics['metric'] = metrics.variable.str.split('.').str[1]
    metrics['member'] = metrics.variable.str.split('.').str[0].str.lower()
    metrics.member = metrics.member.replace({'pro': 'p1', 'sib': 's1'})

    # restrict to parents
    #  metrics = metrics.loc[(metrics.member != 'pro') &
    #  (metrics.member != 'sib')].copy()

    # Restrict to parents of candidate de novo variants
    metrics = metrics.loc[(~metrics.AffectedMember.str.contains('fa')) &
                          (~metrics.AffectedMember.str.contains('mo'))].copy()

    metrics['SampleIDs'] = (metrics['Family'].astype(str) +
                            '.' +
                            metrics['member'])

    # Restrict to ID and metric cols
    cols = 'chr Start End CNVID SampleIDs Type metric value'.split()
    metrics = metrics[cols].copy()

    # Pivot out to standard RdTest format
    id_vars = 'chr Start End CNVID SampleIDs Type'.split()
    metrics = metrics.pivot_table(values='value',
                                  index=id_vars,
                                  columns='metric').reset_index()

    # Rename and add standard RdTest columns
    repl = dict(secMaxP='2ndMaxP', rank='Median_Rank', Sep='Median_Separation')
    metrics = metrics.rename(columns=repl)
    metrics['Median_Power'] = np.nan
    metrics['Model'] = 'denovo'

    # Re-order and return
    cols = ('chr Start End CNVID SampleIDs Type Median_Power P 2ndMaxP Model '
            'Median_Rank Median_Separation').split()

    return metrics[cols].copy()


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('dn_metrics', help='De novo RdTest metrics')
    parser.add_argument('fout', help='Reformatted metrics')
    args = parser.parse_args()

    dn_metrics = pd.read_table(args.dn_metrics)
    metrics = reformat_dn_rdtest(dn_metrics)
    metrics.to_csv(args.fout, sep='\t', index=False, na_rep='NA')


if __name__ == '__main__':
    main()
