#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import pandas as pd
import pysam
import svtk.utils as svu
from collections import namedtuple


def select_mosaic_candidates(vcf, metrics, mode, cutoffs):
    candidates = metrics.loc[metrics.svtype.isin('DEL DUP'.split())]
    candidates = candidates.loc[candidates.svsize >= 5000]
    candidates = candidates.loc[candidates.poor_region_cov < 0.3]

    if mode == 'pesr':
        candidates = candidates.loc[candidates.RD_log_pval >= cutoffs.min_pval]
        candidates = candidates.loc[candidates.RD_log_2ndMaxP >=
                                    cutoffs.min_secondp]
    else:
        pval_filter = (((candidates.svtype == 'DEL') & (candidates.RD_log_pval >= cutoffs.del_min_pval)) |
                       ((candidates.svtype == 'DUP') & (candidates.RD_log_pval >= cutoffs.dup_min_pval)))

        candidates = candidates.loc[pval_filter]

    candidate_IDs = candidates.name.values

    for record in vcf:
        called = svu.get_called_samples(record)
        if len(called) == 1 and record.id in candidate_IDs and record.info['SVTYPE'] in 'DEL DUP'.split():
            yield record


def choose_cutoffs(cutoffs, mode):
    cutoffs = cutoffs.loc[cutoffs['test'] == 'RD']

    if mode == 'depth':
        cutoffs = cutoffs.loc[cutoffs.algtype == 'Depth']
        cutoffs = cutoffs.loc[cutoffs.metric == 'RD_log_pval']
        min_del = cutoffs.loc[cutoffs.svtype == 'DEL', 'cutoff'].iloc[0]
        min_dup = cutoffs.loc[cutoffs.svtype == 'DUP', 'cutoff'].iloc[0]
        Cutoffs = namedtuple('Cutoffs', 'del_min_pval dup_min_pval'.split())
        cutoffs = Cutoffs(min_del, min_dup)
    else:
        cutoffs = cutoffs.loc[cutoffs.algtype == 'PESR']
        cutoffs = cutoffs.loc[cutoffs.min_svsize == 1000]
        min_pval = cutoffs.loc[cutoffs.metric ==
                               'RD_log_pval', 'cutoff'].iloc[0]
        min_secondp = cutoffs.loc[cutoffs.metric ==
                                  'RD_log_2ndMaxP', 'cutoff'].iloc[0]
        Cutoffs = namedtuple('Cutoffs', 'min_pval min_secondp'.split())
        cutoffs = Cutoffs(min_pval, min_secondp)

    return cutoffs


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('metrics')
    parser.add_argument('cutoffs')
    parser.add_argument('fout')
    parser.add_argument('--mode', choices='depth pesr'.split())
    args = parser.parse_args()

    metrics = pd.read_table(args.metrics)
    vcf = pysam.VariantFile(args.vcf)
    fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    cutoffs = pd.read_table(args.cutoffs)
    cutoffs = choose_cutoffs(cutoffs, args.mode)

    for record in select_mosaic_candidates(vcf, metrics, args.mode, cutoffs):
        fout.write(record)


if __name__ == '__main__':
    main()
