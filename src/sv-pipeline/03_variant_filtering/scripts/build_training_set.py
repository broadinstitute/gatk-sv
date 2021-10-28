#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import os
import pandas as pd


def build_baf_training_set(metrics):
    MIN_SIZE = 5000
    MAX_BLACKLIST_COV = 0.3
    MAX_FAIL_SEP = 0.15
    MIN_PASS_SEP = 0.40

    # Remove BAF failures (variants where BAF can not be assessed due to
    # insufficient informative SNPs or a ROH)
    del_cols = 'BAF_snp_ratio BAF_del_loglik'.split()
    del_mask = ((metrics.svtype == 'DEL') &
                (~metrics[del_cols].isnull().any(axis=1)))

    dup_cols = 'BAF_KS_stat BAF_KS_log_pval'.split()
    dup_mask = ((metrics.svtype == 'DUP') &
                (~metrics[dup_cols].isnull().any(axis=1)))

    # Restrict BAF testing to variants >5 kb
    size_mask = (metrics['size'] >= MIN_SIZE)

    # Filter to all variants eligible for classification
    all_metrics = metrics[(del_mask | dup_mask) & size_mask].copy()

    # Training variants must have <30% blacklist coverage and
    # be below 15% or above 40% separation
    cov_mask = (metrics['poor_region_cov'] < MAX_BLACKLIST_COV)
    sep_mask = ((metrics['RD_Median_Separation'] < MAX_FAIL_SEP) |
                (metrics['RD_Median_Separation'] >= MIN_PASS_SEP))

    # Filter to training set
    train = metrics[(del_mask | dup_mask) & size_mask &
                    sep_mask & cov_mask].copy()

    # Label training classes based on separation
    train.loc[train.RD_Median_Separation >= MIN_PASS_SEP, 'Status'] = 'Pass'
    train.loc[train.RD_Median_Separation < MAX_FAIL_SEP, 'Status'] = 'Fail'

    # TEMPORARY: write to file directly
    os.makedirs('baf_rf', exist_ok=True)
    id_cols = 'name Status'.split()

    del_train = train.loc[train.svtype == 'DEL', id_cols + del_cols]
    del_train.to_csv('baf_metrics/Phase1.del.train.metrics',
                     sep='\t', index=False, na_rep='NA')
    del_all = all_metrics.loc[all_metrics.svtype == 'DEL']
    del_all.to_csv('baf_metrics/Phase1.del.all.metrics',
                   sep='\t', index=False, na_rep='NA')

    dup_train = train.loc[train.svtype == 'DUP', id_cols + dup_cols]
    dup_train.to_csv('baf_metrics/Phase1.dup.train.metrics',
                     sep='\t', index=False, na_rep='NA')
    dup_all = all_metrics.loc[all_metrics.svtype == 'DUP']
    dup_all.to_csv('baf_metrics/Phase1.dup.all.metrics',
                   sep='\t', index=False, na_rep='NA')


def build_training_set(metrics):
    build_baf_training_set(metrics)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('metrics')
    args = parser.parse_args()

    metrics = pd.read_table(args.metrics)

    build_training_set(metrics)


if __name__ == '__main__':
    main()
