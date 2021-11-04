#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
from collections import deque
import numpy as np
import scipy.stats as ss
import pandas as pd
import pysam


def process_rdtest(rdtest):
    """Standardize rdtest column names"""

    # Drop metadata columns (available from VCF) and rename CNVID
    skip_cols = 'chr Start End Type'.split()
    repl = {'SampleIDs': 'sample', 'CNVID': 'name'}
    rdtest = rdtest.drop(skip_cols, axis=1).rename(columns=repl)

    numeric_cols = 'Median_Power P 2ndMaxP Median_Rank Median_Separation'
    numeric_cols = numeric_cols.split()

    # Replace strings with NA
    for col in numeric_cols:
        repl = ['All_samples_called_CNV_no_analysis',
                'No_samples_for_analysis']
        rdtest[col] = rdtest[col].replace(repl, np.nan).astype(np.float)

    rdtest['log_pval'] = -np.log10(rdtest.P)
    rdtest['log_2ndMaxP'] = -np.log10(rdtest['2ndMaxP'])

    maxp = rdtest.loc[rdtest.log_pval != np.inf, 'log_pval'].max()
    max2p = rdtest.loc[rdtest.log_2ndMaxP != np.inf, 'log_2ndMaxP'].max()

    rdtest.loc[rdtest.log_pval == np.inf, 'log_pval'] = maxp + 5
    rdtest.loc[rdtest.log_2ndMaxP == np.inf, 'log_2ndMaxP'] = max2p + 5

    rdtest.log_pval = rdtest.log_pval.abs()
    rdtest.log_2ndMaxP = rdtest.log_2ndMaxP.abs()

    return rdtest


def process_srtest(srtest):
    # posB-posA dist is in the sum pos column
    dists = srtest.loc[srtest.coord == 'sum', 'name sample pos'.split()].copy()
    dists = dists.rename(columns=dict(pos='dist'))

    metrics = 'log_pval called_median bg_median'.split()

    # remove -0.0 (temporary, should fix in SR-test)
    srtest.log_pval = srtest.log_pval.abs()

    # force one-sided (temporary, should fix in SR-test)
    srtest.loc[srtest.bg_median > srtest.called_median, 'log_pval'] = 0

    srtest = srtest.pivot_table(index='name sample'.split(),
                                values=metrics, columns='coord')
    srtest.columns = ['_'.join(col[::-1]).strip()
                      for col in srtest.columns.values]
    srtest = srtest.reset_index()

    srtest = pd.merge(srtest, dists, on='name sample'.split(), how='left')

    return srtest


def process_petest(petest):
    # remove -0.0 (temporary, should fix in PE-test)
    petest.log_pval = petest.log_pval.abs()

    # force one-sided (temporary, should fix in PE-test)
    petest.loc[petest.bg_median > petest.called_median, 'log_pval'] = 0

    return petest


def preprocess(df, dtype):
    if dtype == 'RD':
        return process_rdtest(df)
    elif dtype == 'SR':
        return process_srtest(df)
    elif dtype == 'PE':
        return process_petest(df)
    else:
        raise Exception('Invalid dtype {0}'.format(dtype))


def add_pesr(evidence):
    evidence['PESR_called_median'] = (evidence['PE_called_median'] +
                                      evidence['SR_sum_called_median'])
    evidence['PESR_bg_median'] = (evidence['PE_bg_median'] +
                                  evidence['SR_sum_bg_median'])

    def calc_p(row):
        pval = ss.poisson.cdf(row.PESR_bg_median, row.PESR_called_median)
        return np.abs(-np.log10(pval))

    evidence['PESR_log_pval'] = evidence.apply(calc_p, axis=1)

    one_sided_mask = (evidence.PESR_bg_median > evidence.PESR_called_median)
    evidence.loc[one_sided_mask, 'PESR_log_pval'] = 0

    return evidence


def make_columns():
    PE_names = ('log_pval called_median bg_median').split()
    PESR_names = ['PESR_' + name for name in PE_names]
    PE_names = ['PE_' + name for name in PE_names]

    SR_names = ('posA_log_pval posB_log_pval sum_log_pval dist '
                'posA_called_median posB_called_median sum_called_median '
                'posA_bg_median posB_bg_median sum_bg_median').split()
    SR_names = ['SR_' + name for name in SR_names]

    RD_names = ('Median_Power P 2ndMaxP Model Median_Rank Median_Separation '
                'log_pval log_2ndMaxP').split()
    RD_names = ['RD_' + name for name in RD_names]

    metadata_names = 'name sample svtype svsize sources'.split()

    return (metadata_names + PE_names + SR_names + PESR_names + RD_names)


def process_metadata(vcf):
    def _get_metadata(record):
        return (record.id,
                record.info['SVTYPE'],
                record.stop - record.pos,
                ','.join(record.info['ALGORITHMS']))

    metadata = [_get_metadata(record) for record in vcf]
    cols = 'name svtype svsize sources'.split()
    metadata = pd.DataFrame(metadata, columns=cols)

    return metadata


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('fout')
    parser.add_argument('-r', '--RDtest')
    parser.add_argument('-s', '--SRtest')
    parser.add_argument('-p', '--PEtest')
    parser.add_argument('-d', '--depth', action='store_true', default=False)
    args = parser.parse_args()

    if args.depth:
        dtypes = 'RD'.split()
    else:
        dtypes = 'PE SR RD'.split()

    evidence = deque()

    for dtype in dtypes:
        dtable = getattr(args, dtype + 'test')
        if dtable is None:
            continue

        df = pd.read_table(dtable)
        df = preprocess(df, dtype)
        df = df.set_index('name sample'.split())
        df = df.rename(columns=lambda c: dtype + '_' + c)
        evidence.append(df)

    evidence = list(evidence)
    evidence = evidence[0].join(evidence[1:], how='left', sort=True)
    evidence = evidence.reset_index()

    # Add SV types
    vcf = pysam.VariantFile(args.vcf)
    metadata = process_metadata(vcf)
    evidence = pd.merge(evidence, metadata, on='name', how='left')

    has_petest = (getattr(args, 'PEtest') is not None)
    has_srtest = (getattr(args, 'SRtest') is not None)
    if not args.depth and has_petest and has_srtest:
        evidence = add_pesr(evidence)

    # Replace infinite log-pvals
    LOG_CEIL = 300
    evidence = evidence.replace(np.inf, LOG_CEIL)

    evidence = evidence.reindex(columns=make_columns())
    evidence.to_csv(args.fout, index=False, sep='\t', na_rep='NA')


if __name__ == '__main__':
    main()
