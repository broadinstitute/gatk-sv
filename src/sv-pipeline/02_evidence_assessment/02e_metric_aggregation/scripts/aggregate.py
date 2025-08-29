#!/usr/bin/env python
# -*- coding: utf-8 -*-
# freq

"""
"""

import argparse
from collections import deque, defaultdict
import numpy as np
import scipy.stats as ss
import pandas as pd
import pysam
import pybedtools as pbt
import svtk.utils as svu


def process_rdtest(rdtest):
    """Standardize rdtest column names"""

    # Drop metadata columns (available from VCF) and rename CNVID
    skip_cols = 'chr Start End SampleIDs Type'.split()
    rdtest = rdtest.drop(skip_cols, axis=1).rename(columns={'CNVID': 'name'})

    numeric_cols = 'Median_Power P 2ndMaxP Median_Rank Median_Separation'
    numeric_cols = numeric_cols.split()

    # Replace strings with NA
    for col in numeric_cols:
        repl = ['All_samples_called_CNV_no_analysis',
                'No_samples_for_analysis',
                'coverage_failure']
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
    metrics = 'log_pval called_median bg_median bg_frac pos'.split()

    # remove -0.0 (temporary, should fix in SR-test)
    srtest.log_pval = srtest.log_pval.abs()

    srtest.pos = srtest.pos.astype(int)

    # force one-sided (temporary, should fix in SR-test)
    srtest.loc[srtest.bg_median > srtest.called_median, 'log_pval'] = 0

    srtest = srtest.pivot_table(index='name', values=metrics, columns='coord')
    srtest.columns = ['_'.join(col[::-1]).strip()
                      for col in srtest.columns.values]
    srtest = srtest.reset_index()

    return srtest


def process_petest(petest):
    # remove -0.0 (temporary, should fix in PE-test)
    petest.log_pval = petest.log_pval.abs()

    # force one-sided (temporary, should fix in PE-test)
    petest.loc[petest.bg_median > petest.called_median, 'log_pval'] = 0

    return petest


def process_baftest(baftest):
    skip_cols = 'chrom start end samples svtype'.split()
    baftest = baftest.drop(skip_cols, axis=1)

    baftest['KS_log_pval'] = (- np.log10(baftest.KS_pval)).abs()
    baftest['del_loglik'] = -baftest.del_loglik

    repl = 'Potential ROHregion or reference error'
    baftest.delstat = baftest.delstat.replace(repl, 'Ref_error')

    return baftest


def preprocess(df, dtype):
    if dtype == 'RD':
        return process_rdtest(df)
    elif dtype == 'SR':
        return process_srtest(df)
    elif dtype == 'PE':
        return process_petest(df)
    elif dtype == 'BAF':
        return process_baftest(df)
    else:
        return df


def _is_parent(s):
    return s.endswith('fa') or s.endswith('mo')


def _is_child(s):
    return s.endswith('p1') or s.endswith('s1') or s.endswith('pb')


def fam_info_readin(fam_file):
    fin = open(fam_file)
    # samp_pedi_hash = {}
    [fam, samp, fa, mo] = [[], [], [], []]
    for line in fin:
        pin = line.strip().split()
        fam.append(pin[0])
        samp.append(pin[1])
        fa.append(pin[2])
        mo.append(pin[3])
    fin.close()
    return [fam, samp, fa, mo]


def process_metadata(variants, bed=False, batch_list=None, outlier_sample_ids=None):
    if bed:
        samples = [s.strip() for s in batch_list.readlines()]
    else:
        samples = list(variants.header.samples)

    outlier_set = set()
    if outlier_sample_ids:
        with open(outlier_sample_ids, 'r') as f:
            outlier_set = set(line.strip() for line in f)

    # parents = [s for s in samples if _is_parent(s)]
    # children = [s for s in samples if _is_child(s)]
    # n_parents = len(parents)
    # n_children = len(children)

    called_counts = dict()
    called_samples = dict()
    for svtype in 'DEL DUP INV BND INS'.split():
        # Counts of variants per sample
        called_counts[svtype] = defaultdict(int)

        # List of variants specific to each sample
        called_samples[svtype] = defaultdict(list)

    metadata = deque()
    for variant in variants:
        # bed record
        if bed:
            if variant.startswith('#'):
                continue
            data = variant.strip().split()
            chrom = data[0]
            start = int(data[1])
            end = int(data[2])
            called = data[4].split(',')
            name = data[3]
            svtype = data[5]
            svlen = int(data[2]) - int(data[1])
        # VCF record
        else:
            chrom = variant.chrom
            start = variant.pos
            end = variant.info['END2'] if variant.info['SVTYPE'] == 'BND' else variant.stop
            called = svu.get_called_samples(variant)
            name = variant.id
            svtype = variant.info['SVTYPE']
            svlen = variant.info['SVLEN']

        # Only use start/end for seg dup coverage. if it's a tloc,
        # we don't care so we can just set its "END" to pos + 1
        if end <= start:
            end = start + 1

        # Calculate VF
        vf = len(called) / len(samples)

        # Increment counts of variants per sample
        for s in called:
            called_counts[svtype][s] += 1

        # Track called samples for outlier filtering
        called_samples[svtype][name] = set(called)

        dat = [chrom, start, end, name, svtype, svlen, vf]
        metadata.append(dat)

    metadata = np.array(metadata)
    cols = 'chrom start end name svtype svsize vf'.split()
    metadata = pd.DataFrame(metadata, columns=cols)

    # Flag variants specific to outlier samples
    metadata['is_outlier_specific'] = False
    if len(outlier_set) > 0:
        for variants in called_samples.values():
            for name, called in variants.items():
                if called and called.issubset(outlier_set):
                    metadata.loc[metadata.name == name, 'is_outlier_specific'] = True

    for col in 'start end svsize'.split():
        metadata[col] = metadata[col].astype(int)

    return metadata


def add_pesr(evidence):
    evidence['PESR_called_median'] = (evidence['PE_called_median'] +
                                      evidence['SR_sum_called_median'])
    evidence['PESR_bg_median'] = (evidence['PE_bg_median'] +
                                  evidence['SR_sum_bg_median'])
    evidence['PESR_bg_frac'] = (evidence['PESR_bg_median'] /
                                (evidence['PESR_bg_median'] + evidence['PESR_called_median']))

    def calc_p(row):
        pval = ss.poisson.cdf(row.PESR_bg_median, row.PESR_called_median)
        return np.abs(-np.log10(pval))

    evidence['PESR_log_pval'] = evidence.apply(calc_p, axis=1)

    one_sided_mask = (evidence.PESR_bg_median > evidence.PESR_called_median)
    evidence.loc[one_sided_mask, 'PESR_log_pval'] = 0

    return evidence


def make_columns():
    PE_names = ('log_pval called_median bg_median bg_frac').split()
    PESR_names = ['PESR_' + name for name in PE_names]
    PE_names = ['PE_' + name for name in PE_names]

    SR_names = ('posA_log_pval posB_log_pval sum_log_pval '
                'posA_called_median posB_called_median sum_called_median '
                'posA_bg_median posB_bg_median sum_bg_median '
                'posA_bg_frac posB_bg_frac sum_bg_frac '
                'posA_pos posB_pos').split()
    SR_names = ['SR_' + name for name in SR_names]

    BAF_names = ('delstat snp_ratio del_loglik dupstat KS_stat KS_log_pval '
                 'total_case_snps total_snps n_nonROH_cases n_samples '
                 'mean_control_snps n_nonROH_controls n_controls').split()
    BAF_names = ['BAF_' + name for name in BAF_names]

    RD_names = ('Median_Power P 2ndMaxP Model Median_Rank Median_Separation '
                'log_pval log_2ndMaxP').split()
    RD_names = ['RD_' + name for name in RD_names]

    #  metadata_names = 'name svtype svsize parental_vf child_vf inh_rate'.split()
    metadata_names = 'name chrom svtype svsize vf poor_region_cov rmsk is_outlier_specific'.split()

    return (metadata_names + PE_names + SR_names + PESR_names + RD_names +
            BAF_names)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--variants', required=True, help='Default VCF')
    parser.add_argument('-r', '--RDtest')
    parser.add_argument('-b', '--BAFtest')
    parser.add_argument('-s', '--SRtest')
    parser.add_argument('-p', '--PEtest')
    parser.add_argument('-o', '--outlier-sample-ids')
    parser.add_argument('--batch-list', type=argparse.FileType('r'))
    parser.add_argument('--segdups', required=True)
    parser.add_argument('--rmsk', required=True)
    parser.add_argument('--fam')
    parser.add_argument('-d', '--bed', action='store_true', default=False)
    parser.add_argument('fout')
    args = parser.parse_args()

    if args.bed:
        if not hasattr(args, 'batch_list'):
            raise Exception('batch list must be specified when passing a bed')
        variants = open(args.variants)
        dtypes = 'RD BAF'.split()
    else:
        variants = pysam.VariantFile(args.variants)
        dtypes = 'PE SR RD BAF'.split()

    outlier_sample_ids = None
    if args.outlier_sample_ids:
        outlier_sample_ids = args.outlier_sample_ids

    metadata = process_metadata(variants, args.bed, args.batch_list, outlier_sample_ids)

    # Calculate segdup coverage
    bt = pbt.BedTool.from_dataframe(metadata['chrom start end'.split()])
    segdups = pbt.BedTool(args.segdups)
    cov = bt.coverage(segdups).to_dataframe()
    metadata['poor_region_cov'] = cov.thickStart

    # Check if endpoints are in repeat-masked sequence
    starts = metadata['chrom start end name'.split()].copy()
    starts['end'] = starts['start'] + 1
    ends = metadata['chrom start end name'.split()].copy()
    ends['start'] = ends['end'] - 1
    endpoints = pd.concat([starts, ends])
    bt = pbt.BedTool.from_dataframe(endpoints)
    rmsk = pbt.BedTool(args.rmsk)
    sect = bt.intersect(rmsk, u=True)
    rmasked_names = [i.fields[3] for i in sect.intervals]
    metadata['rmsk'] = metadata.name.isin(rmasked_names)

    metadata = metadata.set_index('name')

    evidence = deque()

    for dtype in dtypes:
        dtable = getattr(args, dtype + 'test')
        if dtable is None:
            continue

        df = pd.read_table(dtable)

        df = preprocess(df, dtype)
        df = df.rename(columns=lambda c: dtype + '_' + c if c != 'name' else c)
        df = df.set_index('name')
        evidence.append(df)

    evidence = list(evidence)
    evidence = metadata.join(evidence, how='outer', sort=True)
    evidence = evidence.reset_index().rename(columns={'index': 'name'})

    has_petest = (getattr(args, 'PEtest') is not None)
    has_srtest = (getattr(args, 'SRtest') is not None)
    if not args.bed and has_petest and has_srtest:
        evidence = add_pesr(evidence)

    # Replace infinite log-pvals
    LOG_CEIL = 300
    evidence = evidence.replace(np.inf, LOG_CEIL)

    evidence = evidence.reindex(columns=make_columns())
    evidence.to_csv(args.fout, index=False, sep='\t', na_rep='NA')


if __name__ == '__main__':
    main()
