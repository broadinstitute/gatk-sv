#!/usr/bin/env python
# -*- coding: utf-8 -*-
# freq

"""
"""

import argparse
from collections import deque, defaultdict
import numpy as np
import pandas as pd
import pysam
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
        rdtest[col] = rdtest[col].replace(repl, np.nan).astype(float)

    rdtest['log_pval'] = -np.log10(rdtest.P)
    rdtest['log_2ndMaxP'] = -np.log10(rdtest['2ndMaxP'])

    maxp = rdtest.loc[rdtest.log_pval != np.inf, 'log_pval'].max()
    max2p = rdtest.loc[rdtest.log_2ndMaxP != np.inf, 'log_2ndMaxP'].max()

    rdtest.loc[rdtest.log_pval == np.inf, 'log_pval'] = maxp + 5
    rdtest.loc[rdtest.log_2ndMaxP == np.inf, 'log_2ndMaxP'] = max2p + 5

    rdtest.log_pval = rdtest.log_pval.abs()
    rdtest.log_2ndMaxP = rdtest.log_2ndMaxP.abs()

    rdtest = rdtest.rename(columns=lambda c: 'RD_' + c if c != 'name' else c)
    rdtest.set_index('name', inplace=True)

    return rdtest


def process_metadata(vcf):

    n_samples = len(vcf.header.samples)
    called_counts = dict()
    called_samples = dict()

    for svtype in 'DEL DUP INV BND INS'.split():
        # Counts of variants per sample
        called_counts[svtype] = defaultdict(int)

        # List of variants specific to each sample
        called_samples[svtype] = defaultdict(list)

    stats_int = ['BAF_KS_Q', 'SR1Q', 'SR1CS', 'SR2Q', 'SR2CS', 'SRQ', 'SRCS', 'SR1POS', 'SR2POS', 'PEQ', 'PECS',
                 'PESRQ', 'PESRCS']
    stats_float = ['BAF_HET_RATIO', 'BAF_KS_STAT']
    metadata = deque()
    for variant in vcf:
        chrom = variant.chrom
        start = variant.pos
        end = variant.stop
        name = variant.id
        svtype = variant.info['SVTYPE']
        called = svu.get_called_samples(variant)
        if svtype == 'BND':
            svlen = -1
        elif svtype == 'INS':
            svlen = variant.info.get('SVLEN', -1)
        else:
            svlen = end - start

        # Only use start/end for seg dup coverage. if it's a tloc,
        # we don't care so we can just set its "END" to pos + 1
        if end <= start:
            end = start + 1

        # Calculate VF
        vf = len(called) / n_samples

        # Increment counts of variants per sample
        for s in called:
            called_counts[svtype][s] += 1

        # Track called samples for outlier filtering
        called_samples[svtype][name] = set(called)

        # Repeatmasker / segdup track overlap
        rmsk = variant.info['OVERLAP_FRAC_RMSK'] > 0
        segdup = variant.info['OVERLAP_FRAC_SEGDUP']

        dat = [chrom, start, end, name, svtype, svlen, vf, rmsk, segdup]
        dat.extend([int(variant.info[stat]) if stat in variant.info.keys() and variant.info[stat] is not None else None for stat in stats_int])
        dat.extend([float(variant.info[stat]) if stat in variant.info.keys() and variant.info[stat] is not None else None for stat in stats_float])
        metadata.append(dat)

    metadata = np.array(metadata)
    cols = 'chrom start end name svtype svsize vf rmsk poor_region_cov'.split()
    cols.extend(stats_int)
    cols.extend(stats_float)
    metadata = pd.DataFrame(metadata, columns=cols)

    # Flag variants specific to outlier samples
    metadata['is_outlier_specific'] = False
    for svtype in 'DEL DUP INV BND INS'.split():
        counts = pd.DataFrame.from_dict(called_counts[svtype], orient='index')\
                             .reset_index()\
                             .rename(columns={'index': 'sample', 0: 'var_count'})
        if counts.shape[0] == 0:
            continue

        q1 = counts.var_count.quantile(0.25)
        q3 = counts.var_count.quantile(0.75)
        thresh = q3 + 1.5 * (q3 - q1)
        outliers = counts.loc[counts.var_count >= thresh, 'sample'].values

        flagged = []
        for var_name, samples in called_samples[svtype].items():
            if samples.issubset(outliers):
                flagged.append(var_name)
        metadata.loc[metadata.name.isin(flagged), 'is_outlier_specific'] = True

    for col in 'start end svsize'.split():
        metadata[col] = metadata[col].astype(int)

    metadata.set_index(keys='name', inplace=True)
    return metadata


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--variants', required=True, help='Input VCF')
    parser.add_argument('-r', '--rdtest', help='RDtest table')
    parser.add_argument('fout')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.variants)
    evidence = process_metadata(vcf)

    # Parse RDTest table
    rd_path = getattr(args, 'rdtest')
    if rd_path is not None:
        rd_df = pd.read_table(rd_path)
        rd_df = process_rdtest(rd_df)
        evidence = evidence.join(rd_df, how='outer', sort=True)

    # Replace infinite log-pvals
    LOG_CEIL = 300
    evidence = evidence.replace(np.inf, LOG_CEIL)

    evidence.to_csv(args.fout, index=True, sep='\t', na_rep='NA')


if __name__ == '__main__':
    main()
