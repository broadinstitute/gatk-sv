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


def process_metadata(vcf, outlier_sample_ids):

    n_samples = len(vcf.header.samples)
    called_counts = dict()
    called_samples = dict()

    outlier_set = set()
    if outlier_sample_ids:
        with open(outlier_sample_ids, 'r') as f:
            outlier_set = set(line.strip() for line in f)

    for svtype in 'DEL DUP INV BND INS'.split():
        # Counts of variants per sample
        called_counts[svtype] = defaultdict(int)

        # List of variants specific to each sample
        called_samples[svtype] = defaultdict(list)

    stats_int = ['SR1POS', 'SR2POS']
    stats_float = ['BAF_KS_Q', 'SR1Q', 'SR1CS', 'SR2Q', 'SR2CS', 'SRQ', 'SRCS', 'BAF_HET_RATIO', 'BAF_KS_STAT',
                   'PEQ', 'PECS', 'PESRQ', 'PESRCS', 'RDQ', 'RD_P2', 'RD_MEDIAN_SEPARATION']
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
        rmsk = variant.info['OVERLAP_FRAC_RMSK'] > 0 if 'OVERLAP_FRAC_RMSK' in variant.info.keys() else None
        segdup = variant.info['OVERLAP_FRAC_SEGDUP'] if 'OVERLAP_FRAC_SEGDUP' in variant.info.keys() else None

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
    outlier_specific_vids = set()
    if len(outlier_set) > 0:
        for variants in called_samples.values():
            for name, called in variants.items():
                if called and called.issubset(outlier_set):
                    outlier_specific_vids.add(name)
    metadata['is_outlier_specific'] = metadata.name.isin(outlier_specific_vids)

    for col in 'start end svsize'.split():
        metadata[col] = metadata[col].astype(int)

    metadata.set_index(keys='name', inplace=True)
    return metadata


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--variants', required=True, help='Input VCF')
    parser.add_argument('-o', '--outlier-sample-ids', help='Path to file containing outlier sample IDs')
    parser.add_argument('fout')
    args = parser.parse_args()

    outlier_sample_ids = None
    if args.outlier_sample_ids:
        outlier_sample_ids = args.outlier_sample_ids

    vcf = pysam.VariantFile(args.variants)
    evidence = process_metadata(vcf, outlier_sample_ids=outlier_sample_ids)
    evidence.to_csv(args.fout, index=True, sep='\t', na_rep='NA')


if __name__ == '__main__':
    main()
