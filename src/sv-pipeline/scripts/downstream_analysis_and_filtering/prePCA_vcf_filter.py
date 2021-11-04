#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Filter vcf to clean, autosomal, biallelic sites prior to cohort-wide PCA
"""


import argparse
import sys
import pysam


NULL_GTs = [(0, 0), (None, None), (0, ), (None, ), (None, 2)]


def get_call_rate(record):
    n_samples = len(record.samples)
    n_non_null = len(
        [s for s in record.samples if record.samples[s]['GT'] not in NULL_GTs])

    callrate = n_non_null / n_samples

    return callrate


def filter_vcf(vcf, fout_common, minAF=0.01, fout_all=None, minCallRate=0.99):
    for record in vcf:
        # #Do not include UNRESOLVED variants
        # if 'UNRESOLVED' in record.info.keys() \
        # or 'UNRESOLVED_TYPE' in record.info.keys() \
        # or 'UNRESOLVED' in record.filter:
        #     continue

        # #Do not include multiallelic variants
        # if 'MULTIALLELIC' in record.info.keys() \
        # or 'MULTIALLELIC' in record.filter \
        # or len(record.alts) > 1:
        #     continue

        # Do not include variants on X or Y
        allosomes = ['X', 'Y', 'chrX', 'chrY']
        if record.chrom in allosomes:
            continue

        # Only include PASS variants
        if 'PASS' not in record.filter:
            continue

        # Only include variants with ≥minCallRate
        if get_call_rate(record) < minCallRate:
            continue

        # Only keep common variants
        if 'AF' in record.info.keys():
            if record.info['AF'][0] >= minAF:
                fout_common.write(record)

        # Write AF-unfiltered variants, if optioned
        if fout_all is not None:
            fout_all.write(record)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('fout')
    parser.add_argument('--minAF', type=float, default=0.01,
                        help='Minimum allele frequency. [0.01]')
    parser.add_argument('--minCallRate', type=float, default=0.99,
                        help='Minimum call rate. [0.99]')
    parser.add_argument('--noAFoutput', default=None,
                        help='Output file for all variants unfiltered on AF.')

    args = parser.parse_args()

    # Open input VCF
    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    header = vcf.header

    # Open outut VCFs
    fout_common = pysam.VariantFile(args.fout, 'w', header=header)
    if args.noAFoutput is not None:
        fout_all = pysam.VariantFile(args.noAFoutput, 'w', header=header)
    else:
        fout_all = None

    # Filter VCF
    filter_vcf(vcf, fout_common, args.minAF, fout_all, args.minCallRate)


if __name__ == '__main__':
    main()
