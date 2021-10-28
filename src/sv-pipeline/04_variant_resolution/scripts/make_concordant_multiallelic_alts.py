#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
When joining VCFs after genotyping, some batches will have different maximum
copy states at multiallelic sites. This script ensures that the maximum
allele count is chosen so all genotype indices are present in the record alts.
"""

import argparse
import pysam


def make_alts(max_alts, is_bca):
    if is_bca:
        alts = tuple(['<CN1>'] +
                     ['<CN{0}>'.format(i) for i in range(2, max_alts + 1)])
    else:
        alts = tuple(['<CN0>'] +
                     ['<CN{0}>'.format(i) for i in range(2, max_alts + 1)])

    return alts


def make_concordant_multiallelic_alts(joined_record, batch_records):
    """
    If any of batch records are multiallelic, update joined record alts and
    genotypes to be consistent with site with highest alternate allele count
    """

    is_multiallelic = False
    is_bca = joined_record.info['SVTYPE'] not in 'DEL DUP'.split()

    max_alts = 0
    non_multiallelic_idxs = []

    for i, record in enumerate(batch_records):
        if record.alts[0].startswith('<CN'):
            is_multiallelic = True
            max_alts = max(max_alts, len(record.alts))
        else:
            non_multiallelic_idxs.append(i)

    # If any records were multiallelic, use maximum alt count
    if is_multiallelic:
        alts = make_alts(max_alts, is_bca)

        stop = joined_record.stop
        joined_record.alts = alts
        joined_record.stop = stop

        # then update genotypes of all non-multiallelic records
        # (only necessary for CNV where first alt is CN0)
        if not is_bca:
            for idx in non_multiallelic_idxs:
                record = batch_records[idx]
                for sample in record.samples:
                    gt = record.samples[sample]['GT']
                    if gt == (0, 0):
                        continue
                    elif gt == (0, 1):
                        joined_record.samples[sample]['GT'] = (0, 2)
                    elif gt == (1, 1):
                        joined_record.samples[sample]['GT'] = (2, 2)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('joined_vcf')
    parser.add_argument('batch_vcfs_list', type=argparse.FileType('r'))
    parser.add_argument('fout')
    args = parser.parse_args()

    joined_vcf = pysam.VariantFile(args.joined_vcf)

    batch_vcfnames = [l.strip() for l in args.batch_vcfs_list.readlines()]
    batch_vcfs = [pysam.VariantFile(fname) for fname in batch_vcfnames]

    fout = pysam.VariantFile(args.fout, 'w', header=joined_vcf.header)

    for joined_record, batch_records in zip(joined_vcf, zip(*batch_vcfs)):
        make_concordant_multiallelic_alts(joined_record, batch_records)
        fout.write(joined_record)


if __name__ == '__main__':
    main()
