#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
All the hard-coded, SFARI-specific filtering
"""

import argparse
import sys
import pysam
import svtk.utils as svu


# Samples with allosomal aneuploidies
ANEUPLOIDIES = set([
    '13706.p1',
    '14100.p1',
    '11554.s1',
    '12412.p1',
    '12860.p1'
])

# Dosage outliers
OUTLIERS = set([
    '11572.s1',
    '12175.s1',
    '12568.s1',
    '12680.fa',
    '13023.mo'
])

# Samples with Robertsonian translocations
ROBERTSONIANS = set([
    '11219.s1',
    '13424.fa',
    '14005.s1',
])

DEPTH_EXCLUDED = OUTLIERS.union(ROBERTSONIANS)


def sfari_filters(vcf):
    for record in vcf:
        called = set(svu.get_called_samples(record))

        # Remove allosomal calls in aneuploidy samples
        if record.chrom in 'X Y'.split():
            if called.issubset(ANEUPLOIDIES):
                continue
            else:
                for sample in called.intersection(ANEUPLOIDIES):
                    svu.set_null(record, sample)

        # Remove depth-only events in Robertsonian translocation cases
        if record.info['SOURCES'] == ('depth',):
            if called.issubset(DEPTH_EXCLUDED):
                continue

        # Check variant wasn't only in aneuploidies and Robertsonians
        #  called = set(svu.get_called_samples(record))
        #  if len(called) == 0:
            #  continue

        # Remove variants specific to dosage outliers
        #  if called.issubset(OUTLIERS):
            #  continue

        yield record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('fout')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    for record in sfari_filters(vcf):
        fout.write(record)


if __name__ == '__main__':
    main()
