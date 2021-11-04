#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Clean up FILTER fields & remove MEMBERS tag at the end of cleanVCF
"""

import argparse
import sys
import pysam


def clean_filt(vcf, fout):

    for record in vcf:

        filters = record.filter

        if 'MULTIALLELIC' in filters \
                or record.info['SVTYPE'] == 'MCNV':
            record.filter.add('MULTIALLELIC')

        if 'UNRESOLVED' in filters \
                or 'UNRESOLVED' in record.info.keys() \
                or 'UNRESOLVED_TYPE' in record.info.keys():
            record.filter.add('UNRESOLVED')
            if 'UNRESOLVED' in record.info.keys():
                record.info.pop('UNRESOLVED')

        if 'HIGH_SR_BACKGROUND' in filters \
                or 'HIGH_SR_BACKGROUND' in record.info.keys():
            record.filter.add('HIGH_SR_BACKGROUND')
            if 'HIGH_SR_BACKGROUND' in record.info.keys():
                record.info.pop('HIGH_SR_BACKGROUND')

        # Pop members field
        if 'MEMBERS' in record.info.keys():
            record.info.pop('MEMBERS')

        # Write record to file
        fout.write(record)


def main():

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('fout')

    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    # No longer needed -- FILTERs assigned in cleanVCF step 5

    # #Set eligible filters
    # VCF_FILTERS = [
    # '##FILTER=<ID=MULTIALLELIC,Description="Multiallelic site">',
    # '##FILTER=<ID=UNRESOLVED,Description="Variant is incompletely resolved">',
    # '##FILTER=<ID=HIGH_SR_BACKGROUND,Description="An abundance of split reads in reference samples indicates a noisy locus">'
    # ]

    # # Add metadata lines for annotations
    header = vcf.header
    # for f in VCF_FILTERS:
    #     header.add_line(f)

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=header)

    clean_filt(vcf, fout)


if __name__ == '__main__':
    main()
