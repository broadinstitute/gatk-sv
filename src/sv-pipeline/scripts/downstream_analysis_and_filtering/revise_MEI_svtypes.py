#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 The Broad Institute of M.I.T. and Harvard
# Distributed under terms of the MIT license.
# Contact: Ryan Collins <rlcollins@g.harvard.edu>

"""
Revise SVTYPE fields of INS:ME from INS to MEI (or vice-versa)
"""


import argparse
import pysam
import sys


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('invcf', help='Input vcf. Also accepts "stdin" and "-".')
    parser.add_argument('outvcf', help='Output vcf. Also accepts "stdout" and "-".')
    parser.add_argument('--reverse', help='Instead of assigning SVTYPE=MEI, assign ' +
                        'SVTYPE=INS', action='store_true', default=False)
    args = parser.parse_args()

    # Open connection to input VCF
    if args.invcf in '- stdin'.split():
        invcf = pysam.VariantFile(sys.stdin) 
    else:
        invcf = pysam.VariantFile(args.invcf)

    # Open connections to output vcf
    if args.outvcf in '- stdin'.split():
        outvcf = pysam.VariantFile(sys.stdout, 'w', header=invcf.header) 
    else:
        outvcf = pysam.VariantFile(args.outvcf, 'w', header=invcf.header)

    # Iterate over records in input VCF and modify as needed before writing to output vcf
    for record in invcf:
        if any(['INS:ME' in a for a in record.alts]):
            if args.reverse:
                new_svtype = 'INS'
            else:
                new_svtype = 'MEI'
            record.info['SVTYPE'] = new_svtype
        # Add SVLEN if missing
        if 'SVLEN' not in record.info.keys():
            record.info['SVLEN'] = 0
        outvcf.write(record)

    # Close output VCF file handle
    outvcf.close()


if __name__ == '__main__':
    main()

