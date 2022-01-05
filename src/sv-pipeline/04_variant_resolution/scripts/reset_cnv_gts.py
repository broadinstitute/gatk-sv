#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Sets CNV GT fields to "."
This is needed following any HailMerge step for VCFs containing CNVs
"""

import argparse
import sys
import pysam


def reset_cnv_gts(vcf, fout):

    for record in vcf:
        if record.info['SVTYPE'] == 'CNV':
            for sample in record.samples:
                record.samples[sample]['GT'] = (None,)
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

    header = vcf.header

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=header)

    reset_cnv_gts(vcf, fout)


if __name__ == '__main__':
    main()
