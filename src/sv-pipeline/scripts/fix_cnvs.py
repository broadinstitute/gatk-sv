#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Pegs all multiallelic CNV genotypes to '.' and rescales CNQ from 0-999 to 0-99
"""

import argparse
import math
import sys
import pysam


def fix_cnvs(vcf, fout):
    for record in vcf.fetch():
        if record.info['SVTYPE'] != 'CNV':
            continue
        for gt in record.samples.values():
            gt['GT'] = None
            gt['CNQ'] = int(math.floor(gt['CNQ'] / 10.))
        fout.write(record)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input vcf (supports "stdin").')
    parser.add_argument('fout', help='Output file (supports "stdout").')

    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    fix_cnvs(vcf, fout)

    fout.close()


if __name__ == '__main__':
    main()
