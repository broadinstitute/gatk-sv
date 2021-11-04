#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Drop records with zero non-reference alleles in AF-annotated VCF
"""

import argparse
import sys
import pysam
from svtk.utils import is_biallelic


NULL_GTs = [(0, 0), (None, None), (0, ), (None, ), (None, 2)]


def prune_allrefs(vcf, fout):
    for record in vcf.fetch():
        if is_biallelic(record):
            if record.info['AC'][0] > 0:
                fout.write(record)
        else:
            if record.info['CN_NONREF_FREQ'] > 0:
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

    prune_allrefs(vcf, fout)

    fout.close()


if __name__ == '__main__':
    main()
