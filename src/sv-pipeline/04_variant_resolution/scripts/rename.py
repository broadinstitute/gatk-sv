#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Rename variants & do final cleanup after complex resolution
"""

import argparse
import sys
from collections import defaultdict
import pysam


def rename(vcf, fout, chrom=None, prefix='SV_pipeline'):
    indexes = defaultdict(int)
    fmt = '{prefix}_{svtype}_{chrom}_{idx}'

    for record in vcf:
        if 0 <= record.info['SVLEN'] < 50:
            continue

        if chrom is None:
            chrom = record.chrom
        svtype = record.info['SVTYPE']
        indexes[svtype] += 1
        idx = indexes[svtype]
        prefix = prefix

        record.id = fmt.format(**locals())

        # Clean metadata
        record.ref = 'N'

        # Clean metadata for CNVs sucked into complex resolution
        # and not appropriately cleaned
        if svtype in 'DEL DUP'.split():
            for info in 'EVENT UNRESOLVED UNRESOLVED_TYPE'.split():
                if info in record.info.keys():
                    record.info.pop(info)

        # Assign all unresolved inversions as BNDs
        if svtype == 'INV':
            if 'UNRESOLVED' in record.info.keys():
                record.info['SVTYPE'] = 'BND'

        if svtype != 'BND':
            if 'STRANDS' in record.info.keys():
                record.info.pop('STRANDS')

        # Write record to file
        fout.write(record)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('fout')
    parser.add_argument('--chrom')
    parser.add_argument('--prefix', default='SV_pipeline',
                        help='Tag prepended to all variant IDs')

    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    rename(vcf, fout, args.chrom, args.prefix)


if __name__ == '__main__':
    main()
