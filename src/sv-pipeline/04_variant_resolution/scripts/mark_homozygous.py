#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import sys
import pysam


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('homozygotes', type=argparse.FileType('r'))
    parser.add_argument('fout')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    homozygotes = [l.strip().split() for l in args.homozygotes.readlines()]
    homozygotes = {name: samples.split(',') for name, samples in homozygotes}

    for record in vcf:
        if record.id in homozygotes:
            for sample in homozygotes[record.id]:
                record.samples[sample]['GT'] = (1, 1)
        fout.write(record)


if __name__ == '__main__':
    main()
