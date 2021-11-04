#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import pysam


def label_mosaics(vcf, mosaics):
    for record in vcf:
        if record.id in mosaics.keys():
            record.info['MOSAIC'] = mosaics[record.id]
        yield record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('mosaics', type=argparse.FileType('r'))
    parser.add_argument('fout')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    vcf.header.add_line(
        '##INFO=<ID=MOSAIC,Number=.,Type=String,Description="Samples predicted to harbor somatic or germline mosaicism">')
    fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    mosaics = {}
    for line in args.mosaics:
        name, sample, sep = line.strip().split()
        mosaics[name] = [sample]

    for record in label_mosaics(vcf, mosaics):
        fout.write(record)


if __name__ == '__main__':
    main()
