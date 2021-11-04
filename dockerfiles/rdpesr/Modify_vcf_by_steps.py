#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
"""

import pysam
import argparse


def modify_vcf(vcf_in_file, vcf_out_file, step_size, contig):
    vcf_in = pysam.VariantFile(vcf_in_file)
    vcf_out = pysam.VariantFile(vcf_out_file, 'w', header=vcf_in.header)
    for rec in vcf_in:
        rec.pos += step_size
        rec.stop += step_size
        if rec.pos > 0 and not rec.stop > contig[rec.contig]:
            vcf_out.write(rec)
    vcf_out.close()
    vcf_in.close()


def contig_readin(contig):
    out = {}
    fin = open(contig)
    for line in fin:
        pin = line.strip().split()
        out[pin[0]] = int(pin[1])
    fin.close()
    return out


def main():
    parser = argparse.ArgumentParser(
        description='Shift each variants in vcf by a fixed step.')
    parser.add_argument('vcf_in', metavar='', type=str,
                        help='name of vcf file to be modified')
    parser.add_argument('vcf_out', metavar='', type=str,
                        help='name of output vcf file')
    parser.add_argument('-s', '--step_size', type=int,
                        help='size of step to be shifted.')
    parser.add_argument('-c', '--contig', type=str,
                        help='contig files, or reference index.')

    args = parser.parse_args()
    contig = contig_readin(args.contig)
    modify_vcf(args.vcf_in, args.vcf_out, args.step_size, contig)


if __name__ == "__main__":
    main()
