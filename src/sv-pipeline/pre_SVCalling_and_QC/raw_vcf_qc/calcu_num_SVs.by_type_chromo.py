#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
"""

import os
import argparse
import gzip


def chr_extract(pin, chr):
    out = ''
    for i in pin[7].split(';'):
        if i.split('=')[0] == chr:
            out = i.split('=')[1]
    return out


def vcf_stat_readin(filename):
    fin = gzip.open(filename, 'rb')
    out = {}
    sv_types = []
    for line in fin:
        pin = line.strip().decode().split()
        if not pin[0][0] == '#':
            if not pin[0] in out.keys():
                out[pin[0]] = {}
            svtype = chr_extract(pin, 'SVTYPE')
            if svtype not in out[pin[0]].keys():
                out[pin[0]][svtype] = 0
                if svtype not in sv_types:
                    sv_types.append(svtype)
            out[pin[0]][svtype] += 1
    fin.close()
    return [out, sv_types]


def write_output(output_name, stat_hash, sv_types, sample_name):
    if not os.path.isfile(output_name):
        fo = open(output_name, 'w')
        print('\t'.join(['#CHROM', 'SVTYPE', 'NUM', 'SAMPLE']), file=fo)
        fo.close()
    fo = open(output_name, 'a')
    for i in sv_types:
        for j in stat_hash.keys():
            if i in stat_hash[j].keys():
                print(
                    '\t'.join([str(k) for k in [j, i, stat_hash[j][i], sample_name]]), file=fo)
            else:
                print('\t'.join([str(k)
                                 for k in [j, i, 0, sample_name]]), file=fo)
    fo.close()


def main():
    parser = argparse.ArgumentParser("calcu_num_SVs.by_type_chromo.py")
    parser.add_argument(
        "vcf", type=str, help="name of input file in vcf.gz format")
    parser.add_argument("out", type=str, help="name of output file")
    args = parser.parse_args()
    sample_name = args.vcf.replace('.vcf.gz', '')
    [stat_hash, sv_types] = vcf_stat_readin(args.vcf)
    write_output(args.out, stat_hash, sv_types, sample_name)


if __name__ == '__main__':
    main()
