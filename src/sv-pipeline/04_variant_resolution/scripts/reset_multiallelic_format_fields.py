#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import csv
import os
import pysam
import sys


def read_vid_list(vid_list):
    ids = []
    with open(vid_list) as fin:
        reader = csv.reader(fin, delimiter='\t')
        for row in reader:
            vid = row[-1]
            if vid not in ids:
                ids.append(vid)

    return ids


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('vid_list')
    parser.add_argument('-o', '--outfile',
                        help='Output file [default: stdout]')

    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)

    # Read list of IDs to label
    ids = read_vid_list(args.vid_list)

    # Set eligible filters
    VCF_FORMAT_LINES = [
        '##FORMAT=<ID=CNQ,Number=1,Type=Integer,Description="Read-depth genotype quality">',
        '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Predicted copy state">'
    ]

    # # Add metadata lines for annotations
    header = vcf.header
    for f in VCF_FORMAT_LINES:
        header.add_line(f)

    # Open connection to outfile
    if args.outfile is None:
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        out = args.outfile
        if '.gz' in out or '.bgz' in out:
            out = os.path.splitext(out)[0]
        fout = pysam.VariantFile(out, 'w', header=header)

    for record in vcf:
        if record.id in ids:
            for j, sample in enumerate(record.samples):
                record.samples[sample]['GT'] = None
                record.samples[sample]['GQ'] = None
                record.samples[sample]['CN'] = record.samples[sample]['RD_CN']
                record.samples[sample]['CNQ'] = record.samples[sample]['RD_GQ']
            fout.write(record)

    fout.close()


if __name__ == '__main__':
    main()
