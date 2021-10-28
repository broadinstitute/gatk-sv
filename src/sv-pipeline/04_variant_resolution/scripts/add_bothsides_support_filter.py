#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Label VCF with BOTHSIDES_SUPPORT based on input list
"""

import argparse
import sys
import pysam
from os import path
import csv
import subprocess


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
    parser.add_argument('-o', '--outfile',
                        help='Output file [default: stdout]')
    parser.add_argument('-z', '--bgzip', action='store_true',
                        help='Compress output file')
    parser.add_argument('vcf')
    parser.add_argument('vid_list')

    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    # Set eligible filters
    VCF_FILTERS = [
        '##FILTER=<ID=BOTHSIDES_SUPPORT,Description="Variant has read-level support for both sides of breakpoint">'
    ]

    # # Add metadata lines for annotations
    header = vcf.header
    for f in VCF_FILTERS:
        header.add_line(f)

    # Read list of IDs to label
    ids = read_vid_list(args.vid_list)

    # Open connection to outfile
    if args.outfile is None:
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        out = args.outfile
        if '.gz' in out or '.bgz' in out:
            out = path.splitext(out)[0]
        fout = pysam.VariantFile(out, 'w', header=header)

    # Process VCF
    for record in vcf:
        if record.id in ids:
            record.filter.add('BOTHSIDES_SUPPORT')
        fout.write(record)

    fout.close()  # not closing here can cause truncated output files if used with bzip below

    # Bgzip, if optioned
    if args.bgzip \
            and args.outfile is not None:
        subprocess.run(['bgzip', '-f', out])


if __name__ == '__main__':
    main()
