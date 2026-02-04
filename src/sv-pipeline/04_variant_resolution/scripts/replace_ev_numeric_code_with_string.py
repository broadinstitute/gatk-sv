#!/usr/bin/env python

import argparse
import gzip
import sys

ev_numeric_code_to_string_map = {
    '1': 'RD',
    '2': 'PE',
    '3': 'RD,PE',
    '4': 'SR',
    '5': 'RD,SR',
    '6': 'PE,SR',
    '7': 'RD,PE,SR'
}


def replace_ev_numeric_codes(vcf, fout):
    while True:
        line = vcf.readline().rstrip()
        if not line:
            break
        line_fields = line.split("\t")
        format_def = line_fields[8]
        format_def_fields = format_def.split(":")
        if 'EV' in format_def_fields:
            ev_idx = format_def_fields.index('EV')
            for idx in range(9, len(line_fields)):
                samp_gt_fields = line_fields[idx].split(":")
                ev_numeric = samp_gt_fields[ev_idx]
                ev_string = ev_numeric_code_to_string_map[ev_numeric]
                samp_gt_fields[ev_idx] = ev_string
                new_samp_gt = ":".join(samp_gt_fields)
                line_fields[idx] = new_samp_gt
        new_line = "\t".join(line_fields)
        print(new_line, file=fout)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('fout')

    args = parser.parse_args()

    # we don't use pysam for this due to bugs in our version: https://github.com/pysam-developers/pysam/issues/554
    # if we ever upgrade pysam to version 0.15.3 or later we should consider replacing this with a pysam implementation

    if args.vcf in '- stdin'.split():
        vcf = sys.stdin
    else:
        if args.vcf.endswith(".gz"):
            vcf = gzip.open(args.vcf, 'rt')
        else:
            vcf = open(args.vcf, 'r')

    new_ev_header_line = '##FORMAT=<ID=EV,Number=1,Type=String,Description="Classes of evidence supporting final genotype">'

    if args.fout in '- stdout'.split():
        fout = sys.stdout
    else:
        if args.fout.endswith(".gz"):
            fout = gzip.open(args.fout, 'wt')
        else:
            fout = open(args.fout, 'w')

    while True:
        line = vcf.readline().rstrip()
        if not line.startswith("##"):
            break
        if line.startswith('##FORMAT=<ID=EV,'):
            print(new_ev_header_line, file=fout)
        else:
            print(line, file=fout)

    # print the #CHROM line
    print(line, file=fout)

    replace_ev_numeric_codes(vcf, fout)


if __name__ == '__main__':
    main()
