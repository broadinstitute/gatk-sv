#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Replaces EV field numeric codes with strings and fixes alleles for multi-allelic CNVs
"""

import argparse
import pysam
import sys

EVIDENCE_LIST = [
    '',           # 0
    'RD',         # 1
    'PE',         # 2
    'RD,PE',      # 3
    'SR',         # 4
    'RD,SR',      # 5
    'PE,SR',      # 6
    'RD,PE,SR'    # 7
]


def get_new_header(header):
    header_list = str(header).split('\n')
    new_header_lines = list()
    for line in header_list:
        if line.startswith('##FORMAT=<ID=EV,'):
            # Updates Number/Type
            new_header_lines.append('##FORMAT=<ID=EV,Number=.,Type=String,Description="Classes of evidence supporting final genotype">')
        elif line.startswith('##INFO=<ID=MULTIALLELIC,'):
            # Remove MULTIALLELIC field (legacy)
            continue
        elif line.startswith('#CHROM'):
            # Exclude samples line
            continue
        elif not line:
            # Skip empty lines
            continue
        else:
            new_header_lines.append(line)
    new_header = pysam.VariantHeader()
    new_header.add_samples(header.samples)
    for line in new_header_lines:
        print(line)
        new_header.add_line(line)
    return new_header


def process_record(record):
    # Fix multi-allelic CNV alts (legacy)
    if record.alts[0].startswith('<CN'):
        record.alts = ('<CNV>',)
    # Remove MULTIALLELIC field (legacy)
    record.info.pop('MULTIALLELIC')
    # Since pysam messes with some of the formatting (i.e. END limitation) we parse the string and replace
    max_ev = len(EVIDENCE_LIST) - 1
    record_tokens = str(record).strip().split('\t')
    format_keys = record_tokens[8].split(':')
    ev_index = format_keys.index('EV')
    new_record_tokens = record_tokens[:9]
    for gt in record_tokens[9:]:
        gt_tokens = gt.split(':')
        ev = gt.split(':')[ev_index]
        if ev == '.' or not ev:
            new_ev = '.'
        else:
            ev_int = int(ev)
            if ev_int < 0 or ev_int > max_ev:
                raise ValueError(f"Invalid EV value {ev_int} in record {record.id}")
            new_ev = EVIDENCE_LIST[ev_int]
        gt_tokens[ev_index] = new_ev
        new_record_tokens.append(':'.join(gt_tokens))
    return '\t'.join(new_record_tokens)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    args = parser.parse_args()

    with pysam.VariantFile(args.vcf) as vcf:
        new_header = get_new_header(vcf.header)
        sys.stdout.write(str(new_header))
        for record in vcf:
            new_record = process_record(record)
            sys.stdout.write(new_record + "\n")


if __name__ == '__main__':
    main()
