#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Identify and classify duplicated variants from an input VCF.
"""

from typing import List, Text, Optional
from collections import defaultdict
from itertools import groupby

import argparse
import sys
import pysam


def process_duplicates(vcf, fout):
    # Initialize counters and buffers
    counts = defaultdict(int)
    exact_buffer = []
    ins_buffer = []
    current_chrom = None
    current_pos = None

    # Create output files
    with open(f"{fout}_duplicate_records.tsv", 'w') as f_records, open(f"{fout}_duplicate_counts.tsv", 'w') as f_counts:
        f_records.write("TYPE\tDUP_RECORDS\n")
        f_counts.write("TYPE\tDUP_COUNTS\n")

        # Iterate through all records
        for record in vcf.fetch():
            # Process current buffers if we've reached a new chrom or pos
            if record.chrom != current_chrom or record.pos != current_pos:
                process_buffers(exact_buffer, ins_buffer, counts, f_records)
                exact_buffer = []
                ins_buffer = []
                current_chrom = record.chrom
                current_pos = record.pos

            # Update buffers with new record
            exact_key = (
                record.chrom,
                record.pos,
                record.stop,
                record.info.get('SVTYPE'),
                record.info.get('SVLEN'),
                record.info.get('CHR2', ""),
                record.info.get('END2', 0),
                record.info.get('STRANDS', ""),
                record.info.get('CPX_TYPE', ""),
                record.info.get('CPX_INTERVALS', "")
            )
            exact_buffer.append((exact_key, record.id))

            if record.info.get('SVTYPE') == 'INS':
                insert_key = (
                    record.id,
                    record.info.get('SVLEN'),
                    record.alts[0]
                )
                ins_buffer.append(insert_key)

        # Process remaining records in the buffer
        process_buffers(exact_buffer, ins_buffer, counts, f_records)

        # Write counts to file
        for match_type in sorted(counts.keys()):
            f_counts.write(f"{match_type}\t{counts[match_type]}\n")


def process_buffers(exact_buffer, ins_buffer, counts, f_records):
    # Process exact matches
    sorted_buffer = sorted(exact_buffer, key=lambda x: x[0])
    exact_matches = {
        key: [record for _, record in group] for key, group in groupby(sorted_buffer, key=lambda x: x[0])
    }
    for records in exact_matches.values():
        if len(records) > 1:
            counts['exact'] += 1
            f_records.write(f"exact\t{','.join(sorted(records))}\n")

    # Process insert matches
    for i in range(len(ins_buffer)):
        for j in range(i + 1, len(ins_buffer)):
            rec1, rec2 = ins_buffer[i], ins_buffer[j]

            # Size comparison
            if rec1[1] == rec2[1]:
                counts['ins_size_similarity_100'] += 1
                f_records.write(f"ins_size_similarity_100\t{rec1[0]},{rec2[0]}\n")
            elif abs(rec1[1] - rec2[1]) <= 0.5 * max(rec1[1], rec2[1]):
                counts['ins_size_similarity_50'] += 1
                f_records.write(f"ins_size_similarity_50\t{rec1[0]},{rec2[0]}\n")
            else:
                counts['ins_size_similarity_0'] += 1
                f_records.write(f"ins_size_similarity_0\t{rec1[0]},{rec2[0]}\n")

            # ALT comparison
            if rec1[2] == rec2[2]:
                counts['ins_alt_identical'] += 1
                f_records.write(f"ins_alt_identical\t{rec1[0]},{rec2[0]}\n")
            elif ('<INS>' in (rec1[2], rec2[2])) and ('<INS:' in (rec1[2] + rec2[2])):
                counts['ins_alt_same_subtype'] += 1
                f_records.write(f"ins_alt_same_subtype\t{rec1[0]},{rec2[0]}\n")
            elif rec1[2].startswith('<INS:') and rec2[2].startswith('<INS:') and rec1[2] != rec2[2]:
                counts['ins_alt_different_subtype'] += 1
                f_records.write(f"ins_alt_different_subtype\t{rec1[0]},{rec2[0]}\n")


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Identify duplicated records from a sorted input VCF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-v', '--vcf', required=True, help='Input VCF.')
    parser.add_argument('-f', '--fout', required=True, help='Output file name.')
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    args = _parse_arguments(argv)

    vcf = pysam.VariantFile(args.vcf)
    process_duplicates(vcf, args.fout)


if __name__ == '__main__':
    main()
