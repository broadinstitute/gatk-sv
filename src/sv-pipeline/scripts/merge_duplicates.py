#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Merge duplicate records and counts files across multiple TSVs.
"""

from typing import List, Text, Optional
from collections import defaultdict

import argparse
import sys
import csv


def merge_duplicates(record_files: List[str], count_files: List[str], fout: str):
    # Merge records
    with open(f"{fout}_duplicate_records.tsv", 'w', newline='') as out_records:
        writer = csv.writer(out_records, delimiter='\t')
        # Write header
        writer.writerow(['TYPE', 'DUP_RECORDS'])
        # Append records from each file
        for record_file in record_files:
            with open(record_file, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                next(reader)
                for row in reader:
                    writer.writerow(row)

    # Sum counts
    counts = defaultdict(int)
    for count_file in count_files:
        with open(count_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader)
            for row in reader:
                counts[row[0]] += int(row[1])

    # Merge counts
    with open(f"{fout}_duplicate_counts.tsv", 'w', newline='') as out_counts:
        writer = csv.writer(out_counts, delimiter='\t')

        # Write header
        writer.writerow(['TYPE', 'DUP_COUNTS'])

        # Append each row from merged counts
        for category, count in counts.items():
            writer.writerow([category, count])


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge duplicate records and counts files across multiple TSVs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-r', '--records', nargs='+', required=True, help='Input duplicated record TSV files.')
    parser.add_argument('-c', '--counts', nargs='+', required=True, help='Input duplicated counts TSV files.')
    parser.add_argument('-f', '--fout', required=True, help='Output file name.')
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    args = _parse_arguments(argv)

    merge_duplicates(args.records, args.counts, args.fout)


if __name__ == '__main__':
    main()
