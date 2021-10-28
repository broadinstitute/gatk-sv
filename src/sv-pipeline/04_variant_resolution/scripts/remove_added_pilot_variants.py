#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Remove calls from Pilot pass VCF that were added in the raw overlap with Phase1
"""

import argparse
import pysam


def records_match(A, B):
    return (A.chrom == B.chrom and
            A.pos == B.pos and
            A.stop == B.stop and
            A.info['STRANDS'] == B.info['STRANDS'])


def filter_vcf(recordsA, recordsB):
    """Remove records in B from A."""

    currB = next(recordsB)

    for record in recordsA:
        # While A is lagging, all variants are OK
        if record.pos < currB.pos:
            yield record
        else:
            # Skip until B matches or is ahead
            while record.pos > currB.pos:
                currB = next(recordsB)
            # Skip if the record is in B
            if records_match(record, currB):
                continue
            else:
                yield record


def remove_overlap(overlap, raw_pilot, pass_pilot):
    # First, get IDs of raw variants from overlap VCF
    raw_IDs = [m for record in overlap for m in record.info['MEMBERS']]
    raw_IDs = [m for m in raw_IDs if m.startswith('Pilot')]

    # Next, get corresponding raw variant records
    raw_records = (record for record in raw_pilot if record.id in raw_IDs)

    # Finally, if a record in the passing VCF matches, skip it
    for record in filter_vcf(pass_pilot, raw_records):
        yield record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('overlap_vcf', help="Output of overlap_pass")
    parser.add_argument('raw_pilot_vcf', help="Raw vcfclustered Pilot vcf.")
    parser.add_argument('pass_pilot_vcf', help="De novo filtered, passing "
                        "pilot VCF.")
    parser.add_argument('fout')
    args = parser.parse_args()

    overlap_vcf = pysam.VariantFile(args.overlap_vcf)
    raw_pilot_vcf = pysam.VariantFile(args.raw_pilot_vcf)
    pass_pilot_vcf = pysam.VariantFile(args.pass_pilot_vcf)

    fout = pysam.VariantFile(args.fout, 'w', header=pass_pilot_vcf.header)

    for record in remove_overlap(overlap_vcf, raw_pilot_vcf, pass_pilot_vcf):
        fout.write(record)


if __name__ == '__main__':
    main()
