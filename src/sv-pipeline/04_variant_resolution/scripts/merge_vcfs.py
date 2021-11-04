#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
"""

import argparse
import heapq
from collections import deque
import pysam
import svtk.utils as svu


def records_match(record, other):
    """
    Test if two records are same SV: check chromosome, position, stop, SVTYPE, SVLEN (for insertions),
    STRANDS (for BNDS and INVs), and (if they exist) CHR2/END2 for multi-chromosomal events
    """
    return (record.chrom == other.chrom and
            record.pos == other.pos and
            record.stop == other.stop and
            record.info['SVTYPE'] == other.info['SVTYPE'] and
            record.info['SVLEN'] == other.info['SVLEN'] and
            record.info['STRANDS'] == other.info['STRANDS'] and
            (('CHR2' not in record.info and 'CHR2' not in other.info) or ('CHR2' in record.info and 'CHR2' in other.info and record.info['CHR2'] == other.info['CHR2'])) and
            (('END2' not in record.info and 'END2' not in other.info) or ('END2' in record.info and 'END2' in other.info and record.info['END2'] == other.info['END2'])))


def merge_key(record):
    """
    Sort records by all fields that records_match will use to check for duplicates, in sequence,
    so that all identical records according to records_match will be adjacent
    """
    chr2 = record.info['CHR2'] if 'CHR2' in record.info else None
    end2 = record.info['END2'] if 'END2' in record.info else None
    return (record.pos, record.stop, record.info['SVTYPE'], record.info['SVLEN'], chr2, end2, record.info['STRANDS'], record.id)


def dedup_records(records):
    """Take unique subset of records"""

    records = sorted(records, key=merge_key)

    curr_record = records[0]
    for record in records[1:]:
        if records_match(curr_record, record):
            # keep more informative ALT field, assumed to be the one with more colons
            # ex: <INS:ME:ALU> kept over <INS>
            curr_alt = curr_record.alts[0]
            new_alt = record.alts[0]
            if (curr_alt.startswith('<') and curr_alt.endswith('>') and new_alt.startswith('<') and new_alt.endswith('>') and
                    len(new_alt.split(':')) > len(curr_alt.split(':'))):
                curr_record = record
            continue
        else:
            yield curr_record
            curr_record = record

    yield curr_record


class VariantRecordComparison:
    def __init__(self, record):
        self.record = record

    def __lt__(self, other):
        if self.record.chrom == other.record.chrom:
            return self.record.pos < other.record.pos
        else:
            return svu.is_smaller_chrom(self.record.chrom, other.record.chrom)


def merge_records(vcfs):
    """
    Take unique set of VCF records
    Strategy: Merge & roughly sort records from all VCFs by chrom & pos, then gather records that share the same chrom & pos and remove duplicates.
    Note: The output from heapq.merge cannot be directly used to remove duplicates because it is not sufficiently sorted, so duplicates may not be
        adjacent. It is also not sufficient to alter the comparator function to take more than chrom & pos into account, because heapq.merge assumes
        that each VCF is already sorted and will make no attempt to further sort them according to the comparator function. Re-sorting all records
        that share a chrom & pos by all necessary comparison fields is more efficient than re-sorting each entire VCF.
    """

    merged_vcfs = heapq.merge(*vcfs, key=lambda r: VariantRecordComparison(r))

    record = next(merged_vcfs)
    curr_records = deque([record])
    curr_chrom = record.chrom
    curr_pos = record.pos

    for record in merged_vcfs:
        if record.chrom == curr_chrom and record.pos == curr_pos:
            curr_records.append(record)
        else:
            for rec in dedup_records(curr_records):
                yield rec

            curr_records = deque([record])
            curr_pos = record.pos
            curr_chrom = record.chrom

    for rec in dedup_records(curr_records):
        yield rec


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcflist', type=argparse.FileType('r'))
    parser.add_argument('fout', type=argparse.FileType('w'))
    args = parser.parse_args()

    # VCFs from other batches
    fnames = [l.strip() for l in args.vcflist.readlines()]
    vcfs = [pysam.VariantFile(f) for f in fnames]

    # Copy base VCF header without samples
    args.fout.write('\t'.join(str(vcfs[0].header).split('\t')[:11]) + '\n')

    # Write out sites-only records for dedupped variants + 2 dummy GTs
    # including one 0/1 so svtk bedcluster doesn't break & clusters only on variants not samples
    for record in merge_records(vcfs):
        base = '\t'.join(str(record).split('\t')[:8])
        args.fout.write(base + '\tGT\t0/1\t0/0\n')

    args.fout.close()


if __name__ == '__main__':
    main()
