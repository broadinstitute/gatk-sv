#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2018 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import heapq
from collections import deque
import pysam
import svtk.utils as svu


# def records_match(record, other):
#     """Test if two records are same SV: check chromosome, position, stop, SVTYPE, and (if they exist) CHR2/END2 for multi-chromosomal events"""
#     return (record.chrom == other.chrom and
#             (('CHR2' not in record.info and 'CHR2' not in other.info) or ('CHR2' in record.info and 'CHR2' in other.info and record.info['CHR2'] == other.info['CHR2'])) and
#             record.pos == other.pos and
#             record.stop == other.stop and
#             (('END2' not in record.info and 'END2' not in other.info) or ('END2' in record.info and 'END2' in other.info and record.info['END2'] == other.info['END2'])) and
#             record.info['SVTYPE'] == other.info['SVTYPE'])


# def merge_key(record):
#     """Sort records by coordinate then svtype"""
#     return (record.pos, record.stop, record.info['SVTYPE'], record.id)


# def dedup_records(records):
#     """Take unique subset of records"""

#     records = sorted(records, key=merge_key)

#     curr_record = records[0]
#     for record in records[1:]:
#         if records_match(curr_record, record):
#             print("Match!: ", end="") # debugging
#             print(curr_record, end=", ")
#             print(record)
#             continue
#         else:
#             yield curr_record
#             curr_record = record

#     yield curr_record


class VariantRecordComparison:
    def __init__(self, record):
        self.pos = record.pos
        self.stop = record.stop
        self.svtype = record.info['SVTYPE']
        self.id = record.id
        self.chrom = record.chrom
        self.chr2 = record.info['CHR2'] if 'CHR2' in record.info else None
        self.end2 = record.info['END2'] if 'END2' in record.info else None


    def __lt__(self, other):
        if self.chrom == other.chrom:
            return (self.pos, self.stop, self.svtype, self.chr2, self.end2, self.id) < (other.pos, other.stop, other.svtype, other.chr2, other.end2, other.id)
            # return (self.pos, self.stop, self.svtype, self.id) < (other.pos, other.stop, other.svtype, other.id)
        else:
            return svu.is_smaller_chrom(self.chrom, other.chrom)

    def __gt__(self, other):
        if self.chrom == other.chrom:
            return (self.pos, self.stop, self.svtype, self.chr2, self.end2, self.id) > (other.pos, other.stop, other.svtype, other.chr2, other.end2, other.id)
            # return (self.pos, self.stop, self.svtype, self.id) < (other.pos, other.stop, other.svtype, other.id)
        else:
            return (not svu.is_smaller_chrom(self.chrom, other.chrom))

    def __le__(self, other):
        return (((self.chrom, self.pos, self.stop, self.chr2, self.end2, self.svtype) == (other.chrom, other.pos, other.stop, other.chr2, other.end2, other.svtype)) or
            svu.is_smaller_chrom(self.chrom, other.chrom) or
            (self.pos, self.stop, self.svtype, self.chr2, self.end2, self.id) < (other.pos, other.stop, other.svtype, other.chr2, other.end2, other.id))

    def __ge__(self, other):
        return (((self.chrom, self.pos, self.stop, self.chr2, self.end2, self.svtype) == (other.chrom, other.pos, other.stop, other.chr2, other.end2, other.svtype)) or
            (not svu.is_smaller_chrom(self.chrom, other.chrom)) or
            (self.pos, self.stop, self.svtype, self.chr2, self.end2, self.id) > (other.pos, other.stop, other.svtype, other.chr2, other.end2, other.id))

    def __eq__(self, other):
        return (self.chrom, self.pos, self.stop, self.chr2, self.end2, self.svtype) == (other.chrom, other.pos, other.stop, other.chr2, other.end2, other.svtype)

def merge_records(vcfs):
    """
    Take unique set of VCFs
    """

    merged_vcfs = heapq.merge(*vcfs, key=lambda r: VariantRecordComparison(r))

    record = next(merged_vcfs)
    curr_record = record

    for record in merged_vcfs:
        # if records_match(curr_record, record):
        if VariantRecordComparison(curr_record) == VariantRecordComparison(record):
            print("Match!: ", end="") # debugging
            print(curr_record, end=", ")
            print(record)
            continue
        else:
            yield curr_record
            curr_record = record

    yield curr_record


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

    # Copy base VCF
    args.fout.write(str(vcfs[0].header))
    n_samples = len(vcfs[0].header.samples)

    # Write out null records for dedupped variants
    for record in merge_records(vcfs):
        base = '\t'.join(str(record).split('\t')[:8])
        null_gts = '\t'.join(['0/0' for i in range(n_samples - 1)])
        args.fout.write(base + '\tGT\t0/1\t' + null_gts + '\n')

    args.fout.close()


if __name__ == '__main__':
    main()
