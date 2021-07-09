#!/bin/python

# Converts gCNV segments vcf to .del.bed and .dup.bed files
# Author: Mark Walker (markw@broadinstitute.org)

import sys
import argparse
from pysam import VariantFile

COPY_NUMBER_FIELD = 'CN'
SCORE_FIELD = 'QS'
DEFAULT_CUTOFF = 30

# Writes record in bed format


def write_bed_entry(file, record, score, type, sample):
    bed_entries = (record.chrom, str(record.start), str(
        record.stop), str(score), sample, type, 'gcnv')
    file.write('\t'.join(bed_entries) + '\n')

# Writes record in bed format if it is a del or dup of sufficient quality


def convert_record(record, del_file, dup_file, input_sample, cutoff, neutral_copy_state, output_sample):
    CN = record.samples[input_sample][COPY_NUMBER_FIELD]
    if CN != neutral_copy_state:
        QS = record.samples[input_sample][SCORE_FIELD]
        if QS >= cutoff:
            if CN < neutral_copy_state:
                write_bed_entry(del_file, record, QS, 'DEL', output_sample)
            else:
                write_bed_entry(dup_file, record, QS, 'DUP', output_sample)

# Gets dictionary of contig ploidy from gcnv ploidy call tsv file


def get_contig_ploidy(file):
    ploidy = {}
    for line in file:
        if line.startswith("@") or line.startswith("CONTIG\t"):
            continue
        tokens = line.strip().split('\t')
        ploidy[tokens[0]] = int(tokens[1])
    return ploidy

# Main function


def main():
    parser = argparse.ArgumentParser(
        description="Converts gCNV segments vcf to .del.bed and .dup.bed files")
    parser.add_argument(
        "--input_sample", help="sample id in the vcf (defaults to first sample in vcf header)")
    parser.add_argument("--cutoff", help="QS score cutoff")
    parser.add_argument(
        "ploidy_calls", help="gCNV ploidy calls file (contig_ploidy.tsv)")
    parser.add_argument("segments_vcf", help="gCNV segments vcf")
    parser.add_argument(
        "output_sample", help="sample id to use in output and base name for output files")
    args = parser.parse_args()

    input_sample = args.input_sample
    if args.cutoff:
        cutoff = int(args.cutoff)
    else:
        cutoff = DEFAULT_CUTOFF

    del_name = args.output_sample + ".del.bed"
    dup_name = args.output_sample + ".dup.bed"

    with open(args.ploidy_calls, 'r') as f:
        contig_ploidy = get_contig_ploidy(f)

    with open(del_name, 'w') as del_bed, open(dup_name, 'w') as dup_bed:
        vcf = VariantFile(args.segments_vcf)
        if not input_sample:
            input_sample = vcf.header.samples[0]
        for rec in vcf.fetch():
            if rec.chrom in contig_ploidy:
                record_ploidy = contig_ploidy[rec.chrom]
                convert_record(record=rec, del_file=del_bed, dup_file=dup_bed, input_sample=input_sample, cutoff=cutoff,
                               neutral_copy_state=record_ploidy, output_sample=args.output_sample)
            else:
                sys.stderr.write(
                    "Warning: no ploidy call found for contig: \"" + rec.chrom + "\", record omitted.\n")


if __name__ == "__main__":
    main()
