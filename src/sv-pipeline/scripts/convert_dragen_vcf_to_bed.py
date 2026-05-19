#!/usr/bin/env python3

# Converts single sample SV/CNV vcf.gz (e.g., DRAGEN) to .del.bed and .dup.bed files
# Uses both explicit SVTYPE/ALT tags AND Copy Number (CN) fields to infer DEL/DUP

import sys
import argparse
from pysam import VariantFile

COPY_NUMBER_FIELD = 'CN'
SCORE_FIELD = 'QS'
DEFAULT_CUTOFF = 30


# Writes record in bed format
def write_bed_entry(file, record, score, sv_type, sample):
    bed_entries = (
        record.chrom,
        str(record.start),
        str(record.stop),
        str(score),
        sample,
        sv_type,
        'gcnv'
    )
    file.write('\t'.join(bed_entries) + '\n')


# Writes record in bed format if it is a del or dup of sufficient quality
def convert_record(record, del_file, dup_file, input_sample, qs_cutoff, output_sample, xy_ploidy):
    sample_data = record.samples[input_sample]

    # Extract Quality Score
    QS = sample_data.get(SCORE_FIELD)
    if QS is None:
        QS = record.qual

    if isinstance(QS, tuple):
        QS = QS[0]

    if QS is not None and QS >= qs_cutoff:
        alt_allele = record.alts[0] if record.alts else ""
        svtype = record.info.get('SVTYPE', '')

        # 1. Try explicit SVTYPE or ALT string (Standard for SV callers)
        if svtype == 'DEL' or '<DEL>' in alt_allele:
            write_bed_entry(del_file, record, QS, 'DEL', output_sample)
            return
        elif svtype == 'DUP' or '<DUP>' in alt_allele:
            write_bed_entry(dup_file, record, QS, 'DUP', output_sample)
            return

        # 2. Fallback to Copy Number (CN) field (Standard for DRAGEN CNV callers)
        CN = sample_data.get(COPY_NUMBER_FIELD)
        if CN is not None:
            if isinstance(CN, tuple):
                CN = CN[0]

            # Determine baseline ploidy
            chrom = record.chrom.lower()
            if 'x' in chrom or 'y' in chrom:
                baseline_ploidy = xy_ploidy
            else:
                baseline_ploidy = 2  # Autosomes default to diploid

            # Separate based on copy number vs baseline
            if CN < baseline_ploidy:
                write_bed_entry(del_file, record, QS, 'DEL', output_sample)
            elif CN > baseline_ploidy:
                write_bed_entry(dup_file, record, QS, 'DUP', output_sample)


# Main function
def main():
    parser = argparse.ArgumentParser(
        description="Converts CNV/SV vcf.gz to .del.bed and .dup.bed files")
    parser.add_argument(
        "--input_sample", help="sample id in the vcf (defaults to first sample in vcf header)")
    parser.add_argument(
        "--cutoff", help="QS score cutoff", type=int, default=DEFAULT_CUTOFF)
    parser.add_argument(
        "--xy_ploidy", help="Baseline ploidy for sex chromosomes (1 for XY males, 2 for XX females). Defaults to 2.",
        type=int, default=2)
    parser.add_argument(
        "segments_vcf", help="Input .vcf or .vcf.gz file (e.g., DRAGEN output)")
    parser.add_argument(
        "output_sample", help="sample id to use in output and base name for output files")

    args = parser.parse_args()

    input_sample = args.input_sample
    qs_cutoff = args.cutoff

    del_name = args.output_sample + ".del.bed"
    dup_name = args.output_sample + ".dup.bed"

    with open(del_name, 'w') as del_bed, open(dup_name, 'w') as dup_bed:
        vcf = VariantFile(args.segments_vcf)

        if not input_sample:
            input_sample = list(vcf.header.samples)[0]

        for rec in vcf.fetch():
            convert_record(
                record=rec,
                del_file=del_bed,
                dup_file=dup_bed,
                input_sample=input_sample,
                qs_cutoff=qs_cutoff,
                output_sample=args.output_sample,
                xy_ploidy=args.xy_ploidy
            )


if __name__ == "__main__":
    main()