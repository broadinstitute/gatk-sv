#!/usr/bin/env python3

# Converts single sample SV/CNV vcf.gz (e.g., DRAGEN) to .del.bed and .dup.bed files

import argparse
from pysam import VariantFile

DEFAULT_QUAL_CUTOFF = 30


# Writes record in bed format
def write_bed_entry(file, record, score, sv_type, sample):
    bed_entries = (
        record.chrom,
        str(record.start),
        str(record.stop),
        str(score),
        sample,
        sv_type,
        'dragen'
    )
    file.write('\t'.join(bed_entries) + '\n')


# Writes record in bed format if it is a del or dup of sufficient quality
def convert_record(record, del_file, dup_file, input_sample, qual_cutoff, output_sample):
    # Extract Quality Score directly from the VCF QUAL column
    qual = record.qual

    # If QUAL is missing or below the cutoff, skip
    if qual is None or qual < qual_cutoff:
        return

    # pysam returns None for alts if it's a reference block ('.'), so we handle that safely
    alt_allele = record.alts[0] if record.alts else ""

    # Check ALT tags directly emitted by DRAGEN (<DEL> / <DUP>)
    if '<DEL>' in alt_allele:
        write_bed_entry(del_file, record, qual, 'DEL', output_sample)
    elif '<DUP>' in alt_allele:
        write_bed_entry(dup_file, record, qual, 'DUP', output_sample)


# Main function
def main():
    parser = argparse.ArgumentParser(
        description="Converts DRAGEN CNV/SV vcf.gz to .del.bed and .dup.bed files")
    parser.add_argument(
        "--input_sample", help="sample id in the vcf (defaults to first sample in vcf header)")
    parser.add_argument(
        "--cutoff", help="QUAL score cutoff", type=int, default=DEFAULT_QUAL_CUTOFF)
    parser.add_argument(
        "segments_vcf", help="Input .vcf or .vcf.gz file (e.g., DRAGEN output)")
    parser.add_argument(
        "output_sample", help="sample id to use in output and base name for output files")

    args = parser.parse_args()

    input_sample = args.input_sample
    qual_cutoff = args.cutoff

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
                qual_cutoff=qual_cutoff,
                output_sample=args.output_sample
            )


if __name__ == "__main__":
    main()
