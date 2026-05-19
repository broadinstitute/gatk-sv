#!/usr/bin/env python3

# Converts single sample SV/CNV vcf.gz (e.g., DRAGEN) to .del.bed and .dup.bed files

import argparse
from pysam import VariantFile

COPY_NUMBER_FIELD = 'CN'
DEFAULT_QUAL_CUTOFF = 30
DEFAULT_MIN_SVLEN = 5000


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


# Writes record in bed format if it is a del or dup of sufficient quality and length
def convert_record(record, del_file, dup_file, input_sample, qual_cutoff, min_svlen, output_sample, xy_ploidy):
    # Extract Quality Score directly from the VCF QUAL column
    qual = record.qual

    # If QUAL is missing or below the cutoff, skip
    if qual is None or qual < qual_cutoff:
        return

    # Extract SVLEN from the INFO field.
    # pysam returns INFO fields as tuples (e.g., (-1036,)), so we grab the first element.
    svlen = record.info.get('SVLEN')
    if svlen is not None:
        if isinstance(svlen, tuple):
            svlen = svlen[0]
        # Use absolute value because deletions have negative lengths (e.g., -1036)
        if abs(svlen) < min_svlen:
            return
    else:
        # If no SVLEN is present, it's likely a reference (REF) block. Skip.
        return

    sample_data = record.samples[input_sample]
    alt_allele = record.alts[0] if record.alts else ""

    # Check ALT tags first
    if '<DEL>' in alt_allele:
        write_bed_entry(del_file, record, qual, 'DEL', output_sample)
        return
    elif '<DUP>' in alt_allele:
        write_bed_entry(dup_file, record, qual, 'DUP', output_sample)
        return

    # Fallback to Copy Number (CN) if ALT tag is generic (e.g., <CNV>)
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
            write_bed_entry(del_file, record, qual, 'DEL', output_sample)
        elif CN > baseline_ploidy:
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
        "--min_svlen", help="Minimum absolute length of SV to keep", type=int, default=DEFAULT_MIN_SVLEN)
    parser.add_argument(
        "--xy_ploidy", help="Baseline ploidy for sex chromosomes (1 for XY males, 2 for XX females). Defaults to 2.",
        type=int, default=2)
    parser.add_argument(
        "segments_vcf", help="Input .vcf or .vcf.gz file (e.g., DRAGEN output)")
    parser.add_argument(
        "output_sample", help="sample id to use in output and base name for output files")

    args = parser.parse_args()

    input_sample = args.input_sample
    qual_cutoff = args.cutoff
    min_svlen = args.min_svlen

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
                min_svlen=min_svlen,
                output_sample=args.output_sample,
                xy_ploidy=args.xy_ploidy
            )


if __name__ == "__main__":
    main()