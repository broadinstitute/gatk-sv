#!/bin/env python

import argparse
import pysam
import sys
from typing import List, Text, Optional


INFO_FIELDS = ["RU", "PERIOD", "LOCUS"]
FORMAT_FIELDS = ["REPCN", "REPCI1", "REPCI2", "OUTLIER", "ZSCORE", "DP", "STR_FILTER"]


def update_header(header: pysam.VariantHeader) -> None:
    """
    Adds needed STRipy header lines. These should correspond to the FORMAT_FIELDS and INFO_FIELDS variables.
    """
    header.add_line('##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit in the reference orientation">')
    header.add_line('##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Length of the repeat unit">')
    header.add_line('##INFO=<ID=DISEASES,Number=.,Type=String,Description="Associated disease symbols for this STR locus (| separated)">')
    header.add_line('##INFO=<ID=LOCUS,Number=1,Type=String,Description="Gene/locus identifier from STRipy">')
    header.add_line('##FORMAT=<ID=REPCN,Number=2,Type=Float,Description="Number of repeat units spanned by each allele">')
    header.add_line('##FORMAT=<ID=REPCI1,Number=2,Type=Integer,Description="95% CI min,max on repeat counts of first allele">')
    header.add_line('##FORMAT=<ID=REPCI2,Number=2,Type=Integer,Description="95% CI min,max on repeat counts of second allele">')
    header.add_line('##FORMAT=<ID=OUTLIER,Number=2,Type=Integer,Description="Allelic population outlier flags (0/1) assigned by STRipy">')
    header.add_line('##FORMAT=<ID=ZSCORE,Number=2,Type=Float,Description="Allelic population Z-scores assigned by STRipy">')
    header.add_line('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth at STR site">')
    header.add_line('##FORMAT=<ID=STR_FILTER,Number=.,Type=String,Description="Filter status assigned by STRipy">')


def validate_sample_id(main_header: pysam.VariantHeader,
                       stripy_header: pysam.VariantHeader):
    if len(stripy_header.samples) != 1:
        raise ValueError(f"Expected exactly 1 sample in STRipy header but got {len(stripy_header.samples)}")
    sample_raw = stripy_header.samples[0]
    # TODO temp fix; we should expect the correct ID here
    if not sample_raw.endswith(".final"):
        raise ValueError(f"Expected STRipy sample ID to end with .final but got {sample_raw}")
    sample = sample_raw.rsplit(".", 1)[0]
    if sample not in main_header.samples:
        raise ValueError(f"Sample {sample} (raw ID {sample_raw}) not found in single-sample VCF header")
    return sample


def _process(sample: Text,
             main_vcf: pysam.VariantFile,
             stripy_vcf: pysam.VariantFile,
             out_vcf: pysam.VariantFile) -> None:
    """"
    Master function for processing the given input vcf and writing output
    """
    for record in main_vcf:
        out_vcf.write(record)
    for stripy_record in stripy_vcf:
        record = out_vcf.new_record(contig=stripy_record.contig,
                                    start=stripy_record.start,
                                    stop=stripy_record.stop, alleles=stripy_record.alleles,
                                    id=stripy_record.id, qual=stripy_record.qual, filter=stripy_record.filter)
        for key in INFO_FIELDS:
            if key in stripy_record.info:
                record.info[key] = stripy_record.info[key]
        record.samples[sample]['GT'] = (None, None)
        for key in FORMAT_FIELDS:
            s = stripy_record.samples[sample + ".final"]
            if key in s:
                record.samples[sample][key] = s[key]
        out_vcf.write(record)


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Merges a single-sample STRipy VCF into a GATK-SV call set",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--main-vcf", type=str, required=True,
                        help="GATK-SV vcf, must contain the sample provided by --stripy-vcf")
    parser.add_argument("--stripy-vcf", type=str, required=True,
                        help="Single-sample STRipy vcf, with sample ID ending with \".final\"")
    parser.add_argument("--out", type=str, required=True,
                        help="Output vcf (unsorted)")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = _parse_arguments(argv)

    # convert vcf header and records
    with pysam.VariantFile(arguments.main_vcf) as main_vcf, pysam.VariantFile(arguments.stripy_vcf) as str_vcf:
        sample = validate_sample_id(main_header=main_vcf.header, stripy_header=str_vcf.header)
        update_header(header=main_vcf.header)
        with pysam.VariantFile(arguments.out, mode='w', header=main_vcf.header) as out_vcf:
            _process(sample=sample, main_vcf=main_vcf, stripy_vcf=str_vcf, out_vcf=out_vcf)


if __name__ == "__main__":
    main()
