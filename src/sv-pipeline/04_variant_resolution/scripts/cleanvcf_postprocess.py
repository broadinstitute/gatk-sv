#!/bin/python

import argparse
import pysam
import gzip

DUP_SVTYPE = 'DUP'


def process_record(record):
    record = process_svtype(record)
    return record


def process_svtype(record):
    if record.info.get('SVTYPE') == DUP_SVTYPE:
        record.alts = ('<' + record.info.get('SVTYPE') + '>',)
    return record


if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(description='CleanVcf postprocessing.')
    parser.add_argument('-V', '--input', dest='input_vcf', required=True, help='Input VCF file')
    parser.add_argument('-O', '--output', dest='output_vcf', required=True, help='Output VCF file')
    args = parser.parse_args()

    # Read input files
    if args.input_vcf.endswith('.gz'):
        vcf_in = pysam.VariantFile(gzip.open(args.input_vcf, 'rt'))
    else:
        vcf_in = pysam.VariantFile(args.input_vcf)

    # Open output file
    if args.output_vcf.endswith('.gz'):
        vcf_out = pysam.VariantFile(args.output_vcf, 'wz', header=vcf_in.header)
    else:
        vcf_out = pysam.VariantFile(args.output_vcf, 'w', header=vcf_in.header.copy())

    # Process records
    for record in vcf_in:
        record = process_record(record)
        vcf_out.write(record)

    # Close files
    vcf_in.close()
    vcf_out.close()
