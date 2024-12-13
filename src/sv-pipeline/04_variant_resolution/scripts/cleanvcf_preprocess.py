#!/bin/python

import argparse
import pysam
import gzip


VAR_GQ = 'varGQ'
MULTIALLELIC = 'MULTIALLELIC'
UNRESOLVED = 'UNRESOLVED'
HIGH_SR_BACKGROUND = 'HIGH_SR_BACKGROUND'
BOTHSIDES_SUPPORT = 'BOTHSIDES_SUPPORT'
REVISED_EVENT = 'REVISED_EVENT'
EV_VALUES = ['SR', 'PE', 'SR,PE', 'RD', 'BAF', 'RD,BAF']


def read_last_column(file_path):
    result_set = set()
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip():
                columns = line.strip().split()
                result_set.add(columns[-1])
    return result_set


def process_record(record, fail_set, pass_set):
    record = process_varGQ(record)
    record = process_multiallelic(record)
    record = process_unresolved(record)
    record = process_noisy(record, fail_set)
    record = process_bothsides_support(record, pass_set)
    return record


def process_varGQ(record):
    if VAR_GQ in record.info:
        var_gq = record.info[VAR_GQ]
        if isinstance(var_gq, list):
            var_gq = var_gq[0]
        del record.info[VAR_GQ]
        record.qual = var_gq
    return record


def process_multiallelic(record):
    if MULTIALLELIC in record.info:
        del record.info[MULTIALLELIC]
    return record


def process_unresolved(record):
    if UNRESOLVED in record.info:
        del record.info[UNRESOLVED]
        record.filter.add(UNRESOLVED)
    return record


def process_noisy(record, fail_set):
    if record.id in fail_set:
        record.info[HIGH_SR_BACKGROUND] = True
    return record


def process_bothsides_support(record, pass_set):
    if record.id in pass_set:
        record.info[BOTHSIDES_SUPPORT] = True
    return record


if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(description='CleanVcf preprocessing.')
    parser.add_argument('-V', '--input', dest='input_vcf', required=True, help='Input VCF file')
    parser.add_argument('-O', '--output', dest='output_vcf', required=True, help='Output VCF file')
    parser.add_argument('--fail-list', required=True, help='File with variants failing the background test')
    parser.add_argument('--pass-list', required=True, help='File with variants passing both sides')
    args = parser.parse_args()

    # Read input files
    fail_set = read_last_column(args.fail_list)
    pass_set = read_last_column(args.pass_list)
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
        record = process_record(record, fail_set, pass_set)
        vcf_out.write(record)
    
    # Close files
    vcf_in.close()
    vcf_out.close()
