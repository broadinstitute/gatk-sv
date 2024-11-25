#!/bin/python

import argparse
import pysam

# Constants
EV = 'EV'
VAR_GQ = 'VAR_GQ'
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

def add_header_lines(header):
    header.add_line('##FILTER=<ID=UNRESOLVED,Description="Variant is unresolved">')
    header.add_line('##INFO=<ID=HIGH_SR_BACKGROUND,Number=0,Type=Flag,Description="High number of SR splits in background samples indicating messy region">')
    header.add_line('##INFO=<ID=BOTHSIDES_SUPPORT,Number=0,Type=Flag,Description="Variant has read-level support for both sides of breakpoint">')
    header.add_line('##INFO=<ID=REVISED_EVENT,Number=0,Type=Flag,Description="Variant has been revised due to a copy number mismatch">')

def process_record(record, fail_set, pass_set):
    record = process_EV(record)
    record = process_varGQ(record)
    record = process_multiallelic(record)
    record = process_unresolved(record)
    record = process_noisy(record, fail_set)
    record = process_bothsides_support(record, pass_set)
    return record

def process_EV(record):
    for sample in record.samples:
        genotype = record.samples[sample]
        if EV in genotype and genotype[EV] is not None:
            ev_attribute = genotype[EV]
            try:
                ev_index = int(ev_attribute)
                if 0 <= ev_index < len(EV_VALUES):
                    genotype[EV] = EV_VALUES[ev_index]
            except ValueError:
                pass
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
    parser = argparse.ArgumentParser(description='Process VCF variants.')
    parser.add_argument('--chr-X', dest='chrX', default='chrX', help='chrX column name')
    parser.add_argument('--chr-Y', dest='chrY', default='chrY', help='chrY column name')
    parser.add_argument('--fail-list', required=True, help='File with variants failing the background test')
    parser.add_argument('--pass-list', required=True, help='File with variants passing both sides')
    parser.add_argument('-O', '--output', dest='output_vcf', required=True, help='Output VCF name')
    parser.add_argument('input_vcf', help='Input VCF file')
    args = parser.parse_args()

    # Read noisy and bothsides support events into sets
    fail_set = read_last_column(args.fail_list)
    pass_set = read_last_column(args.pass_list)

    # Open input VCF
    vcf_in = pysam.VariantFile(args.input_vcf)

    # Modify header
    header = vcf_in.header.copy()
    add_header_lines(header)

    # Open output VCF
    vcf_out = pysam.VariantFile(args.output_vcf, 'w', header=header)

    # Process and write variants
    for record in vcf_in:
        record = process_record(record, fail_set, pass_set)
        vcf_out.write(record)

    # Close files
    vcf_in.close()
    vcf_out.close()
