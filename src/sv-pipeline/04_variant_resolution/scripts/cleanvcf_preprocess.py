#!/bin/python

import argparse
import pysam
import gzip


ME_ALT = ':ME'
DUP_SVTYPE = 'DUP'
VAR_GQ = 'varGQ'
MULTIALLELIC = 'MULTIALLELIC'
UNRESOLVED = 'UNRESOLVED'
HIGH_SR_BACKGROUND = 'HIGH_SR_BACKGROUND'
BOTHSIDES_SUPPORT = 'BOTHSIDES_SUPPORT'
REVISED_EVENT = 'REVISED_EVENT'
EV_VALUES = ['SR', 'PE', 'SR,PE', 'RD', 'BAF', 'RD,BAF']
MIN_ALLOSOME_EVENT_SIZE = 50


def read_last_column(file_path):
    result_set = set()
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip():
                columns = line.strip().split()
                result_set.add(columns[-1])
    return result_set


def process_record(record, chrX, chrY, fail_set, pass_set):
    record = process_varGQ(record)
    record = process_multiallelic(record)
    record = process_unresolved(record)
    record = process_svtype(record)
    record = process_noisy(record, fail_set)
    record = process_bothsides_support(record, pass_set)
    record = process_allosomes(record, chrX, chrY)
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


def process_svtype(record):
    if not any(ME_ALT in alt for alt in record.alts) and not record.info.get('SVTYPE') == DUP_SVTYPE:
        record.alts = ('<' + record.info.get('SVTYPE') + '>',)
    return record


def process_noisy(record, fail_set):
    if record.id in fail_set:
        record.info[HIGH_SR_BACKGROUND] = True
    return record


def process_bothsides_support(record, pass_set):
    if record.id in pass_set:
        record.info[BOTHSIDES_SUPPORT] = True
    return record


def process_allosomes(record, chrX, chrY):
    chromosome = record.chrom
    if chromosome not in (chrX, chrY):
        return record

    updated_samples = []
    sv_type = record.info.get('SVTYPE', '')
    sv_len = record.info.get('SVLEN', 0)

    if sv_type in ('DEL', 'DUP') and sv_len >= MIN_ALLOSOME_EVENT_SIZE:
        is_y = (chromosome == chrY)

        for sample in record.samples:
            genotype = record.samples[sample]
            sex = genotype.get('EXPECTED_COPY_NUMBER_FORMAT', None)

            if sex == 1:  # Male
                if is_revisable_event(record, is_y, sex):
                    record.info[REVISED_EVENT] = True
                    adjust_male_genotype(genotype, sv_type)
            elif sex == 2 and is_y:  # Female - NO_CALL if chrY
                genotype['GT'] = (None, None)
            elif sex == 0:  # Unknown sex - NO_CALL
                genotype['GT'] = (None, None)

            updated_samples.append(sample)

    return record


def is_revisable_event(record, is_y, sex):
    genotypes = record.samples.values()
    male_counts = [0, 0, 0, 0]
    female_counts = [0, 0, 0, 0]

    for genotype in genotypes:
        rd_cn = genotype.get('RD_CN', -1)
        rd_cn_val = min(rd_cn, 3) if rd_cn != -1 else -1
        if rd_cn_val == -1:
            continue

        if sex == 1:  # Male
            male_counts[rd_cn_val] += 1
        elif sex == 2:  # Female
            female_counts[rd_cn_val] += 1

    male_median = calc_median_distribution(male_counts)
    female_median = calc_median_distribution(female_counts)

    return male_median == 1 and (female_median == 0 if is_y else female_median == 2)


def adjust_male_genotype(genotype, sv_type):
    rd_cn = genotype.get('RD_CN', 0)
    genotype['RD_CN'] = rd_cn + 1
    ref_allele, alt_allele = genotype['alleles']

    if sv_type == 'DEL':
        if rd_cn >= 1:
            genotype['GT'] = (ref_allele, ref_allele)
        elif rd_cn == 0:
            genotype['GT'] = (ref_allele, alt_allele)
    elif sv_type == 'DUP':
        if rd_cn <= 1:
            genotype['GT'] = (ref_allele, ref_allele)
        elif rd_cn == 2:
            genotype['GT'] = (ref_allele, alt_allele)
        else:
            genotype['GT'] = (alt_allele, alt_allele)


def calc_median_distribution(counts):
    total = sum(counts)
    if total == 0:
        return -1

    target = total // 2
    running_total = 0
    for i, count in enumerate(counts):
        running_total += count
        if running_total >= target:
            return i * 2 if running_total > target else i * 2 + 1


if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(description='CleanVcf preprocessing.')
    parser.add_argument('-V', '--input', dest='input_vcf', required=True, help='Input VCF file')
    parser.add_argument('-O', '--output', dest='output_vcf', required=True, help='Output VCF file')
    parser.add_argument('--chrX', required=True, help='Chromosome X representation in VCF')
    parser.add_argument('--chrY', required=True, help='Chromosome Y representation in VCF')
    parser.add_argument('--fail-list', required=True, help='File with variants failing the background test')
    parser.add_argument('--pass-list', required=True, help='File with variants passing both sides')
    args = parser.parse_args()

    # Read input files
    chrX = args.chrX
    chrY = args.chrY
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
        record = process_record(record, chrX, chrY, fail_set, pass_set)
        vcf_out.write(record)

    # Close files
    vcf_in.close()
    vcf_out.close()
