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
MIN_ALLOSOME_EVENT_SIZE = 5000


def load_unknown_sex_samples(ped_path):
    unknown_samples = set()
    with open(ped_path, 'r') as ped:
        for line in ped:
            columns = line.strip().split()
            if not columns:
                continue
            sample = columns[1]
            sex = columns[4]
            if sex not in ('1', '2'):
                unknown_samples.add(sample)
    return unknown_samples


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


def process_allosomes(record, chrX, chrY, unknown_sex_samples):
    chromosome = record.chrom
    if chromosome not in (chrX, chrY):
        return record

    is_y = (chromosome == chrY)
    sv_type = record.info.get('SVTYPE', '')
    sv_len = record.info.get('SVLEN', 0)

    if sv_type in ('DEL', 'DUP') and sv_len >= MIN_ALLOSOME_EVENT_SIZE:
        if is_revisable_event(record, is_y, unknown_sex_samples):
            record.info[REVISED_EVENT] = True

    for sample in record.samples:
        genotype = record.samples[sample]
        ecn = genotype.get('ECN')
        if sample in unknown_sex_samples:
            record.samples[sample]['GT'] = (None, None)
        elif ecn == 0:
            if is_y:
                clear_genotype_fields(genotype)
        elif ecn == 1:
            if record.info[REVISED_EVENT]:
                adjust_male_genotype(genotype, sv_type)

    return record


def clear_genotype_fields(genotype):
    for key in genotype.keys():
        genotype[key] = None
    genotype['GT'] = (None, None)


def is_revisable_event(record, is_y, unknown_sex_samples):
    male_counts = [0, 0, 0, 0]
    female_counts = [0, 0, 0, 0]

    for sample in record.samples:
        sample_ecn = record.samples[sample]['ECN']
        rd_cn = record.samples[sample]['RD_CN']
        if rd_cn is None:
            continue

        rd_cn_val = min(int(rd_cn), 3)
        if sample in unknown_sex_samples:
            continue
        elif sample_ecn == 1:
            male_counts[rd_cn_val] += 1
        elif is_y and sample_ecn == 0:
            female_counts[rd_cn_val] += 1
        elif not is_y and sample_ecn == 2:
            female_counts[rd_cn_val] += 1

    male_median = calc_median_distribution(male_counts)
    female_median = calc_median_distribution(female_counts)
    return male_median == 1.0 and (female_median == 0.0 if is_y else female_median == 2.0)


def adjust_male_genotype(genotype, sv_type):
    rd_cn = genotype.get('RD_CN')
    if rd_cn is None:
        return

    genotype['RD_CN'] = rd_cn + 1
    if sv_type == 'DEL':
        if rd_cn >= 1:
            genotype['GT'] = (0, 0)
        elif rd_cn == 0:
            genotype['GT'] = (0, 1)
    elif sv_type == 'DUP':
        if rd_cn <= 1:
            genotype['GT'] = (0, 0)
        elif rd_cn == 2:
            genotype['GT'] = (0, 1)
        else:
            genotype['GT'] = (1, 1)


def calc_median_distribution(counts):
    total = sum(counts)
    if total == 0:
        return -1

    target = total / 2.0
    running_total = 0
    for i, count in enumerate(counts):
        running_total += count
        if running_total == target:
            return i + 0.5
        elif running_total > target:
            return float(i)
    return -1


def read_last_column(file_path):
    result_set = set()
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip():
                columns = line.strip().split()
                result_set.add(columns[-1])
    return result_set


def process_record(record, chrX, chrY, fail_set, pass_set, unknown_sex_samples):
    record = process_varGQ(record)
    record = process_multiallelic(record)
    record = process_unresolved(record)
    record = process_svtype(record)
    record = process_noisy(record, fail_set)
    record = process_bothsides_support(record, pass_set)
    record = process_allosomes(record, chrX, chrY, unknown_sex_samples)
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='CleanVcf preprocessing.')
    parser.add_argument('-V', '--input', dest='input_vcf', required=True, help='Input VCF file')
    parser.add_argument('-O', '--output', dest='output_vcf', required=True, help='Output VCF file')
    parser.add_argument('--chrX', required=True, help='Chromosome X representation in VCF')
    parser.add_argument('--chrY', required=True, help='Chromosome Y representation in VCF')
    parser.add_argument('--fail-list', required=True, help='File with variants failing the background test')
    parser.add_argument('--pass-list', required=True, help='File with variants passing both sides')
    parser.add_argument('--ped-file', required=True, help='PED file containing sample sex annotations')
    args = parser.parse_args()

    chrX = args.chrX
    chrY = args.chrY
    fail_set = read_last_column(args.fail_list)
    pass_set = read_last_column(args.pass_list)
    unknown_sex_samples = load_unknown_sex_samples(args.ped_file)

    if args.input_vcf.endswith('.gz'):
        vcf_in = pysam.VariantFile(gzip.open(args.input_vcf, 'rt'))
    else:
        vcf_in = pysam.VariantFile(args.input_vcf)

    if args.output_vcf.endswith('.gz'):
        vcf_out = pysam.VariantFile(args.output_vcf, 'wz', header=vcf_in.header)
    else:
        vcf_out = pysam.VariantFile(args.output_vcf, 'w', header=vcf_in.header.copy())

    for record in vcf_in:
        record = process_record(record, chrX, chrY, fail_set, pass_set, unknown_sex_samples)
        vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()
