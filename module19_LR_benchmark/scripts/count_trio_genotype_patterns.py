#!/usr/bin/env python3

import argparse
import csv
import gzip
from collections import Counter


def open_text(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def normalize_gt(sample_field):
    gt = sample_field.split(':', 1)[0].replace('|', '/')
    if gt in {'.', './.', '.|.'}:
        return './.'

    alleles = gt.split('/')
    if len(alleles) != 2:
        return gt

    if '.' in alleles:
        non_missing = sorted([allele for allele in alleles if allele != '.'])
        missing = ['.'] * alleles.count('.')
        return '/'.join(non_missing + missing)

    try:
        return '/'.join(sorted(alleles, key=lambda allele: int(allele)))
    except ValueError:
        return '/'.join(sorted(alleles))


def read_trios(trio_file):
    trios = []
    with open(trio_file, 'r', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'fam', 'father', 'mother', 'child'}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError('Trio file must be tab-delimited with columns: fam, father, mother, child')

        for row in reader:
            trios.append((row['fam'], row['father'], row['mother'], row['child']))

    if not trios:
        raise ValueError('Trio file is empty')

    return trios


def read_vcf_header(vcf_path):
    with open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith('#CHROM'):
                header = line.rstrip('\n').split('\t')
                return {sample: index for index, sample in enumerate(header[9:], start=9)}
    raise ValueError('VCF header line starting with #CHROM was not found')


def count_patterns(vcf_path, trios, sample_to_idx, pass_only=False):
    pattern_counts = {family_id: Counter() for family_id, _, _, _ in trios}

    with open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue

            fields = line.rstrip('\n').split('\t')
            if pass_only and fields[6] != 'PASS':
                continue

            for family_id, father, mother, child in trios:
                father_gt = normalize_gt(fields[sample_to_idx[father]])
                mother_gt = normalize_gt(fields[sample_to_idx[mother]])
                child_gt = normalize_gt(fields[sample_to_idx[child]])
                pattern_counts[family_id][(father_gt, mother_gt, child_gt)] += 1

    return pattern_counts


def write_output(output_path, pattern_counts):
    with open(output_path, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['father_gt', 'mo_gt', 'child_gt', 'count_variants', 'family_ID'])

        for family_id in sorted(pattern_counts):
            family_counts = pattern_counts[family_id]
            for (father_gt, mother_gt, child_gt), count in sorted(
                family_counts.items(), key=lambda item: (-item[1], item[0])
            ):
                writer.writerow([father_gt, mother_gt, child_gt, count, family_id])


def parse_args():
    parser = argparse.ArgumentParser(
        description='Count father/mother/child genotype combinations per family from a VCF'
    )
    parser.add_argument('--vcf', required=True, help='Input VCF or VCF.GZ file')
    parser.add_argument(
        '--trios',
        required=True,
        help='Tab-delimited trio file with columns: fam, father, mother, child',
    )
    parser.add_argument('--output', required=True, help='Output TSV path')
    parser.add_argument(
        '--pass-only',
        action='store_true',
        help='If set, count only variants with FILTER=PASS',
    )
    return parser.parse_args()


def main():
    args = parse_args()
    trios = read_trios(args.trios)
    sample_to_idx = read_vcf_header(args.vcf)

    missing_samples = []
    for family_id, father, mother, child in trios:
        for sample in (father, mother, child):
            if sample not in sample_to_idx:
                missing_samples.append((family_id, sample))

    if missing_samples:
        missing_text = ', '.join(f'{family}:{sample}' for family, sample in missing_samples)
        raise ValueError(f'Samples missing from VCF header: {missing_text}')

    pattern_counts = count_patterns(args.vcf, trios, sample_to_idx, pass_only=args.pass_only)
    write_output(args.output, pattern_counts)


if __name__ == '__main__':
    main()