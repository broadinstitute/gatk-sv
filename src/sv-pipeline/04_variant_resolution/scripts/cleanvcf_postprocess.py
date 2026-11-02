#!/bin/python

import argparse
import gzip

import pysam
from svtk.utils import get_called_samples


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


def process_revised_allosomes(record, unknown_sex_samples):
    if record.info['REVISED_EVENT']:
        for sample in record.samples:
            genotype = record.samples[sample]
            ecn = genotype.get('ECN')
            if sample not in unknown_sex_samples and ecn == 1:
                cn = record.samples[sample].get('RD_CN')
                if cn is not None and int(cn) > 0:
                    cn = int(cn)
                    record.samples[sample]['RD_CN'] = cn - 1
                    if 'CN' in record.samples[sample]:
                        record.samples[sample]['CN'] = cn - 1
    return record


def process_svtype(record):
    if record.info.get('SVTYPE') == 'DUP':
        record.alts = ('<' + record.info.get('SVTYPE') + '>',)
    return record


def process_uncalled_genotypes(record):
    if len(get_called_samples(record)) == 0:
        return None
    return record


def process_record(record, unknown_sex_samples):
    record = process_revised_allosomes(record, unknown_sex_samples)
    record = process_svtype(record)
    record = process_uncalled_genotypes(record)
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='CleanVcf postprocessing.')
    parser.add_argument('-V', '--input', dest='input_vcf', required=True, help='Input VCF file')
    parser.add_argument('-O', '--output', dest='output_vcf', required=True, help='Output VCF file')
    parser.add_argument('--ped-file', required=True, help='PED file containing sample sex annotations')
    args = parser.parse_args()

    if args.input_vcf.endswith('.gz'):
        vcf_in = pysam.VariantFile(gzip.open(args.input_vcf, 'rt'))
    else:
        vcf_in = pysam.VariantFile(args.input_vcf)

    if args.output_vcf.endswith('.gz'):
        vcf_out = pysam.VariantFile(args.output_vcf, 'wz', header=vcf_in.header)
    else:
        vcf_out = pysam.VariantFile(args.output_vcf, 'w', header=vcf_in.header.copy())

    unknown_sex_samples = load_unknown_sex_samples(args.ped_file)

    for record in vcf_in:
        processed = process_record(record, unknown_sex_samples)
        if processed:
            vcf_out.write(processed)

    vcf_in.close()
    vcf_out.close()
