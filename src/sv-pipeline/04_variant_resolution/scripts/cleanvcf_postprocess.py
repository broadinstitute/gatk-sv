#!/bin/python

import argparse
import gzip

import pysam
from svtk.utils import get_called_samples

DUP_SVTYPE = 'DUP'


def process_svtype(record):
    if record.info.get('SVTYPE') == DUP_SVTYPE:
        record.alts = ('<' + record.info.get('SVTYPE') + '>',)
    return record


def process_uncalled_genotypes(record):
    if len(get_called_samples(record)) == 0:
        return None
    return record


def process_record(record):
    record = process_svtype(record)
    record = process_uncalled_genotypes(record)
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='CleanVcf postprocessing.')
    parser.add_argument('-V', '--input', dest='input_vcf', required=True, help='Input VCF file')
    parser.add_argument('-O', '--output', dest='output_vcf', required=True, help='Output VCF file')
    args = parser.parse_args()

    if args.input_vcf.endswith('.gz'):
        vcf_in = pysam.VariantFile(gzip.open(args.input_vcf, 'rt'))
    else:
        vcf_in = pysam.VariantFile(args.input_vcf)

    if args.output_vcf.endswith('.gz'):
        vcf_out = pysam.VariantFile(args.output_vcf, 'wz', header=vcf_in.header)
    else:
        vcf_out = pysam.VariantFile(args.output_vcf, 'w', header=vcf_in.header.copy())

    for record in vcf_in:
        processed = process_record(record)
        if processed:
            vcf_out.write(processed)

    vcf_in.close()
    vcf_out.close()
