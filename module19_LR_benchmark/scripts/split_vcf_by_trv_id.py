#!/usr/bin/env python3

import argparse
import gzip
import os


def open_text(path, mode):
    if path.endswith('.gz'):
        return gzip.open(path, mode + 't')
    return open(path, mode)


def output_paths(prefix):
    trv_path = prefix + '.trv_id.vcf.gz'
    non_trv_path = prefix + '.non_trv_id.vcf.gz'
    return trv_path, non_trv_path


def split_vcf(input_vcf, out_prefix):
    trv_out, non_trv_out = output_paths(out_prefix)

    total_records = 0
    trv_records = 0
    non_trv_records = 0

    with open_text(input_vcf, 'r') as fin, gzip.open(trv_out, 'wt') as f_trv, gzip.open(non_trv_out, 'wt') as f_non:
        for line in fin:
            if line.startswith('#'):
                f_trv.write(line)
                f_non.write(line)
                continue

            total_records += 1
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 3:
                raise ValueError('Encountered malformed VCF line with fewer than 3 columns')

            variant_id = fields[2]
            if 'TRV' in variant_id:
                f_trv.write(line)
                trv_records += 1
            else:
                f_non.write(line)
                non_trv_records += 1

    return {
        'input': input_vcf,
        'trv_output': trv_out,
        'non_trv_output': non_trv_out,
        'total_records': total_records,
        'trv_records': trv_records,
        'non_trv_records': non_trv_records,
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description='Split a VCF into two subsets by whether ID (3rd column) contains TRV.'
    )
    parser.add_argument(
        '--input',
        required=True,
        help='Input VCF path (.vcf or .vcf.gz)',
    )
    parser.add_argument(
        '--out-prefix',
        required=True,
        help='Output prefix path for generated files',
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.exists(args.input):
        raise FileNotFoundError(f'Input VCF not found: {args.input}')

    stats = split_vcf(args.input, args.out_prefix)

    print('Split complete')
    print(f"Input: {stats['input']}")
    print(f"TRV subset: {stats['trv_output']} ({stats['trv_records']} variants)")
    print(f"Non-TRV subset: {stats['non_trv_output']} ({stats['non_trv_records']} variants)")
    print(f"Total variants processed: {stats['total_records']}")


if __name__ == '__main__':
    main()
