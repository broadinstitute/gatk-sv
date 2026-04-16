#!/usr/bin/env python3
"""Split TRV variants in a VCF by MOTIFS length.

Keeps only records with "TRV" in ID (3rd column), then writes them to:
1) single-base MOTIFS: any motif length == 1
2) two-base MOTIFS: any motif length == 2 (and no single-base motif)
3) rest TRV: all other TRV records

Usage:
  python3 split_trv_by_motifs.py --input in.vcf.gz --out-prefix out_prefix
"""

import argparse
import gzip
import os


def open_text(path, mode):
    if path.endswith('.gz'):
        return gzip.open(path, mode + 't')
    return open(path, mode)


def parse_info(info_str):
    info = {}
    for item in info_str.split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            info[k] = v
    return info


def classify_trv_by_motifs(info_dict):
    motifs_raw = info_dict.get('MOTIFS', '')
    motifs = [m for m in motifs_raw.split(',') if m]

    # Precedence: single-base first, then 2-base, then rest.
    if any(len(m) == 1 for m in motifs):
        return 'trv_motif_1bp'
    if any(len(m) == 2 for m in motifs):
        return 'trv_motif_2bp'
    return 'trv_motif_rest'


def split_vcf(input_vcf, out_prefix):
    outputs = {
        'trv_motif_1bp': f'{out_prefix}.trv.motif_1bp.vcf.gz',
        'trv_motif_2bp': f'{out_prefix}.trv.motif_2bp.vcf.gz',
        'trv_motif_rest': f'{out_prefix}.trv.motif_rest.vcf.gz',
    }

    out_files = {k: gzip.open(v, 'wt') for k, v in outputs.items()}
    stats = {
        'total_records': 0,
        'trv_records': 0,
        'trv_motif_1bp': 0,
        'trv_motif_2bp': 0,
        'trv_motif_rest': 0,
    }

    with open_text(input_vcf, 'r') as fin:
        for line in fin:
            if line.startswith('#'):
                for fh in out_files.values():
                    fh.write(line)
                continue

            stats['total_records'] += 1
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8:
                raise ValueError('Malformed VCF line with fewer than 8 columns')

            variant_id = fields[2]
            if 'TRV' not in variant_id:
                continue

            stats['trv_records'] += 1
            info_dict = parse_info(fields[7])
            category = classify_trv_by_motifs(info_dict)
            out_files[category].write(line)
            stats[category] += 1

    for fh in out_files.values():
        fh.close()

    return outputs, stats


def main():
    parser = argparse.ArgumentParser(description='Split TRV variants by MOTIFS length')
    parser.add_argument('--input', required=True, help='Input VCF (.vcf or .vcf.gz)')
    parser.add_argument('--out-prefix', required=True, help='Output prefix for split VCFs')
    args = parser.parse_args()

    if not os.path.exists(args.input):
        raise FileNotFoundError(f'Input VCF not found: {args.input}')

    outputs, stats = split_vcf(args.input, args.out_prefix)

    print('Split complete')
    print(f"Input: {args.input}")
    print(f"Total variants read: {stats['total_records']}")
    print(f"TRV variants retained: {stats['trv_records']}")
    print(f"TRV with 1bp MOTIFS: {stats['trv_motif_1bp']} -> {outputs['trv_motif_1bp']}")
    print(f"TRV with 2bp MOTIFS: {stats['trv_motif_2bp']} -> {outputs['trv_motif_2bp']}")
    print(f"TRV other MOTIFS: {stats['trv_motif_rest']} -> {outputs['trv_motif_rest']}")


if __name__ == '__main__':
    main()
