#!/usr/bin/env python3
"""Split a VCF into 6 subsets: TRV and 5 non-TRV categories by type/size."""

import argparse
import gzip
import os


def open_text(path, mode):
    """Open file, handling gzip compression if needed."""
    if path.endswith('.gz'):
        return gzip.open(path, mode + 't')
    return open(path, mode)


def parse_info(info_str):
    """Parse INFO field into dictionary."""
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            info_dict[k] = v
    return info_dict


def classify_variant(variant_id, info_dict):
    """Classify variant into one of 6 categories.
    
    Returns: category (trv, non_trv_snv, non_trv_ins_lt50, 
                       non_trv_del_lt50, non_trv_ins_ge50, non_trv_del_ge50)
             or None if uncategorizable
    """
    if 'TRV' in variant_id:
        return 'trv'
    
    # For non-TRV, get allele_type and allele_length from INFO (case-insensitive)
    allele_type = info_dict.get('allele_type', '').upper()
    allele_length_str = info_dict.get('allele_length', '')
    
    if allele_type == 'SNV':
        return 'non_trv_snv'
    
    # Try to parse length. DEL records in this VCF use negative allele_length,
    # so use absolute size thresholds for both INS and DEL.
    try:
        allele_length = abs(int(allele_length_str))
    except (ValueError, TypeError):
        return None
    
    if allele_type == 'INS':
        return 'non_trv_ins_lt50' if allele_length < 50 else 'non_trv_ins_ge50'
    elif allele_type == 'DEL':
        return 'non_trv_del_lt50' if allele_length < 50 else 'non_trv_del_ge50'
    
    return None


def output_paths(prefix):
    """Generate output file paths for all 6 categories."""
    return {
        'trv': prefix + '.trv.vcf.gz',
        'non_trv_snv': prefix + '.non_trv.snv.vcf.gz',
        'non_trv_ins_lt50': prefix + '.non_trv.ins_lt50.vcf.gz',
        'non_trv_del_lt50': prefix + '.non_trv.del_lt50.vcf.gz',
        'non_trv_ins_ge50': prefix + '.non_trv.ins_ge50.vcf.gz',
        'non_trv_del_ge50': prefix + '.non_trv.del_ge50.vcf.gz',
    }


def split_vcf(input_vcf, out_prefix):
    """Split VCF into 6 subsets by category."""
    output_map = output_paths(out_prefix)
    
    # Initialize counters
    stats = {cat: 0 for cat in output_map.keys()}
    stats['uncategorized'] = 0
    
    # Open all output files
    out_files = {cat: gzip.open(path, 'wt') for cat, path in output_map.items()}
    
    with open_text(input_vcf, 'r') as fin:
        for line in fin:
            if line.startswith('#'):
                # Write header to all outputs
                for f in out_files.values():
                    f.write(line)
                continue
            
            # Parse variant
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8:
                raise ValueError('Encountered malformed VCF line with fewer than 8 columns')
            
            variant_id = fields[2]
            info_str = fields[7]
            info_dict = parse_info(info_str)
            
            # Classify variant
            category = classify_variant(variant_id, info_dict)
            
            if category is None:
                stats['uncategorized'] += 1
            else:
                out_files[category].write(line)
                stats[category] += 1
    
    # Close all output files
    for f in out_files.values():
        f.close()
    
    return {
        'input': input_vcf,
        'outputs': output_map,
        'stats': stats,
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description='Split a VCF into 6 subsets: TRV, non-TRV SNVs, and 4 size/type categories for non-TRV SVs.'
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

    result = split_vcf(args.input, args.out_prefix)
    stats = result['stats']
    outputs = result['outputs']

    print('Split complete')
    print(f"Input: {result['input']}")
    print(f"\nOutput files:")
    print(f"  TRV: {outputs['trv']} ({stats['trv']} variants)")
    print(f"  Non-TRV SNVs: {outputs['non_trv_snv']} ({stats['non_trv_snv']} variants)")
    print(f"  Non-TRV INS <50bp: {outputs['non_trv_ins_lt50']} ({stats['non_trv_ins_lt50']} variants)")
    print(f"  Non-TRV DEL <50bp: {outputs['non_trv_del_lt50']} ({stats['non_trv_del_lt50']} variants)")
    print(f"  Non-TRV INS >=50bp: {outputs['non_trv_ins_ge50']} ({stats['non_trv_ins_ge50']} variants)")
    print(f"  Non-TRV DEL >=50bp: {outputs['non_trv_del_ge50']} ({stats['non_trv_del_ge50']} variants)")
    print(f"  Uncategorized: {stats['uncategorized']} variants")
    
    total = sum(v for k, v in stats.items() if k != 'uncategorized')
    print(f"\nTotal variants: {total}")


if __name__ == '__main__':
    main()
