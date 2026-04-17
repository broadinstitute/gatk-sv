#!/usr/bin/env python3
"""Split a VCF into 6 subsets: TRV and 5 non-TRV categories by type/size.

Non-TRV categories are derived from REF/ALT sequence comparison:
- SNV: len(REF) == len(ALT) == 1
- INS: len(ALT) > len(REF), size = len(ALT) - len(REF)
- DEL: len(REF) > len(ALT), size = len(REF) - len(ALT)

Records with symbolic/breakend ALT alleles (e.g. <DEL>, N]chr:pos]) are left
as uncategorized.
"""

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


def _classify_ref_alt_one(ref, alt):
    """Classify a single REF/ALT allele pair.

    Returns tuple(type, size) where:
      - type in {'SNV', 'INS', 'DEL'}
      - size is None for SNV, integer for INS/DEL
    Or returns (None, None) if allele cannot be categorized.
    """
    if not ref or not alt:
        return None, None

    # Skip symbolic and breakend-like ALT representations.
    if alt.startswith('<') or '[' in alt or ']' in alt or alt == '*':
        return None, None

    rlen = len(ref)
    alen = len(alt)

    if rlen == 1 and alen == 1:
        return 'SNV', None
    if alen > rlen:
        return 'INS', alen - rlen
    if rlen > alen:
        return 'DEL', rlen - alen

    # MNV or complex substitution with no net length change.
    return None, None


def classify_variant(variant_id, ref, alt_field):
    """Classify variant into one of 6 categories.
    
    Returns: category (trv, non_trv_snv, non_trv_ins_lt50, 
                       non_trv_del_lt50, non_trv_ins_ge50, non_trv_del_ge50)
             or None if uncategorizable
    """
    if 'TRV' in variant_id:
        return 'trv'
    
    # For non-TRV, derive type/size from REF/ALT sequence comparison.
    alts = alt_field.split(',')
    per_alt = [_classify_ref_alt_one(ref, a) for a in alts]

    # If any ALT is uncategorizable, keep conservative behavior.
    if any(t is None for t, _ in per_alt):
        return None

    # If multi-allelic record mixes categories, treat as uncategorized.
    cats = {t for t, _ in per_alt}
    if len(cats) != 1:
        return None

    allele_type = next(iter(cats))

    if allele_type == 'SNV':
        return 'non_trv_snv'

    lengths = {size for t, size in per_alt if size is not None}
    if len(lengths) != 1:
        return None
    allele_length = next(iter(lengths))
    
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
            ref = fields[3]
            alt = fields[4]
            
            # Classify variant
            category = classify_variant(variant_id, ref, alt)
            
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
