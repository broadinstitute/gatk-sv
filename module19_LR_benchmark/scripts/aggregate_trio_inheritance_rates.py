#!/usr/bin/env python3
"""
Aggregate trio inheritance rates across multiple summary files.

This script reads all trio_inheritance.*.summary.tsv files and computes:
- Average error rates across families (CHS, PUR, YRI)
- de_novo_rate = de_novo / total_informative
- paternal_error_rate = paternal_inheritance_error / total_informative
- maternal_error_rate = maternal_inheritance_error / total_informative
- either_side_error_rate = error_from_either_side / total_informative

Generates separate rows for "original" (with missing alleles) and
"no_missing_members" (strict) types.

Output: Single condensed TSV with all subsets and their error rates.
"""

import os
import csv
import argparse
from collections import defaultdict


def aggregate_rates(input_dir, subsets, output_file):
    """Aggregate inheritance rates from all summary files."""
    
    results = []
    
    for subset in subsets:
        summary_file = os.path.join(input_dir, f"trio_inheritance.hgsvc_hprc.chr22.{subset}.summary.tsv")
        
        if not os.path.exists(summary_file):
            print(f"⚠ {subset}: Summary file not found at {summary_file}")
            continue
        
        # Group data by summary_type
        by_type = defaultdict(lambda: {
            'de_novo': [],
            'paternal_err': [],
            'maternal_err': [],
            'error_both': [],
            'total': []
        })
        
        with open(summary_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                summary_type = row.get('summary_type')
                
                de_novo = int(row.get('de_novo', 0))
                pat_err = int(row.get('paternal_inheritance_error', 0))
                mat_err = int(row.get('maternal_inheritance_error', 0))
                error_both = int(row.get('error_from_either_side', 0))
                total = int(row.get('total_informative', 0))
                
                by_type[summary_type]['de_novo'].append(de_novo)
                by_type[summary_type]['paternal_err'].append(pat_err)
                by_type[summary_type]['maternal_err'].append(mat_err)
                by_type[summary_type]['error_both'].append(error_both)
                by_type[summary_type]['total'].append(total)
        
        # Calculate rates (average across families)
        for summary_type in sorted(by_type.keys()):
            data = by_type[summary_type]
            
            # Average the raw counts first across families
            avg_de_novo = sum(data['de_novo']) / len(data['de_novo']) if data['de_novo'] else 0
            avg_pat_err = sum(data['paternal_err']) / len(data['paternal_err']) if data['paternal_err'] else 0
            avg_mat_err = sum(data['maternal_err']) / len(data['maternal_err']) if data['maternal_err'] else 0
            avg_error_both = sum(data['error_both']) / len(data['error_both']) if data['error_both'] else 0
            avg_total = sum(data['total']) / len(data['total']) if data['total'] else 0
            
            # Calculate rates
            de_novo_rate = (avg_de_novo / avg_total) if avg_total > 0 else 0
            paternal_rate = (avg_pat_err / avg_total) if avg_total > 0 else 0
            maternal_rate = (avg_mat_err / avg_total) if avg_total > 0 else 0
            either_rate = (avg_error_both / avg_total) if avg_total > 0 else 0
            
            results.append({
                'subset': subset,
                'type': summary_type,
                'de_novo_rate': de_novo_rate,
                'paternal_rate': paternal_rate,
                'maternal_rate': maternal_rate,
                'either_rate': either_rate,
                'avg_total': avg_total
            })
    
    # Write output
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(
            f,
            delimiter='\t',
            fieldnames=['subset', 'type', 'de_novo_rate', 'paternal_rate', 'maternal_rate', 'either_rate', 'avg_total']
        )
        writer.writeheader()
        for row in results:
            writer.writerow(row)
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate trio inheritance error rates across all VCF subsets'
    )
    parser.add_argument(
        '--input-dir',
        required=True,
        help='Directory containing trio_inheritance.*.summary.tsv files'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output TSV file path'
    )
    parser.add_argument(
        '--subsets',
        nargs='+',
        default=[
            'trv',
            'non_trv.snv', 'non_trv.del_lt50', 'non_trv.del_ge50', 'non_trv.ins_lt50', 'non_trv.ins_ge50',
            'non_tr.snv', 'non_tr.del_1_49', 'non_tr.del_50plus', 'non_tr.ins_1_49', 'non_tr.ins_50plus'
        ],
        help='VCF subset names to process'
    )
    
    args = parser.parse_args()
    
    print(f"Reading summary files from: {args.input_dir}")
    print(f"Processing subsets: {', '.join(args.subsets)}")
    
    results = aggregate_rates(args.input_dir, args.subsets, args.output)
    
    print(f"\n✓ Output written to: {args.output}")
    print(f"✓ Processed {len(results)} rows ({len(results)//2} subsets × 2 types)")


if __name__ == '__main__':
    main()
