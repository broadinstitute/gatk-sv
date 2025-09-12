#!/Users/kjaising/miniconda3/envs/venv/bin/python

"""
Extract filter information from VCF files and assign variants to filter categories.
"""

import pandas as pd
import argparse
import pysam
import os

def get_filter_category(filter_val, evidence_val):
    """
    Assign a variant to a filter category based on FILTER and EVIDENCE values.
    
    Categories:
    - 'PASS+MULTIALLELIC': PASS or MULTIALLELIC filter
    - 'Wham-only Non-SR': HIGH_ALGORITHM_FDR filter and EVIDENCE contains non-SR evidence (PE, RD, etc.)
    - 'Wham-only SR': HIGH_ALGORITHM_FDR filter and EVIDENCE contains only "SR"
    - 'Other': Everything else
    """
    # Handle case where filter_val might be a list or string
    if isinstance(filter_val, list):
        filters = filter_val
    else:
        filters = [filter_val] if filter_val else []
    
    # Convert to set for easier checking
    filter_set = set(filters)
    
    # PASS+MULTIALLELIC: variants with PASS or MULTIALLELIC filter
    if 'PASS' in filter_set or 'MULTIALLELIC' in filter_set:
        return 'PASS+MULTIALLELIC'
    
    # Wham-only categories: variants with HIGH_ALGORITHM_FDR filter
    if 'HIGH_ALGORITHM_FDR' in filter_set:
        if evidence_val:
            # Handle different evidence value formats (tuple from pysam vs string)
            if isinstance(evidence_val, (tuple, list)):
                evidence_types = set(evidence_val)
            else:
                # Parse evidence types (comma-separated string)
                evidence_types = set([ev.strip() for ev in str(evidence_val).split(',')])
            
            # If evidence contains only SR, it's Wham-only SR
            if evidence_types == {'SR'}:
                return 'Wham-only SR'
            else:
                # If evidence contains other types (PE, RD, etc.), it's Wham-only Non-SR
                return 'Wham-only Non-SR'
        else:
            # If no evidence info available, default to Wham-only SR
            return 'Wham-only SR'
    
    # Everything else
    return 'Other'

def extract_filter_info(vcf_file, output_file):
    """
    Extract filter information from a VCF file and create a mapping file.
    """
    print(f"Extracting filter info from {os.path.basename(vcf_file)}")
    
    filter_data = []
    
    # Open VCF file
    vcf = pysam.VariantFile(vcf_file)
    
    for record in vcf:
        # Get variant ID (name field)
        variant_id = record.id if record.id else f"{record.chrom}_{record.start}_{record.stop}"
        
        # Get FILTER values
        filter_vals = list(record.filter) if record.filter else ['PASS']
        
        # Get EVIDENCE value from INFO field (may not exist in LR VCFs)
        try:
            evidence_val = record.info.get('EVIDENCE', None)
        except (ValueError, KeyError):
            # If EVIDENCE field doesn't exist or is invalid, set to None
            evidence_val = None
        
        # Assign filter category
        filter_category = get_filter_category(filter_vals, evidence_val)
        
        filter_data.append({
            'variant_id': variant_id,
            'wrapped_variant_id': f'__{variant_id}__',  # Add wrapper for matching
            'filter_category': filter_category,
            'filters': ','.join(filter_vals),
            'evidence': evidence_val
        })
    
    vcf.close()
    
    # Create DataFrame and save
    df = pd.DataFrame(filter_data)
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"  Extracted filter info for {len(df)} variants")
    print(f"  Filter category breakdown:")
    for category, count in df['filter_category'].value_counts().items():
        print(f"    {category}: {count}")
    
    return output_file

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--output', required=True, help='Output filter info file')
    
    args = parser.parse_args()
    
    extract_filter_info(args.vcf, args.output)

if __name__ == '__main__':
    main() 