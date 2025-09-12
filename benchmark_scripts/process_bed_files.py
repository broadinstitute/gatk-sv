#!/Users/kjaising/miniconda3/envs/venv/bin/python3
"""
Process BED files to create query and ref formatted outputs.
Handles both sr.bed (with AF column) and lr.bed (without AF column) formats.
"""

import pandas as pd
import sys
import os
import numpy as np
import argparse

def chromosome_sort_key(chrom):
    """
    Create a sort key for chromosome names to ensure natural ordering.
    Returns a tuple for sorting: (type, number) where type=0 for autosomes, 1 for sex chromosomes
    """
    chrom_str = str(chrom).lower()
    
    if chrom_str.startswith('chr'):
        chrom_str = chrom_str[3:]
    
    # Handle autosomes (1-22)
    if chrom_str.isdigit():
        return (0, int(chrom_str))
    
    # Handle sex chromosomes
    if chrom_str == 'x':
        return (1, 1)
    elif chrom_str == 'y':
        return (1, 2)
    elif chrom_str == 'm' or chrom_str == 'mt':
        return (1, 3)
    
    # Handle any other chromosomes
    return (2, chrom_str)

def sort_dataframe_by_genomic_position(df):
    """
    Sort dataframe by chromosome (natural order) then by start position.
    """
    # Create sort keys for chromosomes
    df['_sort_key'] = df['#chrom'].apply(chromosome_sort_key)
    
    # Sort by chromosome sort key, then by start position
    df_sorted = df.sort_values(['_sort_key', 'start']).drop('_sort_key', axis=1)
    
    return df_sorted

def process_bed_file(input_file, sample_name, is_lr=False, limit=None):
    """
    Process a BED file according to the specified rules.
    
    Args:
        input_file (str): Path to the input BED file
        sample_name (str): Sample name to use in output files
        is_lr (bool): Whether this is a long-read (LR) file
        limit (int): Limit the number of variants processed
    """
    file_type = "LR" if is_lr else "SR"
    print(f"Processing {file_type} BED: {os.path.basename(input_file)}")
    
    # Read the BED file
    df = pd.read_csv(input_file, sep='\t')
    
    # Clean data: remove rows with missing or non-finite SVLEN values
    initial_rows = len(df)
    df = df.dropna(subset=['SVLEN'])
    df = df[np.isfinite(df['SVLEN'])]
    if initial_rows != len(df):
        print(f"  Removed {initial_rows - len(df)} rows with invalid SVLEN")
    
    # Make a copy for processing
    processed_df = df.copy()
    
    # Convert SVLEN to numeric if it's not already
    processed_df['SVLEN'] = pd.to_numeric(processed_df['SVLEN'], errors='coerce')
    
    # Rule 1: For INS variants, set end = end + SVLEN
    ins_mask = processed_df['SVTYPE'] == 'INS'
    processed_df.loc[ins_mask, 'end'] = processed_df.loc[ins_mask, 'end'] + processed_df.loc[ins_mask, 'SVLEN']
    if ins_mask.sum() > 0:
        print(f"  Updated {ins_mask.sum()} INS variants")
    
    # Rule 2: For DEL variants, set SVLEN = abs(SVLEN)
    del_mask = processed_df['SVTYPE'] == 'DEL'
    processed_df.loc[del_mask, 'SVLEN'] = processed_df.loc[del_mask, 'SVLEN'].abs()
    if del_mask.sum() > 0:
        print(f"  Updated {del_mask.sum()} DEL variants")
    
    # Create query file
    base_name = os.path.splitext(input_file)[0]
    query_file = f"{base_name}.query.bed"
    
    query_df = pd.DataFrame({
        '#chrom': processed_df['#chrom'],
        'start': processed_df['start'],
        'end': processed_df['end'],
        'name': '__' + processed_df['name'].astype(str) + '__',
        'SVTYPE': processed_df['SVTYPE'],
        'SVLEN': processed_df['SVLEN'].astype(int)
    })
    
    # Sort query dataframe by genomic position
    query_df = sort_dataframe_by_genomic_position(query_df)
    
    if limit:
        query_df = query_df.head(limit)
    
    query_df.to_csv(query_file, sep='\t', index=False)
    
    # Create ref file
    ref_file = f"{base_name}.ref.bed"
    
    # For LR files, filter ref output to only include variants with length > 50
    ref_processed_df = processed_df.copy()
    if is_lr:
        initial_ref_rows = len(ref_processed_df)
        ref_processed_df = ref_processed_df[ref_processed_df['SVLEN'].abs() > 50]
        if initial_ref_rows != len(ref_processed_df):
            print(f"  LR ref: filtered to {len(ref_processed_df)} variants (length > 50)")
    
    # Check if AF column exists
    has_af = 'AF' in ref_processed_df.columns
    if has_af:
        ref_processed_df['AF'] = ref_processed_df['AF'].replace({None: 0, np.nan: 0, '': 0})
        ref_processed_df['AF'] = pd.to_numeric(ref_processed_df['AF'], errors='coerce').fillna(0)
    else:
        ref_processed_df['AF'] = 0
    
    ref_df = pd.DataFrame({
        '#chrom': ref_processed_df['#chrom'],
        'start': ref_processed_df['start'],
        'end': ref_processed_df['end'],
        'VID': '__' + ref_processed_df['name'].astype(str) + '__',
        'svtype': ref_processed_df['SVTYPE'],
        'length': ref_processed_df['SVLEN'].astype(int),
        'AF': ref_processed_df['AF'],
        'samples': sample_name
    })
    
    # Sort ref dataframe by genomic position
    ref_df = sort_dataframe_by_genomic_position(ref_df)
    
    if limit:
        ref_df = ref_df.head(limit)
    
    ref_df.to_csv(ref_file, sep='\t', index=False)
    print(f"  Created query and ref files")
    
    return query_file, ref_file

def main():
    parser = argparse.ArgumentParser(description='Process BED files to create query and ref formatted outputs')
    parser.add_argument('--sample', required=True, help='Sample name to use in output files')
    parser.add_argument('--sr', required=True, help='Path to short-read (SR) BED file')
    parser.add_argument('--lr', required=True, help='Path to long-read (LR) BED file')
    parser.add_argument('--sr-limit', type=int, help='Limit number of SR variants (for testing)')
    parser.add_argument('--lr-limit', type=int, help='Limit number of LR variants (for testing)')
    
    args = parser.parse_args()
    
    files_to_process = [
        (args.sr, False),  # SR file, not LR
        (args.lr, True)    # LR file, is LR
    ]
    
    for file_path, is_lr in files_to_process:
        if os.path.exists(file_path):
            try:
                # Determine the limit for this file type
                limit = args.lr_limit if is_lr else args.sr_limit
                query_file, ref_file = process_bed_file(file_path, args.sample, is_lr, limit)
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
        else:
            print(f"File not found: {file_path}")

if __name__ == "__main__":
    main() 