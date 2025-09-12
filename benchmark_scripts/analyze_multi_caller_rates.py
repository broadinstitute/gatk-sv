#!/Users/kjaising/miniconda3/envs/venv/bin/python3

"""
Analyze match rates across multiple long-read callers and samples.
Generates coalesced tables showing how many callers each variant was matched in.
"""

import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import os
from collections import defaultdict

def clean_vid(vid):
    """Remove the prefix and suffix "__" from variant IDs"""
    if vid.startswith('__') and vid.endswith('__'):
        return vid[2:-2]
    return vid

def get_length_bucket(length):
    """Assign a length bucket."""
    if length < 50:
        return '<50'
    elif length < 100:
        return '50-100'
    elif length < 500:
        return '100-500'
    elif length < 5000:
        return '500-5000'
    elif length < 10000:
        return '5000-10000'
    else:
        return '>10000'

def load_sample_data(sample_id, input_dirs, query_type):
    """
    Load overlap results for a single sample across multiple callers.
    
    Args:
        sample_id: Sample identifier
        input_dirs: List of input directories (e.g., ['sniffles_results', 'cutesv_results'])
        query_type: Either 'lr_query' or 'sr_query'
    
    Returns:
        Dictionary mapping caller_name -> DataFrame with overlap results
    """
    sample_data = {}
    
    for input_dir in input_dirs:
        if input_dir.endswith('_results'):
            caller_name = input_dir.split('/')[-1][:-8]
        else:
            caller_name = input_dir
        
        query_file = os.path.join(input_dir, sample_id, f"{sample_id}.{query_type}.bed")
        
        if not os.path.exists(query_file):
            continue
        
        df = pd.read_csv(query_file, sep='\t', comment='#', header=None,
                        names=['chr', 'start', 'end', 'VID', 'svtype', 'length', 'AF', 
                                'ovr1a', 'ovr1b', 'ovr2a', 'ovr2b', 'ovr3'])
        df['clean_VID'] = df['VID'].apply(clean_vid)
        df['caller'] = caller_name
        sample_data[caller_name] = df
    
    return sample_data

def load_annotation_data(sample_id, input_dirs, query_type):
    """
    Load annotation and filter data for a sample.
    We'll use the first available directory since these should be consistent.
    """
    for input_dir in input_dirs:
        # Try to load annotated bed file
        if query_type == 'lr_query':
            annotated_file = os.path.join(input_dir, sample_id, f"{sample_id}.sr.annotated.bed")
            filter_file = os.path.join(input_dir, sample_id, f"{sample_id}.sr.filter_info.tsv")
        else:  # sr_query
            annotated_file = os.path.join(input_dir, sample_id, f"{sample_id}.lr.annotated.bed")
            filter_file = os.path.join(input_dir, sample_id, f"{sample_id}.lr.filter_info.tsv")
        
        if os.path.exists(annotated_file) and os.path.exists(filter_file):
            try:
                # Load context mapping
                annotated_df = pd.read_csv(annotated_file, sep='\t', names=['vid', 'context'])
                context_map = dict(zip(annotated_df['vid'], annotated_df['context']))
                
                # Load filter mapping
                filter_df = pd.read_csv(filter_file, sep='\t')
                filter_map = dict(zip(filter_df['wrapped_variant_id'], filter_df['filter_category']))
                
                return context_map, filter_map
            except Exception as e:
                print(f"Error loading annotation data from {input_dir}: {e}")
                continue
    
    print(f"Warning: Could not load annotation data for {sample_id}")
    return {}, {}

def coalesce_variants(sample_data_dict, overlap_methods):
    """
    Coalesce variants across callers and count how many callers each variant was matched in.
    
    Args:
        sample_data_dict: Dictionary mapping caller -> DataFrame
        overlap_methods: List of overlap methods to analyze (e.g., ['ovr1a', 'ovr1b', ...])
    
    Returns:
        DataFrame with coalesced results
    """
    if not sample_data_dict:
        return pd.DataFrame()
    
    # Get all unique variants across all callers
    all_variants = set()
    for caller, df in sample_data_dict.items():
        all_variants.update(df['clean_VID'].tolist())
    
    coalesced_data = []
    
    for vid in all_variants:
        # Get variant info from first caller that has it
        variant_info = None
        for caller, df in sample_data_dict.items():
            variant_rows = df[df['clean_VID'] == vid]
            if not variant_rows.empty:
                variant_info = variant_rows.iloc[0]
                break
        
        if variant_info is None:
            continue
        
        # Count matches across callers for each overlap method
        match_counts = {}
        for method in overlap_methods:
            count = 0
            for caller, df in sample_data_dict.items():
                variant_rows = df[df['clean_VID'] == vid]
                if not variant_rows.empty:
                    # Check if this specific variant was matched in this caller
                    match_value = variant_rows.iloc[0][method]
                    if match_value != 'NO_OVR':
                        count += 1
            match_counts[method] = count
        
        # Create row for this variant
        row = {
            'VID': variant_info['VID'],
            'clean_VID': vid,
            'chr': variant_info['chr'],
            'start': variant_info['start'],
            'end': variant_info['end'],
            'svtype': variant_info['svtype'],
            'length': variant_info['length'],
            'AF': variant_info['AF'],
            **match_counts
        }
        
        coalesced_data.append(row)
    
    return pd.DataFrame(coalesced_data)

def calculate_coalesced_match_rates(df, group_cols=None, num_callers=2):
    """
    Calculate match rates for coalesced data, showing counts for at least N callers.
    """
    
    def _calc_rates(data):
        total = len(data)
        rates = {}
        
        # Calculate match rates for each overlap method
        for col in ['ovr1a', 'ovr1b', 'ovr2a', 'ovr2b', 'ovr3']:
            if col not in data.columns:
                continue
            
            # Count variants matched in at least N callers (>= N)
            for n in range(1, num_callers + 1):
                matches = len(data[data[col] >= n])
                rates[f'{col}_matches_{n}'] = matches
                rates[f'{col}_rate_{n}'] = f"{(matches/total)*100:.2f}%" if total > 0 else "0.00%"
        
        result = pd.Series({
            'total_variants': total,
            **rates
        })
        
        return result
    
    if group_cols:
        # Select only the necessary columns for grouping and calculation
        calc_cols = group_cols + ['ovr1a', 'ovr1b', 'ovr2a', 'ovr2b', 'ovr3']
        calc_df = df[calc_cols].copy()
        result = calc_df.groupby(group_cols, observed=True).apply(_calc_rates)
        if len(group_cols) == 1:
            result.index.name = 'category'
            return result.reset_index()
        else:
            return result.reset_index()
    else:
        return pd.DataFrame([_calc_rates(df)])

def calculate_combined_coalesced_rates(df, num_callers=2):
    """Calculate coalesced match rates for combinations of SV type, length bucket, context, and filter category."""
    # Group by all four categories
    grouped = df.groupby(['svtype', 'length_bucket', 'context', 'filter_category'])
    
    # Only process groups that have variants
    valid_groups = []
    for name, group in grouped:
        if len(group) > 0:
            rates = calculate_coalesced_match_rates(group, num_callers=num_callers)
            rates['svtype'] = name[0]
            rates['length_bucket'] = name[1]
            rates['context'] = name[2]
            rates['filter_category'] = name[3]
            valid_groups.append(rates)
    
    if not valid_groups:
        return pd.DataFrame()
    
    # Combine all valid groups
    combined_rates = pd.concat(valid_groups, ignore_index=True)
    
    # Create category column combining all four dimensions
    combined_rates['category'] = combined_rates.apply(
        lambda x: f"{x['svtype']}_{x['length_bucket']}_{x['context']}", axis=1
    )
    
    # Reorder columns to put category and filter_category first
    cols = ['category', 'filter_category'] + [col for col in combined_rates.columns 
                          if col not in ['category', 'filter_category', 'svtype', 'length_bucket', 'context']]
    
    return combined_rates[cols]

def analyze_sample(sample_id, input_dirs, query_type, output_dir):
    """Analyze a single sample across multiple callers."""
    
    print(f"Processing sample: {sample_id}")
    
    # Load overlap data for all callers
    sample_data = load_sample_data(sample_id, input_dirs, query_type)
    
    if not sample_data:
        print(f"✗ No data found for sample {sample_id}")
        return False
    
    # Load annotation data
    context_map, filter_map = load_annotation_data(sample_id, input_dirs, query_type)
    
    # Coalesce variants across callers
    overlap_methods = ['ovr1a', 'ovr1b', 'ovr2a', 'ovr2b', 'ovr3']
    coalesced_df = coalesce_variants(sample_data, overlap_methods)
    
    if coalesced_df.empty:
        print(f"✗ No coalesced data for sample {sample_id}")
        return False
    
    # Add annotations
    coalesced_df['context'] = coalesced_df['clean_VID'].map(context_map)
    coalesced_df['filter_category'] = coalesced_df['VID'].map(filter_map)
    
    # Handle missing annotations
    coalesced_df['context'] = coalesced_df['context'].fillna('UNKNOWN')
    coalesced_df['filter_category'] = coalesced_df['filter_category'].fillna('Other')
    
    # Calculate length buckets
    coalesced_df['length_bucket'] = coalesced_df['length'].apply(get_length_bucket)
    
    # Ensure proper ordering of length buckets
    coalesced_df['length_bucket'] = pd.Categorical(
        coalesced_df['length_bucket'],
        categories=['<50', '50-100', '100-500', '500-5000', '5000-10000', '>10000'],
        ordered=True
    )
    
    # Create output directory for this sample
    sample_output_dir = os.path.join(output_dir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)
    
    num_callers = len(input_dirs)
    
    # Create main output file
    output_file = os.path.join(sample_output_dir, f'{query_type}.match_rates.coalesced.tsv')
    with open(output_file, 'w') as f:
        # Overall rates
        overall_rates = calculate_coalesced_match_rates(coalesced_df, num_callers=num_callers)
        overall_rates['category'] = 'overall'
        overall_rates['filter_category'] = 'All'
        
        # SV type rates
        svtype_rates = calculate_coalesced_match_rates(coalesced_df, ['svtype'], num_callers=num_callers)
        svtype_rates['filter_category'] = 'All'
        
        # Length bucket rates
        length_rates = calculate_coalesced_match_rates(coalesced_df, ['length_bucket'], num_callers=num_callers)
        length_rates['filter_category'] = 'All'
        
        # Context rates
        context_rates = calculate_coalesced_match_rates(coalesced_df, ['context'], num_callers=num_callers)
        context_rates['filter_category'] = 'All'
        
        # Filter category rates
        filter_rates = calculate_coalesced_match_rates(coalesced_df, ['filter_category'], num_callers=num_callers)
        if len(filter_rates) > 0:
            filter_rates['filter_category'] = 'All'
        else:
            filter_rates = pd.DataFrame()
        
        # Combined rates with filter stratification
        svtype_filter_rates = calculate_coalesced_match_rates(coalesced_df, ['svtype', 'filter_category'], num_callers=num_callers)
        svtype_filter_rates['category'] = svtype_filter_rates['svtype']
        svtype_filter_rates = svtype_filter_rates.drop('svtype', axis=1)
        
        length_filter_rates = calculate_coalesced_match_rates(coalesced_df, ['length_bucket', 'filter_category'], num_callers=num_callers)
        length_filter_rates['category'] = length_filter_rates['length_bucket']
        length_filter_rates = length_filter_rates.drop('length_bucket', axis=1)
        
        context_filter_rates = calculate_coalesced_match_rates(coalesced_df, ['context', 'filter_category'], num_callers=num_callers)
        context_filter_rates['category'] = context_filter_rates['context']
        context_filter_rates = context_filter_rates.drop('context', axis=1)
        
        # Combine all rates into a single DataFrame
        all_rates = pd.concat([
            overall_rates,
            pd.DataFrame({'category': [''], 'filter_category': ['']}),  # Blank row
            svtype_rates,
            svtype_filter_rates,
            pd.DataFrame({'category': [''], 'filter_category': ['']}),  # Blank row
            length_rates,
            length_filter_rates,
            pd.DataFrame({'category': [''], 'filter_category': ['']}),  # Blank row
            context_rates,
            context_filter_rates,
            pd.DataFrame({'category': [''], 'filter_category': ['']}),  # Blank row
            filter_rates
        ], ignore_index=True)
        
        # Reorder columns to put category and filter_category first
        cols = ['category', 'filter_category'] + [col for col in all_rates.columns if col not in ['category', 'filter_category']]
        all_rates = all_rates[cols]
        
        # Write to file
        all_rates.to_csv(f, sep='\t', index=False)
    
    # Create combined analysis output file
    combined_file = os.path.join(sample_output_dir, f'{query_type}.match_rates.combined.coalesced.tsv')
    combined_rates = calculate_combined_coalesced_rates(coalesced_df, num_callers=num_callers)
    combined_rates.to_csv(combined_file, sep='\t', index=False)
    
    print(f"✓ Sample {sample_id} processing completed successfully")
    return True

def analyze_sample_per_caller(sample_id, input_dirs, query_type, output_dir):
    """Analyze a single sample for each caller individually."""
    
    print(f"Processing sample: {sample_id}")
    
    # Load overlap data for all callers
    sample_data = load_sample_data(sample_id, input_dirs, query_type)
    
    if not sample_data:
        print(f"✗ No data found for sample {sample_id}")
        return False
    
    # Load annotation data
    context_map, filter_map = load_annotation_data(sample_id, input_dirs, query_type)
    
    # Process each caller individually
    for caller_name, caller_df in sample_data.items():
        print(f"  Processing caller: {caller_name}")
        
        # Add annotations
        caller_df['context'] = caller_df['clean_VID'].map(context_map)
        caller_df['filter_category'] = caller_df['VID'].map(filter_map)
        
        # Handle missing annotations
        caller_df['context'] = caller_df['context'].fillna('UNKNOWN')
        caller_df['filter_category'] = caller_df['filter_category'].fillna('Other')
        
        # Calculate length buckets
        caller_df['length_bucket'] = caller_df['length'].apply(get_length_bucket)
        
        # Ensure proper ordering of length buckets
        caller_df['length_bucket'] = pd.Categorical(
            caller_df['length_bucket'],
            categories=['<50', '50-100', '100-500', '500-5000', '5000-10000', '>10000'],
            ordered=True
        )
        
        # Create output directory for this sample
        sample_output_dir = os.path.join(output_dir, sample_id)
        os.makedirs(sample_output_dir, exist_ok=True)
        
        # Calculate single-caller match rates (binary: matched or not)
        def calculate_single_caller_match_rates(df, group_cols=None):
            """Calculate match rates for single caller data."""
            
            def _calc_rates(data):
                total = len(data)
                rates = {}
                
                # Calculate match rates for each overlap method
                for col in ['ovr1a', 'ovr1b', 'ovr2a', 'ovr2b', 'ovr3']:
                    if col not in data.columns:
                        continue
                    
                    # Count variants with matches (not 'NO_OVR')
                    matches = len(data[data[col] != 'NO_OVR'])
                    rates[f'{col}_matches_1'] = matches
                    rates[f'{col}_rate_1'] = f"{(matches/total)*100:.2f}%" if total > 0 else "0.00%"
                
                result = pd.Series({
                    'total_variants': total,
                    **rates
                })
                
                return result
            
            if group_cols:
                # Select only the necessary columns for grouping and calculation
                calc_cols = group_cols + ['ovr1a', 'ovr1b', 'ovr2a', 'ovr2b', 'ovr3']
                calc_df = df[calc_cols].copy()
                result = calc_df.groupby(group_cols, observed=True).apply(_calc_rates)
                if len(group_cols) == 1:
                    result.index.name = 'category'
                    return result.reset_index()
                else:
                    return result.reset_index()
            else:
                return pd.DataFrame([_calc_rates(df)])
        
        def calculate_combined_single_caller_rates(df):
            """Calculate single caller match rates for combinations of SV type, length bucket, context, and filter category."""
            # Group by all four categories
            grouped = df.groupby(['svtype', 'length_bucket', 'context', 'filter_category'])
            
            # Only process groups that have variants
            valid_groups = []
            for name, group in grouped:
                if len(group) > 0:
                    rates = calculate_single_caller_match_rates(group)
                    rates['svtype'] = name[0]
                    rates['length_bucket'] = name[1]
                    rates['context'] = name[2]
                    rates['filter_category'] = name[3]
                    valid_groups.append(rates)
            
            if not valid_groups:
                return pd.DataFrame()
            
            # Combine all valid groups
            combined_rates = pd.concat(valid_groups, ignore_index=True)
            
            # Create category column combining all four dimensions
            combined_rates['category'] = combined_rates.apply(
                lambda x: f"{x['svtype']}_{x['length_bucket']}_{x['context']}", axis=1
            )
            
            # Reorder columns to put category and filter_category first
            cols = ['category', 'filter_category'] + [col for col in combined_rates.columns 
                              if col not in ['category', 'filter_category', 'svtype', 'length_bucket', 'context']]
            
            return combined_rates[cols]
        
        # Create main output file
        output_file = os.path.join(sample_output_dir, f'{query_type}.match_rates.{caller_name}.coalesced.tsv')
        with open(output_file, 'w') as f:
            # Overall rates
            overall_rates = calculate_single_caller_match_rates(caller_df)
            overall_rates['category'] = 'overall'
            overall_rates['filter_category'] = 'All'
            
            # SV type rates
            svtype_rates = calculate_single_caller_match_rates(caller_df, ['svtype'])
            svtype_rates['filter_category'] = 'All'
            
            # Length bucket rates
            length_rates = calculate_single_caller_match_rates(caller_df, ['length_bucket'])
            length_rates['filter_category'] = 'All'
            
            # Context rates
            context_rates = calculate_single_caller_match_rates(caller_df, ['context'])
            context_rates['filter_category'] = 'All'
            
            # Filter category rates
            filter_rates = calculate_single_caller_match_rates(caller_df, ['filter_category'])
            if len(filter_rates) > 0:
                filter_rates['filter_category'] = 'All'
            else:
                filter_rates = pd.DataFrame()
            
            # Combined rates with filter stratification
            svtype_filter_rates = calculate_single_caller_match_rates(caller_df, ['svtype', 'filter_category'])
            svtype_filter_rates['category'] = svtype_filter_rates['svtype']
            svtype_filter_rates = svtype_filter_rates.drop('svtype', axis=1)
            
            length_filter_rates = calculate_single_caller_match_rates(caller_df, ['length_bucket', 'filter_category'])
            length_filter_rates['category'] = length_filter_rates['length_bucket']
            length_filter_rates = length_filter_rates.drop('length_bucket', axis=1)
            
            context_filter_rates = calculate_single_caller_match_rates(caller_df, ['context', 'filter_category'])
            context_filter_rates['category'] = context_filter_rates['context']
            context_filter_rates = context_filter_rates.drop('context', axis=1)
            
            # Combine all rates into a single DataFrame
            all_rates = pd.concat([
                overall_rates,
                pd.DataFrame({'category': [''], 'filter_category': ['']}),  # Blank row
                svtype_rates,
                svtype_filter_rates,
                pd.DataFrame({'category': [''], 'filter_category': ['']}),  # Blank row
                length_rates,
                length_filter_rates,
                pd.DataFrame({'category': [''], 'filter_category': ['']}),  # Blank row
                context_rates,
                context_filter_rates,
                pd.DataFrame({'category': [''], 'filter_category': ['']}),  # Blank row
                filter_rates
            ], ignore_index=True)
            
            # Reorder columns to put category and filter_category first
            cols = ['category', 'filter_category'] + [col for col in all_rates.columns if col not in ['category', 'filter_category']]
            all_rates = all_rates[cols]
            
            # Write to file
            all_rates.to_csv(f, sep='\t', index=False)
        
        # Create combined analysis output file
        combined_file = os.path.join(sample_output_dir, f'{query_type}.match_rates.combined.{caller_name}.coalesced.tsv')
        combined_rates = calculate_combined_single_caller_rates(caller_df)
        combined_rates.to_csv(combined_file, sep='\t', index=False)
    
    print(f"✓ Sample {sample_id} processing completed successfully")
    return True

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input-dirs', required=True, nargs='+',
                      help='Input directories containing caller results (e.g., sniffles_results cutesv_results)')
    parser.add_argument('--mapping', required=True,
                      help='Mapping file with sample IDs')
    parser.add_argument('--query-type', required=True, choices=['lr_query', 'sr_query'],
                      help='Type of query to analyze')
    parser.add_argument('--output-dir', required=True,
                      help='Output directory for coalesced results')
    parser.add_argument('--num-samples', type=int, default=None,
                      help='Maximum number of samples to process (default: all)')
    parser.add_argument('--per-caller', action='store_true',
                      help='Analyze each caller individually instead of coalescing across callers')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Read mapping file
    mapping_df = pd.read_csv(args.mapping, sep='\t')
    
    # Determine number of samples to process
    total_samples = len(mapping_df)
    samples_to_process = args.num_samples if args.num_samples else total_samples
    
    print("")
    print("=" * 60)
    if args.per_caller:
        print("PER-CALLER MATCH RATE ANALYSIS")
    else:
        print("MULTI-CALLER MATCH RATE ANALYSIS")
    print("=" * 60)
    print(f"Query type: {args.query_type}")
    print(f"Input directories: {', '.join(args.input_dirs)}")
    print(f"Output directory: {args.output_dir}")
    print(f"Samples to process: {samples_to_process}")
    print("")
    
    # Process samples
    successful_samples = 0
    failed_samples = 0
    total_processed = 0
    
    for _, row in mapping_df.iterrows():
        if args.num_samples and total_processed >= args.num_samples:
            break
        
        sr_id = row['SR_ID']
        total_processed += 1
        
        print(f"SAMPLE {total_processed}/{samples_to_process}")
        
        if args.per_caller:
            success = analyze_sample_per_caller(sr_id, args.input_dirs, args.query_type, args.output_dir)
        else:
            success = analyze_sample(sr_id, args.input_dirs, args.query_type, args.output_dir)
        
        if success:
            successful_samples += 1
        else:
            failed_samples += 1
        
        print("")
    
    print("=" * 60)
    print("PROCESSING COMPLETE")
    print("=" * 60)
    print(f"Successfully processed: {successful_samples} samples")
    print(f"Failed to process: {failed_samples} samples")
    print("=" * 60)

if __name__ == '__main__':
    main() 