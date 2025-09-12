#!/Users/kjaising/miniconda3/envs/venv/bin/python3

"""
Analyze match rates between LR and SR callsets.
Generates tables for overall match rates, rates by SV type, by length bucket, and by genomic context.
Also generates a combined table showing rates by SV type x length bucket x genomic context.
"""

import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import os

def clean_vid(vid):
    # Remove the prefix and suffix "__" from variant IDs
    if vid.startswith('__') and vid.endswith('__'):
        return vid[2:-2]
    return vid

def calculate_match_rates(df, group_cols=None):
    """Calculate match rates for each overlap method."""
    
    def _calc_rates(data):
        total = len(data)
        rates = {}
        
        # Calculate match rates for each overlap method
        for col in ['ovr1a', 'ovr1b', 'ovr2a', 'ovr2b', 'ovr3']:
            matches = len(data[data[col] != 'NO_OVR'])
            rates[f'{col}_matches'] = matches
            rates[f'{col}_rate'] = f"{(matches/total)*100:.2f}%"
        
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

def calculate_combined_match_rates(df):
    """Calculate match rates for combinations of SV type, length bucket, context, and filter category."""
    # Group by all four categories
    grouped = df.groupby(['svtype', 'length_bucket', 'context', 'filter_category'])
    
    # Only process groups that have variants
    valid_groups = []
    for name, group in grouped:
        if len(group) > 0:
            rates = calculate_match_rates(group)
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

def analyze_callset(query_bed, annotated_bed, filter_info_file, output_prefix):
    """Analyze match rates for a callset."""
    # Read input files with proper handling of comments and column names
    query_df = pd.read_csv(query_bed, sep='\t', comment='#', header=None,
                          names=['chr', 'start', 'end', 'VID', 'svtype', 'length', 'AF', 
                                'ovr1a', 'ovr1b', 'ovr2a', 'ovr2b', 'ovr3'])
    annotated_df = pd.read_csv(annotated_bed, sep='\t', names=['vid', 'context'])
    
    # Map wrapped VIDs in query file to unwrapped VIDs for annotation lookup
    # Query file VIDs have format __original_vid__, annotation file has original_vid
    query_df['clean_VID'] = query_df['VID'].apply(clean_vid)
    
    # Create context mapping using cleaned VIDs
    context_map = dict(zip(annotated_df['vid'], annotated_df['context']))
    
    # Load filter info if provided
    filter_map = {}
    if filter_info_file:
        filter_df = pd.read_csv(filter_info_file, sep='\t')
        filter_map = dict(zip(filter_df['wrapped_variant_id'], filter_df['filter_category']))
    
    # Map context and filter info to query variants
    query_df['context'] = query_df['clean_VID'].map(context_map)
    query_df['filter_category'] = query_df['VID'].map(filter_map)
    
    # Calculate mapping success
    n_context_mapped = query_df['context'].notna().sum()
    n_filter_mapped = query_df['filter_category'].notna().sum()
    
    print(f"Mapped context for {n_context_mapped} out of {len(query_df)} variants")
    print(f"Mapped filter info for {n_filter_mapped} out of {len(query_df)} variants")
    
    if n_filter_mapped == 0:
        print("WARNING: No filter mapping found! Defaulting to 'Other'")
    
    # If no context mapping worked, something is wrong
    if query_df['context'].isna().all():
        print("Warning: No context mapping found! This should not happen.")
        query_df['context'] = 'UNKNOWN'
    
    # If no filter mapping worked, assign default
    if query_df['filter_category'].isna().all():
        print("Warning: No filter mapping found! Using 'Other' as default.")
        query_df['filter_category'] = 'Other'
    else:
        # Fill missing filter categories with 'Other'
        query_df['filter_category'] = query_df['filter_category'].fillna('Other')
    
    # Calculate length buckets
    query_df['length_bucket'] = query_df['length'].apply(get_length_bucket)
    
    # Ensure proper ordering of length buckets
    query_df['length_bucket'] = pd.Categorical(
        query_df['length_bucket'],
        categories=['<50', '50-100', '100-500', '500-5000', '5000-10000', '>10000'],
        ordered=True
    )
    
    # Create main output file
    with open(output_prefix + '.match_rates.tsv', 'w') as f:
        # Overall rates - add 'overall' as category
        overall_rates = calculate_match_rates(query_df)
        overall_rates['category'] = 'overall'
        overall_rates['filter_category'] = 'All'
        
        # SV type rates
        svtype_rates = calculate_match_rates(query_df, ['svtype'])
        svtype_rates['filter_category'] = 'All'
        
        # Length bucket rates
        length_rates = calculate_match_rates(query_df, ['length_bucket'])
        length_rates['filter_category'] = 'All'
        
        # Context rates
        context_rates = calculate_match_rates(query_df, ['context'])
        context_rates['filter_category'] = 'All'
        
        # Filter category rates
        filter_rates = calculate_match_rates(query_df, ['filter_category'])
        # When grouping by 'filter_category', the result has 'category' column containing the filter values
        # We need to rename this to match our expected structure
        if len(filter_rates) > 0:
            # The 'category' column contains the filter categories
            # Keep 'category' as is and add 'filter_category' = 'All'
            filter_rates['filter_category'] = 'All'
        else:
            # If no data, create empty rates
            filter_rates = pd.DataFrame()
        
        # Combined rates with filter stratification
        svtype_filter_rates = calculate_match_rates(query_df, ['svtype', 'filter_category'])
        svtype_filter_rates['category'] = svtype_filter_rates['svtype']
        svtype_filter_rates = svtype_filter_rates.drop('svtype', axis=1)
        
        length_filter_rates = calculate_match_rates(query_df, ['length_bucket', 'filter_category'])
        length_filter_rates['category'] = length_filter_rates['length_bucket']
        length_filter_rates = length_filter_rates.drop('length_bucket', axis=1)
        
        context_filter_rates = calculate_match_rates(query_df, ['context', 'filter_category'])
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
    combined_rates = calculate_combined_match_rates(query_df)
    combined_rates.to_csv(output_prefix + '.match_rates.combined.tsv', sep='\t', index=False)

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--query-bed', required=True,
                      help='Query BED file (either LR or SR)')
    parser.add_argument('--annotated-bed', required=True,
                      help='Annotated BED file for context mapping')
    parser.add_argument('--filter-info', required=True,
                      help='Filter information file')
    parser.add_argument('--output-prefix', required=True,
                      help='Output file prefix')
    parser.add_argument('--query-type', required=True, choices=['lr_query', 'sr_query'],
                      help='Type of query (lr_query or sr_query)')
    args = parser.parse_args()

    print(f"Analyzing {args.query_type} match rates...")
    analyze_callset(
        args.query_bed,
        args.annotated_bed,
        args.filter_info,
        args.output_prefix
    )

    print("Match rate analysis completed!")

if __name__ == '__main__':
    main() 