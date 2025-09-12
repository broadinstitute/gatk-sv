#!/Users/kjaising/miniconda3/envs/venv/bin/python3
"""
Analyze multi-caller support rates across multiple samples.
Generates violin plots and heatmap tables for SV support analysis across N callers.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

def parse_rate(rate_str):
    """Convert percentage string to float"""
    if pd.isna(rate_str) or rate_str == '':
        return np.nan
    return float(rate_str.strip('%'))

def parse_category(category):
    """Parse combined category like 'DEL_50-100_RM' into components"""
    if '_' not in category:
        return None, None, category
    
    parts = category.split('_')
    if len(parts) >= 3:
        svtype = parts[0]
        size_bucket = parts[1]
        region = parts[2]
        return svtype, size_bucket, region
    elif len(parts) == 2:
        # Handle cases like 'DEL_PASS+MULTIALLELIC' 
        return parts[0], None, parts[1]
    else:
        return None, None, category

def discover_samples(input_dir, mapping_file):
    """
    Discover all samples with results in the input directory
    
    Args:
        input_dir: Directory containing sample subdirectories
        mapping_file: Path to mapping file with SR_ID and LR_ID columns
    
    Returns:
        List of (sr_id, lr_id) tuples for available samples
    """
    # Read mapping file to get the mapping
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    
    # Find all subdirectories that have the required result files
    available_samples = []
    
    for _, row in mapping_df.iterrows():
        sr_id = row['SR_ID']
        lr_id = row['LR_ID']
        
        # Check if sample directory exists and has the required file
        sample_dir = os.path.join(input_dir, sr_id)
        result_file = os.path.join(sample_dir, 'lr_query.match_rates.combined.coalesced.tsv')
        
        if os.path.exists(result_file):
            available_samples.append((sr_id, lr_id))
    
    return available_samples

def collect_sample_data(mapping_file, input_dir, max_n_matches, svtypes_filter=None, use_combined=True):
    """
    Collect multi-caller match rates data from all available samples
    
    Args:
        mapping_file: Path to mapping file with SR_ID and LR_ID columns
        input_dir: Directory containing multi-caller results
        max_n_matches: Maximum number of caller matches to collect (1, 2, 3, etc.)
        svtypes_filter: List of SV types to include (e.g., ['DEL', 'DUP', 'INS'])
        use_combined: Whether to use combined.coalesced.tsv files (default) or regular .coalesced.tsv files
    
    Returns:
        DataFrame with all sample data
    """
    
    # Discover available samples
    available_samples = discover_samples(input_dir, mapping_file)
    
    if not available_samples:
        raise ValueError("No samples with results found in the input directory")
    
    print(f"Found {len(available_samples)} samples with results")
    
    all_data = []
    processed = 0
    
    for sr_id, lr_id in available_samples:
        # Determine which file to use
        if use_combined:
            file_path = os.path.join(input_dir, sr_id, 'lr_query.match_rates.combined.coalesced.tsv')
        else:
            file_path = os.path.join(input_dir, sr_id, 'lr_query.match_rates.coalesced.tsv')
        
        if not os.path.exists(file_path):
            print(f"Warning: File not found for sample {sr_id}: {file_path}")
            continue
            
        try:
            # Read the file
            sample_df = pd.read_csv(file_path, sep='\t')
            sample_df['sample_id'] = sr_id
            sample_df['lr_id'] = lr_id
            all_data.append(sample_df)
            processed += 1
            print(f"Processed sample {sr_id} ({processed}/{len(available_samples)})")
            
        except Exception as e:
            print(f"Error processing sample {sr_id}: {e}")
            continue
    
    if not all_data:
        raise ValueError("No sample data was successfully loaded")
    
    # Combine all data
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Add parsed columns for combined data
    parsed_data = []
    for _, row in combined_df.iterrows():
        svtype, size_bucket, region = parse_category(row['category'])
        parsed_row = row.copy()
        parsed_row['svtype'] = svtype
        parsed_row['size_bucket'] = size_bucket  
        parsed_row['region'] = region
        parsed_data.append(parsed_row)
    
    result_df = pd.DataFrame(parsed_data)
    
    # Apply SV type filtering if specified
    if svtypes_filter:
        initial_rows = len(result_df)
        result_df = result_df[result_df['svtype'].isin(svtypes_filter)]
        filtered_rows = len(result_df)
        print(f"SV type filtering: kept {filtered_rows} rows out of {initial_rows} (types: {', '.join(svtypes_filter)})")
    
    # Parse rate columns for all n_matches levels
    for n in range(1, max_n_matches + 1):
        rate_col = f'ovr1a_rate_{n}'
        matches_col = f'ovr1a_matches_{n}'
        
        if rate_col in result_df.columns:
            result_df[f'match_rate_{n}'] = result_df[rate_col].apply(parse_rate)
            result_df[f'match_count_{n}'] = result_df[matches_col]
    
    return result_df

def create_violin_plots(df, output_dir, max_n_matches, svtypes_filter=None):
    """Create violin plots for different category combinations with multiple n_matches subplots"""
    
    # Filter for PASS+MULTIALLELIC only
    pass_multi_df = df[df['filter_category'] == 'PASS+MULTIALLELIC'].copy()
    
    if pass_multi_df.empty:
        print("Warning: No PASS+MULTIALLELIC data found")
        return {}
    
    # Create figure with subplots for each n_matches level
    fig, axes = plt.subplots(1, max_n_matches, figsize=(6 * max_n_matches, 8))
    if max_n_matches == 1:
        axes = [axes]  # Make it a list for consistency
    
    all_stats = {}
    all_data_for_ylim = []  # Collect all data to determine common y-axis limits
    
    # First pass: collect all data and calculate stats
    for n_matches in range(1, max_n_matches + 1):
        # Check if we have the required columns
        rate_col = f'match_rate_{n_matches}'
        if rate_col not in pass_multi_df.columns:
            print(f"Warning: Column {rate_col} not found, skipping {n_matches}+ callers")
            continue
        
        # 1. Calculate overall rates (all PASS+MULTIALLELIC calls)
        overall_rates = []
        for sample in pass_multi_df['sample_id'].unique():
            sample_data = pass_multi_df[pass_multi_df['sample_id'] == sample]
            
            # Calculate weighted average rate
            total_variants = sample_data['total_variants'].sum()
            total_matches = (sample_data['total_variants'] * sample_data[rate_col] / 100).sum()
            
            if total_variants > 0:
                overall_rate = (total_matches / total_variants) * 100
                overall_rates.append(overall_rate)
        
        # 2. Calculate US/RM regions rates
        us_rm_categories = ['US', 'RM']
        us_rm_data = []
        
        for sample in pass_multi_df['sample_id'].unique():
            sample_data = pass_multi_df[pass_multi_df['sample_id'] == sample]
            
            # Filter for US/RM regions
            region_data = sample_data[sample_data['region'].isin(us_rm_categories)]
            
            if not region_data.empty:
                total_variants = region_data['total_variants'].sum()
                total_matches = (region_data['total_variants'] * region_data[rate_col] / 100).sum()
                
                if total_variants > 0:
                    rate = (total_matches / total_variants) * 100
                    us_rm_data.append(rate)
        
        # 3. Calculate SD/SR regions rates
        sd_sr_categories = ['SD', 'SR']
        sd_sr_data = []
        
        for sample in pass_multi_df['sample_id'].unique():
            sample_data = pass_multi_df[pass_multi_df['sample_id'] == sample]
            
            # Filter for SD/SR regions
            region_data = sample_data[sample_data['region'].isin(sd_sr_categories)]
            
            if not region_data.empty:
                total_variants = region_data['total_variants'].sum()
                total_matches = (region_data['total_variants'] * region_data[rate_col] / 100).sum()
                
                if total_variants > 0:
                    rate = (total_matches / total_variants) * 100
                    sd_sr_data.append(rate)
        
        # Store stats for this n_matches level
        all_stats[n_matches] = {
            'overall': overall_rates,
            'us_rm': us_rm_data,
            'sd_sr': sd_sr_data
        }
        
        # Collect all data for y-axis calculation
        all_data_for_ylim.extend(overall_rates)
        all_data_for_ylim.extend(us_rm_data)
        all_data_for_ylim.extend(sd_sr_data)
    
    # Use full 0-100% scale for consistency across plots
    common_ylim = (0, 100)
    
    # Second pass: create plots with common y-axis
    for n_idx, n_matches in enumerate(range(1, max_n_matches + 1)):
        ax = axes[n_idx]
        
        if n_matches not in all_stats:
            continue
        
        stats = all_stats[n_matches]
        
        # Create violin plot for this n_matches level
        data_to_plot = []
        labels = []
        positions = []
        
        if stats['overall']:
            data_to_plot.append(stats['overall'])
            labels.append('Overall')
            positions.append(1)
        
        if stats['us_rm']:
            data_to_plot.append(stats['us_rm'])
            labels.append('US/RM')
            positions.append(2)
        
        if stats['sd_sr']:
            data_to_plot.append(stats['sd_sr'])
            labels.append('SD/SR')
            positions.append(3)
        
        if data_to_plot:
            # Create violin plot
            violin_parts = ax.violinplot(data_to_plot, positions=positions, widths=0.6, 
                                        showmeans=True, showmedians=True, showextrema=True)
            
            # Customize violin plot colors
            colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # Blue, Orange, Green
            for i, pc in enumerate(violin_parts['bodies']):
                if i < len(colors):
                    pc.set_facecolor(colors[i])
                    pc.set_alpha(0.7)
            
            # Style the subplot
            ax.set_ylabel('Support Rate (%)', fontsize=12)
            ax.set_xlabel('Genomic Context', fontsize=12)
            
            # Create subplot title
            if n_matches == 1:
                title = '1+ Callers'
            else:
                title = f'{n_matches} Callers'
            ax.set_title(title, fontsize=16, fontweight='bold')
            
            ax.set_xticks(positions)
            ax.set_xticklabels(labels, fontsize=9)
            ax.grid(True, alpha=0.3, axis='y')
            
            # Set common y-axis limits
            ax.set_ylim(common_ylim)
            
            # Add sample counts as text
            for i, (pos, data) in enumerate(zip(positions, data_to_plot)):
                y_min, y_max = ax.get_ylim()
                ax.text(pos, y_min + (y_max - y_min) * 0.02, 
                       f'n={len(data)}', ha='center', va='bottom', fontsize=8, 
                       color='gray', fontweight='bold')
    
    # No overall figure title - removed suptitle
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'violin_combined_support_{max_n_matches}_callers.png'), 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Generated violin plots in {output_dir}")
    for n_matches in range(1, max_n_matches + 1):
        if n_matches in all_stats:
            stats = all_stats[n_matches]
            print(f"  {n_matches}+ callers - Overall: {len(stats['overall'])} samples, US/RM: {len(stats['us_rm'])} samples, SD/SR: {len(stats['sd_sr'])} samples")
    
    return all_stats

def create_heatmap_tables(df, output_dir, max_n_matches, svtypes_filter=None):
    """Create heatmap tables for US/RM and SD/SR regions with multiple n_matches subplots"""
    
    # Filter for PASS+MULTIALLELIC only and data with proper parsing
    pass_multi_df = df[(df['filter_category'] == 'PASS+MULTIALLELIC') & 
                       (df['svtype'].notna()) & 
                       (df['size_bucket'].notna()) & 
                       (df['region'].notna())].copy()
    
    if pass_multi_df.empty:
        print("Warning: No PASS+MULTIALLELIC data found with proper category parsing")
        return
    
    # Define size bucket order
    size_order = ['50-100', '100-500', '500-5000', '5000-10000', '>10000']
    
    # Filter svtype_order based on svtypes_filter if provided
    default_svtype_order = ['DEL', 'DUP', 'INS', 'INV', 'CNV', 'CPX']
    if svtypes_filter:
        svtype_order = [svtype for svtype in default_svtype_order if svtype in svtypes_filter]
    else:
        svtype_order = default_svtype_order
    
    def create_heatmap_for_regions(regions, title_suffix, filename_suffix):
        """Create heatmap for specific regions with multiple n_matches subplots"""
        
        region_data = pass_multi_df[pass_multi_df['region'].isin(regions)]
        
        if region_data.empty:
            print(f"Warning: No data found for regions {regions}")
            return
        
        # Create figure with subplots for each n_matches level
        # Calculate figure size to make cells square
        n_rows = len(svtype_order) if svtypes_filter else len(default_svtype_order)
        n_cols = len(size_order)
        
        # Each cell should be roughly square, aim for ~0.8 inches per cell
        cell_size = 0.8
        subplot_width = n_cols * cell_size + 2  # Add space for labels
        subplot_height = n_rows * cell_size + 1.5  # Add space for labels
        
        fig, axes = plt.subplots(1, max_n_matches, figsize=(subplot_width * max_n_matches, subplot_height))
        if max_n_matches == 1:
            axes = [axes]  # Make it a list for consistency
        
        for n_idx, n_matches in enumerate(range(1, max_n_matches + 1)):
            ax = axes[n_idx]
            
            rate_col = f'match_rate_{n_matches}'
            if rate_col not in region_data.columns:
                print(f"Warning: Column {rate_col} not found, skipping {n_matches}+ callers")
                continue
            
            # Aggregate data across samples
            agg_data = region_data.groupby(['svtype', 'size_bucket']).agg({
                'total_variants': 'sum',
                rate_col: lambda x: np.average(x, weights=region_data.loc[x.index, 'total_variants'])
            }).reset_index()
            
            # Create pivot table
            pivot_data = agg_data.pivot(index='svtype', columns='size_bucket', values=rate_col)
            pivot_counts = agg_data.pivot(index='svtype', columns='size_bucket', values='total_variants')
            
            # Reorder columns and rows
            available_sizes = [size for size in size_order if size in pivot_data.columns]
            available_svtypes = [svtype for svtype in svtype_order if svtype in pivot_data.index]
            
            if available_sizes and available_svtypes:
                pivot_data = pivot_data.reindex(index=available_svtypes, columns=available_sizes)
                pivot_counts = pivot_counts.reindex(index=available_svtypes, columns=available_sizes)
                
                # Create labels with percentages and counts
                labels = pivot_data.copy().astype(str)
                for i in range(len(labels.index)):
                    for j in range(len(labels.columns)):
                        rate = pivot_data.iloc[i, j]
                        count = pivot_counts.iloc[i, j]
                        if pd.notna(rate) and pd.notna(count):
                            labels.iloc[i, j] = f"{rate:.1f}%\n({int(count)})"
                        else:
                            labels.iloc[i, j] = "N/A"
                
                # Create heatmap with square cells
                sns.heatmap(pivot_data, 
                           annot=labels, 
                           fmt='', 
                           cmap='YlOrRd', 
                           cbar_kws={'label': 'Support Rate (%)'},
                           vmin=0, 
                           vmax=100,
                           ax=ax,
                           square=True)  # Make cells square
                
                # Create subplot title
                if n_matches == 1:
                    title = '1+ Callers'
                else:
                    title = f'{n_matches} Callers'
                ax.set_title(title, fontsize=16, fontweight='bold')
                
                ax.set_xlabel('SV Length (bp)', fontsize=12)
                ax.set_ylabel('SV Type', fontsize=12)
                
                # Save individual data table
                pivot_data.to_csv(os.path.join(output_dir, f'support_rates_{filename_suffix}_{n_matches}_callers.csv'))
            else:
                print(f"Warning: No valid data for heatmap {title_suffix} ({n_matches}+ callers)")
        
        # No overall figure title - removed suptitle
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'heatmap_{filename_suffix}_{max_n_matches}_callers.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Generated heatmap for {title_suffix} (1-{max_n_matches} callers)")
    
    # Create heatmaps for both region groups
    create_heatmap_for_regions(['US', 'RM'], 'US/RM Regions', 'us_rm_regions')
    create_heatmap_for_regions(['SD', 'SR'], 'SD/SR Regions', 'sd_sr_regions')

def discover_samples_per_caller(input_dir, mapping_file, caller_name):
    """
    Discover all samples with results for a specific caller in the input directory
    
    Args:
        input_dir: Directory containing sample subdirectories
        mapping_file: Path to mapping file with SR_ID and LR_ID columns
        caller_name: Name of the caller (e.g., 'sniffles', 'cutesv')
    
    Returns:
        List of (sr_id, lr_id) tuples for available samples for this caller
    """
    # Read mapping file to get the mapping
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    
    # Find all subdirectories that have the required result files for this caller
    available_samples = []
    
    for _, row in mapping_df.iterrows():
        sr_id = row['SR_ID']
        lr_id = row['LR_ID']
        
        # Check if sample directory exists and has the required file for this caller
        sample_dir = os.path.join(input_dir, sr_id)
        result_file = os.path.join(sample_dir, f'lr_query.match_rates.combined.{caller_name}.coalesced.tsv')
        
        if os.path.exists(result_file):
            available_samples.append((sr_id, lr_id))
    
    return available_samples

def collect_sample_data_per_caller(mapping_file, input_dir, caller_name, svtypes_filter=None):
    """
    Collect single-caller match rates data from all available samples for a specific caller
    
    Args:
        mapping_file: Path to mapping file with SR_ID and LR_ID columns
        input_dir: Directory containing single-caller results
        caller_name: Name of the caller (e.g., 'sniffles', 'cutesv')
        svtypes_filter: List of SV types to include (e.g., ['DEL', 'DUP', 'INS'])
    
    Returns:
        DataFrame with all sample data for this caller
    """
    
    # Discover available samples for this caller
    available_samples = discover_samples_per_caller(input_dir, mapping_file, caller_name)

    print(input_dir)
    
    if not available_samples:
        raise ValueError(f"No samples with results found for caller {caller_name} in the input directory")
    
    print(f"Found {len(available_samples)} samples with results for caller {caller_name}")
    
    all_data = []
    processed = 0
    
    for sr_id, lr_id in available_samples:
        # Use the caller-specific combined file
        file_path = os.path.join(input_dir, sr_id, f'lr_query.match_rates.combined.{caller_name}.coalesced.tsv')
        
        if not os.path.exists(file_path):
            print(f"Warning: File not found for sample {sr_id}, caller {caller_name}: {file_path}")
            continue
            
        try:
            # Read the file
            sample_df = pd.read_csv(file_path, sep='\t')
            sample_df['sample_id'] = sr_id
            sample_df['lr_id'] = lr_id
            sample_df['caller'] = caller_name
            all_data.append(sample_df)
            processed += 1
            print(f"Processed sample {sr_id} for caller {caller_name} ({processed}/{len(available_samples)})")
            
        except Exception as e:
            print(f"Error processing sample {sr_id} for caller {caller_name}: {e}")
            continue
    
    if not all_data:
        raise ValueError(f"No sample data was successfully loaded for caller {caller_name}")
    
    # Combine all data
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Add parsed columns for combined data
    parsed_data = []
    for _, row in combined_df.iterrows():
        svtype, size_bucket, region = parse_category(row['category'])
        parsed_row = row.copy()
        parsed_row['svtype'] = svtype
        parsed_row['size_bucket'] = size_bucket  
        parsed_row['region'] = region
        parsed_data.append(parsed_row)
    
    result_df = pd.DataFrame(parsed_data)
    
    # Apply SV type filtering if specified
    if svtypes_filter:
        initial_rows = len(result_df)
        result_df = result_df[result_df['svtype'].isin(svtypes_filter)]
        filtered_rows = len(result_df)
        print(f"SV type filtering: kept {filtered_rows} rows out of {initial_rows} (types: {', '.join(svtypes_filter)})")
    
    # Parse rate columns (for per-caller, we only have _1 rates since it's binary)
    rate_col = 'ovr1a_rate_1'
    matches_col = 'ovr1a_matches_1'
    
    if rate_col in result_df.columns:
        result_df['match_rate_1'] = result_df[rate_col].apply(parse_rate)
        result_df['match_count_1'] = result_df[matches_col]
    
    return result_df

def create_violin_plots_per_caller(df, output_dir, caller_name, svtypes_filter=None):
    """Create violin plots for a specific caller"""
    
    # Filter for PASS+MULTIALLELIC only
    pass_multi_df = df[df['filter_category'] == 'PASS+MULTIALLELIC'].copy()
    
    if pass_multi_df.empty:
        print(f"Warning: No PASS+MULTIALLELIC data found for caller {caller_name}")
        return {}
    
    # Create figure with single subplot since we only have binary match/no-match for single caller
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    
    # Calculate rates for single caller (binary: matched or not matched)
    rate_col = 'match_rate_1'
    if rate_col not in pass_multi_df.columns:
        print(f"Warning: Column {rate_col} not found for caller {caller_name}")
        return {}
    
    # 1. Calculate overall rates (all PASS+MULTIALLELIC calls)
    overall_rates = []
    for sample in pass_multi_df['sample_id'].unique():
        sample_data = pass_multi_df[pass_multi_df['sample_id'] == sample]
        
        # Calculate weighted average rate
        total_variants = sample_data['total_variants'].sum()
        total_matches = (sample_data['total_variants'] * sample_data[rate_col] / 100).sum()
        
        if total_variants > 0:
            overall_rate = (total_matches / total_variants) * 100
            overall_rates.append(overall_rate)
    
    # 2. Calculate US/RM regions rates
    us_rm_categories = ['US', 'RM']
    us_rm_data = []
    
    for sample in pass_multi_df['sample_id'].unique():
        sample_data = pass_multi_df[pass_multi_df['sample_id'] == sample]
        
        # Filter for US/RM regions
        region_data = sample_data[sample_data['region'].isin(us_rm_categories)]
        
        if not region_data.empty:
            total_variants = region_data['total_variants'].sum()
            total_matches = (region_data['total_variants'] * region_data[rate_col] / 100).sum()
            
            if total_variants > 0:
                rate = (total_matches / total_variants) * 100
                us_rm_data.append(rate)
    
    # 3. Calculate SD/SR regions rates
    sd_sr_categories = ['SD', 'SR']
    sd_sr_data = []
    
    for sample in pass_multi_df['sample_id'].unique():
        sample_data = pass_multi_df[pass_multi_df['sample_id'] == sample]
        
        # Filter for SD/SR regions
        region_data = sample_data[sample_data['region'].isin(sd_sr_categories)]
        
        if not region_data.empty:
            total_variants = region_data['total_variants'].sum()
            total_matches = (region_data['total_variants'] * region_data[rate_col] / 100).sum()
            
            if total_variants > 0:
                rate = (total_matches / total_variants) * 100
                sd_sr_data.append(rate)
    
    # Create violin plot
    data_to_plot = []
    labels = []
    positions = []
    
    if overall_rates:
        data_to_plot.append(overall_rates)
        labels.append('Overall')
        positions.append(1)
    
    if us_rm_data:
        data_to_plot.append(us_rm_data)
        labels.append('US/RM')
        positions.append(2)
    
    if sd_sr_data:
        data_to_plot.append(sd_sr_data)
        labels.append('SD/SR')
        positions.append(3)
    
    stats = {
        'overall': overall_rates,
        'us_rm': us_rm_data,
        'sd_sr': sd_sr_data
    }
    
    if data_to_plot:
        # Create violin plot
        violin_parts = ax.violinplot(data_to_plot, positions=positions, widths=0.6, 
                                    showmeans=True, showmedians=True, showextrema=True)
        
        # Customize violin plot colors
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # Blue, Orange, Green
        for i, pc in enumerate(violin_parts['bodies']):
            if i < len(colors):
                pc.set_facecolor(colors[i])
                pc.set_alpha(0.7)
        
        # Style the plot
        ax.set_ylabel('Support Rate (%)', fontsize=12)
        ax.set_xlabel('Genomic Context', fontsize=12)
        ax.set_title(f'{caller_name.title()} Caller Support', fontsize=16, fontweight='bold')
        
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=9)
        ax.grid(True, alpha=0.3, axis='y')
        
        # Set y-axis limits to full 0-100%
        ax.set_ylim(0, 100)
        
        # Add sample counts as text
        for i, (pos, data) in enumerate(zip(positions, data_to_plot)):
            y_min, y_max = ax.get_ylim()
            ax.text(pos, y_min + (y_max - y_min) * 0.02, 
                   f'n={len(data)}', ha='center', va='bottom', fontsize=8, 
                   color='gray', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'violin_support_{caller_name}.png'), 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Generated violin plot for {caller_name} in {output_dir}")
    print(f"  Overall: {len(stats['overall'])} samples, US/RM: {len(stats['us_rm'])} samples, SD/SR: {len(stats['sd_sr'])} samples")
    
    return stats

def create_heatmap_tables_per_caller(df, output_dir, caller_name, svtypes_filter=None):
    """Create heatmap tables for US/RM and SD/SR regions for a specific caller"""
    
    # Filter for PASS+MULTIALLELIC only and data with proper parsing
    pass_multi_df = df[(df['filter_category'] == 'PASS+MULTIALLELIC') & 
                       (df['svtype'].notna()) & 
                       (df['size_bucket'].notna()) & 
                       (df['region'].notna())].copy()
    
    if pass_multi_df.empty:
        print(f"Warning: No PASS+MULTIALLELIC data found with proper category parsing for caller {caller_name}")
        return
    
    # Define size bucket order
    size_order = ['50-100', '100-500', '500-5000', '5000-10000', '>10000']
    
    # Filter svtype_order based on svtypes_filter if provided
    default_svtype_order = ['DEL', 'DUP', 'INS', 'INV', 'CNV', 'CPX']
    if svtypes_filter:
        svtype_order = [svtype for svtype in default_svtype_order if svtype in svtypes_filter]
    else:
        svtype_order = default_svtype_order
    
    def create_heatmap_for_regions(regions, title_suffix, filename_suffix):
        """Create heatmap for specific regions for a single caller"""
        
        region_data = pass_multi_df[pass_multi_df['region'].isin(regions)]
        
        if region_data.empty:
            print(f"Warning: No data found for regions {regions} for caller {caller_name}")
            return
        
        # Calculate figure size to make cells square
        n_rows = len(svtype_order) if svtypes_filter else len(default_svtype_order)
        n_cols = len(size_order)
        
        # Each cell should be roughly square, aim for ~0.8 inches per cell
        cell_size = 0.8
        subplot_width = n_cols * cell_size + 2  # Add space for labels
        subplot_height = n_rows * cell_size + 1.5  # Add space for labels
        
        fig, ax = plt.subplots(1, 1, figsize=(subplot_width, subplot_height))
        
        rate_col = 'match_rate_1'  # Single caller only has binary rates
        if rate_col not in region_data.columns:
            print(f"Warning: Column {rate_col} not found for caller {caller_name}")
            return
        
        # Aggregate data across samples
        agg_data = region_data.groupby(['svtype', 'size_bucket']).agg({
            'total_variants': 'sum',
            rate_col: lambda x: np.average(x, weights=region_data.loc[x.index, 'total_variants'])
        }).reset_index()
        
        # Create pivot table
        pivot_data = agg_data.pivot(index='svtype', columns='size_bucket', values=rate_col)
        pivot_counts = agg_data.pivot(index='svtype', columns='size_bucket', values='total_variants')
        
        # Reorder columns and rows
        available_sizes = [size for size in size_order if size in pivot_data.columns]
        available_svtypes = [svtype for svtype in svtype_order if svtype in pivot_data.index]
        
        if available_sizes and available_svtypes:
            pivot_data = pivot_data.reindex(index=available_svtypes, columns=available_sizes)
            pivot_counts = pivot_counts.reindex(index=available_svtypes, columns=available_sizes)
            
            # Create labels with percentages and counts
            labels = pivot_data.copy().astype(str)
            for i in range(len(labels.index)):
                for j in range(len(labels.columns)):
                    rate = pivot_data.iloc[i, j]
                    count = pivot_counts.iloc[i, j]
                    if pd.notna(rate) and pd.notna(count):
                        labels.iloc[i, j] = f"{rate:.1f}%\n({int(count)})"
                    else:
                        labels.iloc[i, j] = "N/A"
            
            # Create heatmap with square cells
            sns.heatmap(pivot_data, 
                       annot=labels, 
                       fmt='', 
                       cmap='YlOrRd', 
                       cbar_kws={'label': 'Support Rate (%)'},
                       vmin=0, 
                       vmax=100,
                       ax=ax,
                       square=True)  # Make cells square
            
            ax.set_title(f'{caller_name.title()} - {title_suffix}', fontsize=16, fontweight='bold')
            ax.set_xlabel('SV Length (bp)', fontsize=12)
            ax.set_ylabel('SV Type', fontsize=12)
            
            # Save individual data table
            pivot_data.to_csv(os.path.join(output_dir, f'support_rates_{filename_suffix}_{caller_name}.csv'))
        else:
            print(f"Warning: No valid data for heatmap {title_suffix} for caller {caller_name}")
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'heatmap_{filename_suffix}_{caller_name}.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Generated heatmap for {title_suffix} for caller {caller_name}")
    
    # Create heatmaps for both region groups
    create_heatmap_for_regions(['US', 'RM'], 'US/RM Regions', 'us_rm_regions')
    create_heatmap_for_regions(['SD', 'SR'], 'SD/SR Regions', 'sd_sr_regions')

def main():
    parser = argparse.ArgumentParser(description='Analyze multi-caller support rates across samples')
    parser.add_argument('--mapping', required=True, help='Mapping file (TSV format)')
    parser.add_argument('--input-dir', required=True, help='Input directory containing multi-caller results')
    parser.add_argument('--n-matches', type=int, default=2, help='Maximum number of caller matches to analyze (generates 1+, 2+, ..., N+ subfigures)')
    parser.add_argument('--svtypes', help='Comma-separated list of SV types to include (e.g., DEL,DUP,INS)')
    parser.add_argument('--output-dir', default='.', help='Output directory for plots and tables')
    parser.add_argument('--use-combined', action='store_true', default=True, 
                       help='Use combined.coalesced.tsv files instead of regular .coalesced.tsv files')
    parser.add_argument('--per-caller', action='store_true',
                       help='Analyze each caller individually instead of multi-caller analysis')
    parser.add_argument('--caller-names', nargs='+', default=['sniffles', 'cutesv'],
                       help='Names of individual callers to analyze when using --per-caller mode')
    
    args = parser.parse_args()
    
    # Parse SV types filter
    svtypes_filter = None
    if args.svtypes:
        svtypes_filter = [svtype.strip().upper() for svtype in args.svtypes.split(',')]
        print(f"Filtering to SV types: {', '.join(svtypes_filter)}")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    if args.per_caller:
        print(f"Running per-caller analysis for callers: {', '.join(args.caller_names)}")
        
        for caller_name in args.caller_names:
            print(f"\nProcessing caller: {caller_name}")
            print(f"Collecting data from all available samples for {caller_name}...")
            
            try:
                # Collect data for this caller
                df = collect_sample_data_per_caller(args.mapping, args.input_dir, caller_name, svtypes_filter)
                
                print(f"Loaded data for {df['sample_id'].nunique()} samples for caller {caller_name}")
                print(f"Total records: {len(df)}")
                
                # Generate visualizations for this caller
                print(f"Generating violin plots for {caller_name}...")
                caller_stats = create_violin_plots_per_caller(df, args.output_dir, caller_name, svtypes_filter)
                
                print(f"Generating heatmap tables for {caller_name}...")
                create_heatmap_tables_per_caller(df, args.output_dir, caller_name, svtypes_filter)
                
                # Save summary statistics for this caller
                if caller_stats:
                    summary_stats = {
                        'overall_mean': np.mean(caller_stats['overall']) if caller_stats['overall'] else 'N/A',
                        'overall_std': np.std(caller_stats['overall']) if caller_stats['overall'] else 'N/A',
                        'us_rm_mean': np.mean(caller_stats['us_rm']) if caller_stats['us_rm'] else 'N/A',
                        'us_rm_std': np.std(caller_stats['us_rm']) if caller_stats['us_rm'] else 'N/A',
                        'sd_sr_mean': np.mean(caller_stats['sd_sr']) if caller_stats['sd_sr'] else 'N/A',
                        'sd_sr_std': np.std(caller_stats['sd_sr']) if caller_stats['sd_sr'] else 'N/A'
                    }
                    
                    with open(os.path.join(args.output_dir, f'summary_stats_{caller_name}.txt'), 'w') as f:
                        svtype_suffix = f" ({', '.join(svtypes_filter)})" if svtypes_filter else ""
                        f.write(f"{caller_name.title()} Caller Support Rate Summary Statistics{svtype_suffix} (50% Reciprocal Overlap)\n")
                        f.write("=" * 80 + "\n\n")
                        f.write(f"Overall (all PASS+MULTIALLELIC):\n")
                        f.write(f"  Mean: {summary_stats['overall_mean']:.2f}%\n" if summary_stats['overall_mean'] != 'N/A' else "  Mean: N/A\n")
                        f.write(f"  Std:  {summary_stats['overall_std']:.2f}%\n" if summary_stats['overall_std'] != 'N/A' else "  Std: N/A\n")
                        f.write(f"\nUS/RM Regions:\n")
                        f.write(f"  Mean: {summary_stats['us_rm_mean']:.2f}%\n" if summary_stats['us_rm_mean'] != 'N/A' else "  Mean: N/A\n")
                        f.write(f"  Std:  {summary_stats['us_rm_std']:.2f}%\n" if summary_stats['us_rm_std'] != 'N/A' else "  Std: N/A\n")
                        f.write(f"\nSD/SR Regions:\n")
                        f.write(f"  Mean: {summary_stats['sd_sr_mean']:.2f}%\n" if summary_stats['sd_sr_mean'] != 'N/A' else "  Mean: N/A\n")
                        f.write(f"  Std:  {summary_stats['sd_sr_std']:.2f}%\n" if summary_stats['sd_sr_std'] != 'N/A' else "  Std: N/A\n")
                
            except Exception as e:
                print(f"Error processing caller {caller_name}: {e}")
                continue
        
        print(f"\nPer-caller analysis complete! Results saved to {args.output_dir}")
        print(f"Generated files for each caller:")
        for caller_name in args.caller_names:
            print(f"  - violin_support_{caller_name}.png")
            print(f"  - heatmap_us_rm_regions_{caller_name}.png")
            print(f"  - heatmap_sd_sr_regions_{caller_name}.png")
            print(f"  - support_rates_us_rm_regions_{caller_name}.csv")
            print(f"  - support_rates_sd_sr_regions_{caller_name}.csv")
            print(f"  - summary_stats_{caller_name}.txt")
        
    else:
        print(f"Collecting data from all available samples...")
        print(f"Analyzing 1+ through {args.n_matches}+ caller matches...")
        
        # Collect data
        df = collect_sample_data(args.mapping, args.input_dir, args.n_matches, svtypes_filter, args.use_combined)
        
        print(f"Loaded data for {df['sample_id'].nunique()} samples")
        print(f"Total records: {len(df)}")
        
        # Generate visualizations
        print("Generating violin plots...")
        all_stats = create_violin_plots(df, args.output_dir, args.n_matches, svtypes_filter)
        
        print("Generating heatmap tables...")
        create_heatmap_tables(df, args.output_dir, args.n_matches, svtypes_filter)
        
        # Save summary statistics for each n_matches level
        for n_matches in range(1, args.n_matches + 1):
            if n_matches in all_stats:
                stats = all_stats[n_matches]
                
                summary_stats = {
                    'overall_mean': np.mean(stats['overall']) if stats['overall'] else 'N/A',
                    'overall_std': np.std(stats['overall']) if stats['overall'] else 'N/A',
                    'us_rm_mean': np.mean(stats['us_rm']) if stats['us_rm'] else 'N/A',
                    'us_rm_std': np.std(stats['us_rm']) if stats['us_rm'] else 'N/A',
                    'sd_sr_mean': np.mean(stats['sd_sr']) if stats['sd_sr'] else 'N/A',
                    'sd_sr_std': np.std(stats['sd_sr']) if stats['sd_sr'] else 'N/A'
                }
                
                with open(os.path.join(args.output_dir, f'summary_stats_{n_matches}_callers.txt'), 'w') as f:
                    svtype_suffix = f" ({', '.join(svtypes_filter)})" if svtypes_filter else ""
                    f.write(f"Multi-Caller Support Rate Summary Statistics ({n_matches}+ Callers{svtype_suffix}, 50% Reciprocal Overlap)\n")
                    f.write("=" * 80 + "\n\n")
                    f.write(f"Overall (all PASS+MULTIALLELIC):\n")
                    f.write(f"  Mean: {summary_stats['overall_mean']:.2f}%\n" if summary_stats['overall_mean'] != 'N/A' else "  Mean: N/A\n")
                    f.write(f"  Std:  {summary_stats['overall_std']:.2f}%\n" if summary_stats['overall_std'] != 'N/A' else "  Std: N/A\n")
                    f.write(f"\nUS/RM Regions:\n")
                    f.write(f"  Mean: {summary_stats['us_rm_mean']:.2f}%\n" if summary_stats['us_rm_mean'] != 'N/A' else "  Mean: N/A\n")
                    f.write(f"  Std:  {summary_stats['us_rm_std']:.2f}%\n" if summary_stats['us_rm_std'] != 'N/A' else "  Std: N/A\n")
                    f.write(f"\nSD/SR Regions:\n")
                    f.write(f"  Mean: {summary_stats['sd_sr_mean']:.2f}%\n" if summary_stats['sd_sr_mean'] != 'N/A' else "  Mean: N/A\n")
                    f.write(f"  Std:  {summary_stats['sd_sr_std']:.2f}%\n" if summary_stats['sd_sr_std'] != 'N/A' else "  Std: N/A\n")
        
        svtype_suffix = f"_{'-'.join(svtypes_filter)}" if svtypes_filter else ""
        print(f"\nAnalysis complete! Results saved to {args.output_dir}")
        print(f"Generated files:")
        print(f"  - violin_combined_support_{args.n_matches}_callers.png")
        print(f"  - heatmap_us_rm_regions_{args.n_matches}_callers.png")
        print(f"  - heatmap_sd_sr_regions_{args.n_matches}_callers.png")
        for n in range(1, args.n_matches + 1):
            print(f"  - support_rates_us_rm_regions_{n}_callers.csv")
            print(f"  - support_rates_sd_sr_regions_{n}_callers.csv")
            print(f"  - summary_stats_{n}_callers.txt")

if __name__ == '__main__':
    main() 