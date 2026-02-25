#!/bin/python

"""
Plot GD CNV Detection Results

This script generates visualization plots for the output of gd_cnv_pyro.py,
showing depth profiles at GD loci with annotations for:
- Segmental duplication regions
- Gene transcripts from GTF file
- Breakpoint intervals
- Copy number calls

Usage:
    python plot_gd_cnv_output.py \
        --calls gd_cnv_calls.tsv.gz \
        --depth-data binned_counts.tsv.gz \
        --gd-table gd_table.tsv \
        --gtf genes.gtf.gz \
        --segdup-bed segdups.bed.gz \
        --output-dir plots/
"""

import argparse
import gzip
import os
import re
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

# Import GDLocus and GDTable from gd_cnv_pyro.py to avoid duplication
from gd_cnv_pyro import GDLocus, GDTable


# =============================================================================
# Data Loading Classes
# =============================================================================


class GTFParser:
    """Parser for GTF gene annotation files."""

    def __init__(self, filepath: str):
        """
        Load gene annotations from GTF file.

        Args:
            filepath: Path to GTF file (can be gzipped)
        """
        self.genes = []
        self.transcripts = []
        self._parse_gtf(filepath)
        print(f"Loaded {len(self.genes)} genes and {len(self.transcripts)} transcripts")

    def _parse_gtf(self, filepath: str):
        """Parse GTF file and extract gene/transcript information."""
        opener = gzip.open if filepath.endswith(".gz") else open
        mode = "rt" if filepath.endswith(".gz") else "r"

        with opener(filepath, mode) as f:
            for line in f:
                if line.startswith("#"):
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 9:
                    continue

                chrom, source, feature, start, end, score, strand, frame, attributes = fields

                if feature not in ["gene", "transcript"]:
                    continue

                # Parse attributes
                attr_dict = self._parse_attributes(attributes)

                record = {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "gene_id": attr_dict.get("gene_id", ""),
                    "gene_name": attr_dict.get("gene_name", ""),
                    "gene_type": attr_dict.get("gene_type", ""),
                }

                if feature == "gene":
                    self.genes.append(record)
                elif feature == "transcript":
                    record["transcript_id"] = attr_dict.get("transcript_id", "")
                    record["transcript_type"] = attr_dict.get("transcript_type", "")
                    self.transcripts.append(record)

    def _parse_attributes(self, attr_string: str) -> dict:
        """Parse GTF attribute string into dictionary."""
        attrs = {}
        # Match key "value" or key value patterns
        pattern = r'(\w+)\s+"([^"]+)"'
        for match in re.finditer(pattern, attr_string):
            key, value = match.groups()
            attrs[key] = value
        return attrs

    def get_genes_in_region(self, chrom: str, start: int, end: int,
                            gene_types: Optional[List[str]] = None) -> List[dict]:
        """
        Get genes overlapping a genomic region.

        Args:
            chrom: Chromosome name
            start: Region start
            end: Region end
            gene_types: Optional list of gene types to filter (e.g., ["protein_coding"])

        Returns:
            List of gene records
        """
        result = []
        for gene in self.genes:
            if gene["chrom"] != chrom:
                continue
            if gene["end"] < start or gene["start"] > end:
                continue
            if gene_types and gene["gene_type"] not in gene_types:
                continue
            result.append(gene)
        return result

    def get_transcripts_in_region(self, chrom: str, start: int, end: int,
                                   transcript_types: Optional[List[str]] = None) -> List[dict]:
        """Get transcripts overlapping a genomic region."""
        result = []
        for tx in self.transcripts:
            if tx["chrom"] != chrom:
                continue
            if tx["end"] < start or tx["start"] > end:
                continue
            if transcript_types and tx["transcript_type"] not in transcript_types:
                continue
            result.append(tx)
        return result


class SegDupAnnotation:
    """Handler for segmental duplication annotations."""

    def __init__(self, filepath: str):
        """Load segmental duplication regions from BED file."""
        self.df = pd.read_csv(
            filepath,
            sep="\t",
            header=None,
            names=["chr", "start", "end", "name", "score", "strand"],
            compression="gzip" if filepath.endswith(".gz") else None,
        )
        self._build_index()
        print(f"Loaded {len(self.df)} segmental duplication regions")

    def _build_index(self):
        """Build index by chromosome."""
        self.by_chrom = {}
        for chrom, group in self.df.groupby("chr"):
            self.by_chrom[chrom] = group[["start", "end"]].values

    def get_regions_in_range(self, chrom: str, start: int, end: int) -> List[Tuple[int, int]]:
        """Get segdup regions overlapping a genomic range."""
        if chrom not in self.by_chrom:
            return []

        regions = self.by_chrom[chrom]
        overlaps = (regions[:, 0] < end) & (regions[:, 1] > start)
        return [(int(r[0]), int(r[1])) for r in regions[overlaps]]


# =============================================================================
# CNV Calling Functions
# =============================================================================

def get_locus_interval_bins(
    bin_mappings_df: pd.DataFrame,
    cluster: str,
) -> Dict[str, List[int]]:
    """
    Get bin indices for each interval in a specific locus.
    
    Args:
        bin_mappings_df: DataFrame with bin-to-interval mappings
        cluster: Cluster name to filter by
    
    Returns:
        Dict mapping interval name to list of bin indices (array_idx from posteriors)
    """
    locus_bins = bin_mappings_df[bin_mappings_df["cluster"] == cluster]
    interval_bins = {}
    for interval_name, group in locus_bins.groupby("interval"):
        # Use array_idx directly - this is the index in the posterior arrays
        interval_bins[interval_name] = group.index.tolist()
    return interval_bins


def compute_interval_cn_stats(
    cn_posteriors_df: pd.DataFrame,
    interval_bins: Dict[str, List[int]],
    sample_id: str,
) -> Dict[str, dict]:
    """
    Compute copy number statistics for each interval for a specific sample.
    
    Args:
        cn_posteriors_df: DataFrame with CN posterior probabilities (indexed by bin position)
        interval_bins: Dict mapping interval name to bin row indices in posteriors DataFrame
        sample_id: Sample identifier
    
    Returns:
        Dict mapping interval name to CN statistics
    """
    result = {}
    
    # Filter posteriors to this sample
    sample_posteriors = cn_posteriors_df[cn_posteriors_df["sample"] == sample_id]
    
    # Get probability columns (prob_cn_0, prob_cn_1, etc.)
    prob_cols = [c for c in sample_posteriors.columns if c.startswith("prob_cn_")]
    n_states = len(prob_cols)
    
    for interval_name, bin_indices in interval_bins.items():
        # Get bins for this interval
        interval_data = sample_posteriors.iloc[bin_indices]
        
        if len(interval_data) == 0:
            result[interval_name] = {
                "n_bins": 0,
                "cn_probs": np.zeros(n_states),
            }
            continue
        
        # Average posterior probabilities across bins
        mean_probs = interval_data[prob_cols].mean(axis=0).values
        
        result[interval_name] = {
            "n_bins": len(interval_data),
            "cn_probs": mean_probs,
        }
    
    return result


def call_gd_cnv(
    locus: GDLocus,
    interval_stats: Dict[str, dict],
    log_prob_threshold: float = -0.5,
    flanking_log_prob_threshold: float = -1.0,
) -> List[dict]:
    """
    Call GD CNVs based on interval copy number statistics.
    
    For each GD entry in the locus (each DEL/DUP definition), check if the
    copy number in the corresponding region supports a CNV call.
    
    Args:
        locus: GDLocus object
        interval_stats: Dict mapping interval name to CN statistics
        log_prob_threshold: Minimum log probability score to call a CNV
        flanking_log_prob_threshold: Minimum log probability score in flanking regions to classify as spanning
    
    Returns:
        List of CNV call dictionaries
    """
    calls = []
    
    for entry in locus.gd_entries:
        gd_id = entry["GD_ID"]
        svtype = entry["svtype"]
        gd_start = entry["start_GRCh38"]
        gd_end = entry["end_GRCh38"]
        bp1 = entry["BP1"]
        bp2 = entry["BP2"]
        
        # Determine which intervals this entry spans using BP1 and BP2
        all_intervals = locus.get_intervals()  # List of (start, end, name) tuples
        
        # Build breakpoint ordering from interval names
        bp_order = []
        for _, _, interval_name in all_intervals:
            parts = interval_name.split("-")
            if len(parts) == 2:
                if parts[0] not in bp_order:
                    bp_order.append(parts[0])
                if parts[1] not in bp_order:
                    bp_order.append(parts[1])
        
        # Find positions of BP1 and BP2 in the ordering
        try:
            pos1 = bp_order.index(bp1)
            pos2 = bp_order.index(bp2)
        except ValueError:
            continue
        
        # Ensure pos1 < pos2
        if pos1 > pos2:
            pos1, pos2 = pos2, pos1
        
        # Collect all intervals between these breakpoints
        covered_intervals = []
        for _, _, interval_name in all_intervals:
            parts = interval_name.split("-")
            if len(parts) == 2:
                try:
                    start_pos = bp_order.index(parts[0])
                    end_pos = bp_order.index(parts[1])
                    if start_pos >= pos1 and end_pos <= pos2:
                        covered_intervals.append(interval_name)
                except ValueError:
                    continue
        
        if len(covered_intervals) == 0:
            continue
        
        # Compute weighted average log probability score
        log_prob_score = 0.0
        total_weight = 0.0
        
        if svtype == "DEL":
            for interval in covered_intervals:
                if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                    probs = interval_stats[interval]["cn_probs"]
                    weight = interval_stats[interval]["n_bins"]
                    p_del = max(probs[0] + probs[1], 1e-5)  # P(CN=0) + P(CN=1)
                    if p_del > 0:
                        log_prob_score += weight * np.log(p_del)
                        total_weight += weight
        
        elif svtype == "DUP":
            for interval in covered_intervals:
                if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                    probs = interval_stats[interval]["cn_probs"]
                    weight = interval_stats[interval]["n_bins"]
                    p_dup = max(probs[3:].sum(), 1e-5)  # P(CN >= 3)
                    if p_dup > 0:
                        log_prob_score += weight * np.log(p_dup)
                        total_weight += weight
        
        # Normalize by total weight
        if total_weight > 0:
            log_prob_score = log_prob_score / total_weight
        else:
            log_prob_score = np.nan
        
        # Determine if this is a CNV
        is_cnv = log_prob_score > log_prob_threshold

        # Check flanking regions for evidence of a large spanning CNV.
        # A true GD only affects the region between its breakpoints; a spanning
        # variant would also show copy number change in the flanking regions.
        flanking_log_prob_score = np.nan
        is_spanning = False
        if is_cnv:
            flank_score_sum = 0.0
            flank_weight_total = 0.0
            for _, _, flank_name in locus.get_flanking_regions():
                if flank_name not in interval_stats:
                    continue
                flank_stat = interval_stats[flank_name]
                if flank_stat["n_bins"] == 0:
                    continue
                flank_probs = flank_stat["cn_probs"]
                flank_weight = flank_stat["n_bins"]
                if svtype == "DEL":
                    p_flank = max(flank_probs[0] + flank_probs[1], 1e-5)
                elif svtype == "DUP":
                    p_flank = max(flank_probs[3:].sum(), 1e-5)
                else:
                    continue
                flank_score_sum += flank_weight * np.log(p_flank)
                flank_weight_total += flank_weight

            if flank_weight_total > 0:
                flanking_log_prob_score = flank_score_sum / flank_weight_total
                # If flanking regions also show a CN change, this is a spanning variant
                if flanking_log_prob_score > flanking_log_prob_threshold:
                    is_spanning = True
                    is_cnv = False

        # Count total bins
        n_bins = sum(interval_stats[iv]["n_bins"] for iv in covered_intervals if iv in interval_stats)

        calls.append({
            "GD_ID": gd_id,
            "cluster": locus.cluster,
            "chrom": locus.chrom,
            "start": gd_start,
            "end": gd_end,
            "svtype": svtype,
            "is_nahr": locus.is_nahr,
            "is_terminal": locus.is_terminal,
            "log_prob_score": log_prob_score,
            "flanking_log_prob_score": flanking_log_prob_score,
            "is_carrier": is_cnv,
            "is_spanning": is_spanning,
            "intervals": covered_intervals,
            "n_bins": n_bins,
        })
    
    return calls


def determine_best_breakpoints(
    locus: GDLocus,
    interval_stats: Dict[str, dict],
    calls: List[dict],
) -> Dict[str, Optional[str]]:
    """
    Determine the most likely breakpoint pair for a GD CNV, separately for DEL and DUP.
    
    For loci with multiple possible breakpoint configurations (e.g., BP1-2, BP1-3),
    determine which best fits the observed data for each svtype.
    
    Args:
        locus: GDLocus object
        interval_stats: Dict mapping interval name to CN statistics
        calls: List of CNV call dictionaries
    
    Returns:
        Dict mapping svtype ("DEL", "DUP") to best matching GD_ID or None
    """
    best_by_svtype = {}
    
    # Get all interval names for this locus
    all_intervals = set(name for _, _, name in locus.get_intervals())
    
    for svtype in ["DEL", "DUP"]:
        # Filter to carrier calls of this svtype
        carrier_calls = [c for c in calls if c["is_carrier"] and c["svtype"] == svtype]
        
        if len(carrier_calls) == 0:
            best_by_svtype[svtype] = None
            continue
        
        if len(carrier_calls) == 1:
            best_by_svtype[svtype] = carrier_calls[0]["GD_ID"]
            continue
        
        # Multiple carrier calls - score each
        best_score = -np.inf
        best_gd_id = None
        
        for call in carrier_calls:
            covered_intervals = set(call["intervals"])
            uncovered_intervals = all_intervals - covered_intervals
            
            score = 0.0
            total_weight = 0.0
            
            if svtype == "DEL":
                # Score affected intervals
                for interval in covered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = 1.  # interval_stats[interval]["n_bins"]
                        p_del = probs[0] + probs[1]
                        if p_del > 0:
                            score += weight * np.log(p_del)
                            total_weight += weight
                
                # Score unaffected intervals
                for interval in uncovered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = 1.  # iinterval_stats[interval]["n_bins"]
                        p_normal = probs[2:].sum()
                        if p_normal > 0:
                            score += weight * np.log(p_normal)
                            total_weight += weight
            
            else:  # DUP
                # Score affected intervals
                for interval in covered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = 1.  # iinterval_stats[interval]["n_bins"]
                        p_dup = probs[3:].sum()
                        if p_dup > 0:
                            score += weight * np.log(p_dup)
                            total_weight += weight
                
                # Score unaffected intervals
                for interval in uncovered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = 1.  # interval_stats[interval]["n_bins"]
                        p_normal = probs[:3].sum()
                        if p_normal > 0:
                            score += weight * np.log(p_normal)
                            total_weight += weight
            
            # Normalize
            if total_weight > 0:
                score = score / total_weight
            
            if score > best_score:
                best_score = score
                best_gd_id = call["GD_ID"]
        
        best_by_svtype[svtype] = best_gd_id
    
    return best_by_svtype


def call_cnvs_from_posteriors(
    cn_posteriors_df: pd.DataFrame,
    bin_mappings_df: pd.DataFrame,
    gd_table: GDTable,
    log_prob_threshold: float = -0.5,
    flanking_log_prob_threshold: float = -1.0,
) -> pd.DataFrame:
    """
    Call CNVs from posterior probabilities.
    
    Args:
        cn_posteriors_df: DataFrame with CN posterior probabilities per bin/sample
        bin_mappings_df: DataFrame with bin-to-interval mappings
        gd_table: GDTable with locus definitions
        log_prob_threshold: Minimum log probability score to call a CNV
        flanking_log_prob_threshold: Minimum log probability score in flanking regions to classify as spanning
    
    Returns:
        DataFrame with CNV calls for all samples and loci
    """
    print("\n" + "=" * 80)
    print("CALLING CNVs FROM POSTERIORS")
    print("=" * 80)
    
    all_results = []
    sample_ids = cn_posteriors_df["sample"].unique()
    
    # Build mapping of bin position to row index for faster lookup
    cn_posteriors_df = cn_posteriors_df.reset_index(drop=True)
    bin_mappings_df = bin_mappings_df.reset_index(drop=True)
    
    # Pre-group posteriors by sample for faster lookups
    print("  Organizing data for fast access...")
    posteriors_by_sample = {}
    for sample_id in sample_ids:
        sample_data = cn_posteriors_df[cn_posteriors_df["sample"] == sample_id]
        posteriors_by_sample[sample_id] = sample_data
    
    for cluster, locus in gd_table.loci.items():
        print(f"\nCalling CNVs for locus: {cluster}")
        
        # Get interval bins for this locus
        interval_bins = get_locus_interval_bins(bin_mappings_df, cluster)
        
        for interval_name, bins in interval_bins.items():
            print(f"  {interval_name}: {len(bins)} bins")
        
        # Process each sample
        for sample_id in sample_ids:
            # Use pre-filtered sample data
            sample_posteriors = posteriors_by_sample[sample_id]
            
            # Get probability columns once
            prob_cols = [c for c in sample_posteriors.columns if c.startswith("prob_cn_")]
            
            # Compute interval statistics efficiently
            interval_stats = {}
            for interval_name, bin_indices in interval_bins.items():
                if len(bin_indices) == 0:
                    interval_stats[interval_name] = {
                        "n_bins": 0,
                        "cn_probs": np.zeros(len(prob_cols)),
                    }
                else:
                    interval_data = sample_posteriors.iloc[bin_indices]
                    mean_probs = interval_data[prob_cols].mean(axis=0).values
                    interval_stats[interval_name] = {
                        "n_bins": len(interval_data),
                        "cn_probs": mean_probs,
                    }
            
            # Call GD CNVs
            calls = call_gd_cnv(locus, interval_stats, log_prob_threshold, flanking_log_prob_threshold)
            
            # Determine best breakpoints
            best_by_svtype = determine_best_breakpoints(locus, interval_stats, calls)
            
            # Record results
            for call in calls:
                svtype = call["svtype"]
                best_gd_for_svtype = best_by_svtype.get(svtype)
                
                # Compute mean depth over covered intervals efficiently
                covered_bin_indices = []
                for interval_name in call["intervals"]:
                    if interval_name in interval_bins:
                        covered_bin_indices.extend(interval_bins[interval_name])
                
                if len(covered_bin_indices) > 0:
                    # Use iloc for fast indexing on pre-filtered data
                    mean_depth = sample_posteriors.iloc[covered_bin_indices]["depth"].mean()
                else:
                    mean_depth = np.nan
                
                result = {
                    "sample": sample_id,
                    "cluster": cluster,
                    "GD_ID": call["GD_ID"],
                    "chrom": call["chrom"],
                    "start": call["start"],
                    "end": call["end"],
                    "svtype": svtype,
                    "is_nahr": call["is_nahr"],
                    "is_terminal": call["is_terminal"],
                    "n_bins": call["n_bins"],
                    "mean_depth": mean_depth,
                    "is_carrier": call["is_carrier"],
                    "is_best_match": call["GD_ID"] == best_gd_for_svtype if best_gd_for_svtype else False,
                    "log_prob_score": call["log_prob_score"],
                }
                all_results.append(result)
    
    return pd.DataFrame(all_results)


# =============================================================================
# Plotting Functions
# =============================================================================

def get_sample_columns(df: pd.DataFrame) -> list:
    """Extract sample column names from DataFrame."""
    metadata_cols = ["Cluster", "Chr", "Start", "End", "source_file", "Bin"]
    return [col for col in df.columns if col not in metadata_cols]


def draw_annotations_panel(
    ax,
    locus: GDLocus,
    region_start: int,
    region_end: int,
    chrom: str,
    title: str,
    gtf: Optional[GTFParser] = None,
    segdup: Optional[SegDupAnnotation] = None,
    show_gd_entries: bool = True,
):
    """
    Draw annotations panel with genes, segdups, and breakpoints.
    
    Args:
        ax: Matplotlib axis to draw on
        locus: GDLocus object
        region_start: Region start position
        region_end: Region end position
        chrom: Chromosome name
        title: Panel title
        gtf: Optional GTFParser for gene annotations
        segdup: Optional SegDupAnnotation
        show_gd_entries: Whether to show GD entry regions and breakpoint intervals
    """
    ax.set_xlim(region_start, region_end)
    ax.set_ylim(0, 0.25)
    ax.set_ylabel("Annotations")
    ax.set_title(title)

    # Shade flanking regions using the actual data-derived extents.
    # region_start/region_end come from min/max of the bins present in depth_df,
    # so these spans exactly match what is plotted — no geometric approximation.
    flank_defs = []
    if region_start < locus.start:
        flank_defs.append((region_start, locus.start, "left flank"))
    if locus.end < region_end:
        flank_defs.append((locus.end, region_end, "right flank"))
    for flank_start, flank_end, flank_name in flank_defs:
        ax.axvspan(flank_start, flank_end, alpha=0.08, color="black", zorder=0)
        ax.text((flank_start + flank_end) / 2, 0.225,
                flank_name, ha="center", va="center",
                fontsize=7, color="gray", style="italic")

    # Draw breakpoint intervals (only if show_gd_entries is True)
    if show_gd_entries:
        intervals = locus.get_intervals()
        colors = plt.cm.Set3(np.linspace(0, 1, len(intervals)))
        for i, (start, end, name) in enumerate(intervals):
            ax.axvspan(start, end, alpha=0.15, color=colors[i], label=name, zorder=0)
            ax.text((start + end) / 2, 0.225, name, ha="center", va="center", fontsize=7)

    # Draw breakpoint ranges as shaded regions (back layer)
    for i, (bp_start, bp_end) in enumerate(locus.breakpoints):
        ax.axvspan(bp_start, bp_end, ymin=0.04, ymax=0.15, alpha=0.3, color="red", zorder=1)
        bp_name = locus.breakpoint_names[i] if i < len(locus.breakpoint_names) else str(i+1)
        ax.text((bp_start + bp_end) / 2, 0.02, f"BP{bp_name}", ha="center", va="center", 
                fontsize=7, fontweight="bold", color="darkred", zorder=10)

    # Draw segdup regions (middle layer)
    if segdup:
        sd_regions = segdup.get_regions_in_range(chrom, region_start, region_end)
        for sd_start, sd_end in sd_regions:
            ax.axvspan(max(sd_start, region_start), min(sd_end, region_end),
                      ymin=0.2, ymax=0.5, alpha=0.2, color="orange", zorder=2)

    # Draw genes (top layer)
    if gtf:
        genes = gtf.get_genes_in_region(chrom, region_start, region_end,
                                         gene_types=["protein_coding"])
        y_pos = 0.2
        min_gene_size = 20000  # Minimum 20kb to show label
        for gene in genes[:10]:  # Limit to 10 genes
            gene_start = max(gene["start"], region_start)
            gene_end = min(gene["end"], region_end)
            ax.hlines(y_pos, gene_start, gene_end, colors="blue", linewidth=4, zorder=3)
            # Only label genes larger than minimum size
            if (gene_end - gene_start) >= min_gene_size:
                ax.text((gene_start + gene_end) / 2, y_pos - 0.05,
                       gene["gene_name"], ha="center", va="bottom", fontsize=7, zorder=3)

    ax.set_yticks([])


def plot_locus_overview(
    locus: GDLocus,
    calls_df: pd.DataFrame,
    depth_df: pd.DataFrame,
    gtf: Optional[GTFParser],
    segdup: Optional[SegDupAnnotation],
    output_dir: str,
    padding: int = 50000,
):
    """
    Create overview plot for a GD locus showing all samples.

    Args:
        locus: GDLocus object
        calls_df: DataFrame with CNV calls for this locus
        depth_df: DataFrame with depth data (from cn_posteriors)
        gtf: Optional GTFParser for gene annotations
        segdup: Optional SegDupAnnotation for segdup regions
        output_dir: Directory to save plots
        padding: Padding around locus boundaries
    """
    # Skip loci with no breakpoints
    if not locus.breakpoints:
        print(f"  Warning: No breakpoints defined for locus {locus.cluster}, skipping")
        return
    
    chrom = locus.chrom
    # Filter strictly to this locus's bins (by cluster name). Using position alone
    # would pull in differently-rebinned bins from neighbouring/overlapping loci.
    mask = (
        (depth_df["Cluster"] == locus.cluster) &
        (depth_df["Chr"] == chrom)
    )
    region_df = depth_df[mask].copy().sort_values("Start")

    if len(region_df) == 0:
        print(f"  Warning: No depth data found for locus {locus.cluster} ({chrom}:{locus.start:,}-{locus.end:,}), skipping plot")
        return

    # Bounds derived from actual data so flanks are always visible
    region_start = int(region_df["Start"].min())
    region_end = int(region_df["End"].max())

    sample_cols = get_sample_columns(region_df)
    n_samples = len(sample_cols)

    # Get carriers for this locus
    carriers = calls_df[
        (calls_df["cluster"] == locus.cluster) &
        (calls_df["is_carrier"])
    ]["sample"].unique()

    # Separate carriers and non-carriers
    carrier_cols = [s for s in sample_cols if s in carriers]
    non_carrier_cols = [s for s in sample_cols if s not in carriers]
    n_carriers = len(carrier_cols)
    n_non_carriers = len(non_carrier_cols)

    # Create figure with 4 panels: annotations, carriers heatmap, non-carriers heatmap, mean depth
    # Height ratios: scale heatmap heights by number of samples (capped at 6 inches each)
    carrier_height = min(6, max(1, n_carriers * 0.15)) if n_carriers > 0 else 0.5
    non_carrier_height = min(6, max(1, n_non_carriers * 0.15)) if n_non_carriers > 0 else 0.5
    fig_height = 4 + carrier_height + non_carrier_height
    
    fig, axes = plt.subplots(4, 1, figsize=(8, fig_height), 
                              gridspec_kw={"height_ratios": [1, carrier_height, non_carrier_height, 1]})

    # X-axis: genomic position
    bin_mids = (region_df["Start"].values + region_df["End"].values) / 2

    # Panel 1: Gene annotations and segdup regions
    ax = axes[0]
    title = f"{locus.cluster} ({chrom}:{locus.start:,}-{locus.end:,})"
    draw_annotations_panel(
        ax, locus, region_start, region_end, chrom, title,
        gtf=gtf, segdup=segdup, show_gd_entries=True
    )

    # Panel 2: Carriers heatmap
    ax = axes[1]
    if n_carriers > 0:
        # Create a mapping from genomic position to bin index
        bin_starts = region_df["Start"].values
        bin_ends = region_df["End"].values
        
        # Create dense matrix with NaN for missing bins
        # Use a fixed bin size for the visualization
        viz_bin_size = int(np.median(bin_ends - bin_starts))
        n_viz_bins = int((region_end - region_start) / viz_bin_size) + 1
        viz_matrix = np.full((n_carriers, n_viz_bins), np.nan)
        
        # Fill in actual data at correct positions
        carrier_data = region_df[carrier_cols].values  # bins x carriers
        for i, (start, end) in enumerate(zip(bin_starts, bin_ends)):
            viz_start_idx = int((start - region_start) / viz_bin_size)
            viz_end_idx = int((end - region_start) / viz_bin_size) + 1
            # Assign data to all visualization bins covered by this real bin
            for viz_idx in range(viz_start_idx, min(viz_end_idx, n_viz_bins)):
                viz_matrix[:, viz_idx] = carrier_data[i, :]
        
        # Create colormap with white for NaN values
        cmap = plt.cm.RdBu_r.copy()
        cmap.set_bad(color='white')
        
        im1 = ax.imshow(viz_matrix, aspect="auto", cmap=cmap,
                       vmin=0, vmax=4, interpolation="nearest",
                       extent=[region_start, region_end, n_carriers, 0])
        
        # Draw breakpoint ranges
        # for i, (bp_start, bp_end) in enumerate(locus.breakpoints):
        #     ax.axvline(bp_start, color="red", linestyle="-", alpha=1, zorder=5)
        #     ax.axvline(bp_end, color="red", linestyle="-", alpha=1, zorder=5)
        
        ax.set_ylabel(f"Carriers (n={n_carriers})")
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_xlim(region_start, region_end)
        ax.set_ylim(0, n_carriers)
    else:
        ax.axis('off')
        ax.text(0.5, 0.5, 'No carriers', ha='center', va='center', transform=ax.transAxes)

    # Panel 3: Non-carriers heatmap
    ax = axes[2]
    if n_non_carriers > 0:
        # Create a mapping from genomic position to bin index
        bin_starts = region_df["Start"].values
        bin_ends = region_df["End"].values
        
        # Create dense matrix with NaN for missing bins
        viz_bin_size = int(np.median(bin_ends - bin_starts))
        n_viz_bins = int((region_end - region_start) / viz_bin_size) + 1
        viz_matrix = np.full((n_non_carriers, n_viz_bins), np.nan)
        
        # Fill in actual data at correct positions
        non_carrier_data = region_df[non_carrier_cols].values  # bins x non_carriers
        for i, (start, end) in enumerate(zip(bin_starts, bin_ends)):
            viz_start_idx = int((start - region_start) / viz_bin_size)
            viz_end_idx = int((end - region_start) / viz_bin_size) + 1
            # Assign data to all visualization bins covered by this real bin
            for viz_idx in range(viz_start_idx, min(viz_end_idx, n_viz_bins)):
                viz_matrix[:, viz_idx] = non_carrier_data[i, :]
        
        # Create colormap with white for NaN values
        cmap = plt.cm.RdBu_r.copy()
        cmap.set_bad(color='white')
        
        im2 = ax.imshow(viz_matrix, aspect="auto", cmap=cmap,
                       vmin=0, vmax=4, interpolation="nearest",
                       extent=[region_start, region_end, n_non_carriers, 0])
        
        # Draw breakpoint ranges
        # for i, (bp_start, bp_end) in enumerate(locus.breakpoints):
        #     ax.axvline(bp_start, color="red", linestyle="-", alpha=1, zorder=5)
        #     ax.axvline(bp_end, color="red", linestyle="-", alpha=1, zorder=5)
        
        ax.set_ylabel(f"Non-carriers (n={n_non_carriers})")
        ax.set_yticks([])
        ax.set_xlim(region_start, region_end)
        ax.set_ylim(0, n_non_carriers)
        ax.set_xlabel(f"Position on {chrom}")
        
        # Add colorbar below non-carriers heatmap
        cbar = plt.colorbar(im2, ax=ax, orientation='horizontal', pad=0.08, aspect=40)
        cbar.set_label("Normalized read depth")
    else:
        ax.axis('off')
        ax.text(0.5, 0.5, 'No non-carriers', ha='center', va='center', transform=ax.transAxes)

    # Panel 4: Mean depth profile with carrier/non-carrier separation
    ax = axes[3]
    
    # Detect gaps in bins and insert NaN to prevent interpolation
    bin_starts = region_df["Start"].values
    bin_ends = region_df["End"].values
    
    # Create arrays with NaN inserted at gaps
    plot_positions = []
    
    # Determine typical bin spacing
    bin_spacing = np.median(bin_starts[1:] - bin_ends[:-1])
    gap_threshold = 3 * bin_spacing  # Consider it a gap if spacing is 3x typical
    
    # Build position array with NaN markers at gaps
    for i in range(len(bin_starts)):
        plot_positions.append(bin_mids[i])
        if i < len(bin_starts) - 1:
            # Check if there's a gap before next bin
            gap_size = bin_starts[i + 1] - bin_ends[i]
            if gap_size > gap_threshold:
                # Insert NaN marker
                plot_positions.append(np.nan)
    
    plot_positions = np.array(plot_positions)

    # Calculate mean depth for carriers and non-carriers
    if n_carriers > 0 and n_carriers < n_samples:
        # Plot non-carriers with gaps
        non_carrier_depths = region_df[[s for s in sample_cols if s not in carriers]].mean(axis=1).values
        plot_depths = []
        j = 0
        for pos in plot_positions:
            if np.isnan(pos):
                plot_depths.append(np.nan)
            else:
                plot_depths.append(non_carrier_depths[j])
                j += 1
        
        ax.plot(plot_positions, plot_depths, "b-", linewidth=1, alpha=0.7,
               label=f"Non-carriers (n={n_samples - n_carriers})")
        
        # Group carriers by SV type and breakpoint combination
        # Only include the best match call for each sample
        carrier_calls = calls_df[
            (calls_df["cluster"] == locus.cluster) &
            (calls_df["is_carrier"]) &
            (calls_df["is_best_match"]) &  # Only the best match for each sample
            (calls_df["sample"].isin(sample_cols))  # Only samples with depth data in this region
        ].copy()
        
        # Get BP1 and BP2 labels from GD table entries by matching coordinates
        def get_bp_labels(row):
            # Find matching GD entry by coordinates and svtype
            for entry in locus.gd_entries:
                if (entry["start_GRCh38"] == row["start"] and 
                    entry["end_GRCh38"] == row["end"] and
                    entry["svtype"] == row["svtype"]):
                    return f"{entry['BP1']}-{entry['BP2']}"
            # Fallback: try to match just coordinates
            for entry in locus.gd_entries:
                if entry["start_GRCh38"] == row["start"] and entry["end_GRCh38"] == row["end"]:
                    return f"{entry['BP1']}-{entry['BP2']}"
            return "unknown"
        
        carrier_calls['bp_interval'] = carrier_calls.apply(get_bp_labels, axis=1)
        
        # Define colors for different groups
        del_colors = ['#FF6B6B', '#FF4757', '#EE5A6F', '#C23B4E']  # Reds for DEL
        dup_colors = ['#6B9BD1', '#4169E1', '#5080D0', '#3A5BB8']  # Blues for DUP
        
        # Group by svtype and interval
        grouped = carrier_calls.groupby(['svtype', 'bp_interval'])
        color_idx = 0
        plotted_any = False
        
        for (svtype, bp_interval), group in grouped:
            group_samples = group['sample'].unique()
            # Filter to samples that exist in this region's data
            valid_samples = [s for s in group_samples if s in sample_cols]
            
            if len(valid_samples) > 0:
                group_depths_raw = region_df[valid_samples].mean(axis=1).values
                
                # Insert NaN at gaps
                plot_depths = []
                j = 0
                for pos in plot_positions:
                    if np.isnan(pos):
                        plot_depths.append(np.nan)
                    else:
                        plot_depths.append(group_depths_raw[j])
                        j += 1
                
                # Only plot if there's actual data
                if len(plot_depths) > 0 and not all(np.isnan(plot_depths)):
                    # Choose color based on svtype
                    if svtype == "DEL":
                        color = del_colors[color_idx % len(del_colors)]
                    else:  # DUP
                        color = dup_colors[color_idx % len(dup_colors)]
                    ax.plot(plot_positions, plot_depths, "-", linewidth=1.5, alpha=0.8,
                           color=color, label=f"{svtype} {bp_interval} (n={len(valid_samples)})")
                    color_idx += 1
                    plotted_any = True
        
        # Only show legend if we actually plotted something
        if plotted_any:
            ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.35), 
                     ncol=6, fontsize=8, frameon=True, fancybox=True)
    else:
        mean_depth_raw = region_df[sample_cols].mean(axis=1).values
        
        # Insert NaN at gaps
        plot_depths = []
        j = 0
        for pos in plot_positions:
            if np.isnan(pos):
                plot_depths.append(np.nan)
            else:
                plot_depths.append(mean_depth_raw[j])
                j += 1
        
        ax.plot(plot_positions, plot_depths, "b-", linewidth=1.5)

    # Draw breakpoint ranges
    for i, (bp_start, bp_end) in enumerate(locus.breakpoints):
        ax.axvline(bp_start, color="red", linestyle="-", alpha=1, zorder=5)
        ax.axvline(bp_end, color="red", linestyle="-", alpha=1, zorder=5)
        ax.axvspan(bp_start, bp_end, alpha=0.1, color="red", zorder=0)

    ax.axhline(2.0, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlim(region_start, region_end)
    ax.set_ylim(0, 4)
    ax.set_xlabel(f"Position on {chrom}")
    ax.set_ylabel("Mean depth")
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()

    # Save plot
    locus_dir = os.path.join(output_dir, "locus_plots")
    os.makedirs(locus_dir, exist_ok=True)
    filename = f"{locus.cluster.replace('/', '_')}_overview.png"
    plt.savefig(os.path.join(locus_dir, filename), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Created: locus_plots/{filename}")


def plot_sample_at_locus(
    sample_id: str,
    locus: GDLocus,
    calls_df: pd.DataFrame,
    depth_df: pd.DataFrame,
    gtf: Optional[GTFParser],
    segdup: Optional[SegDupAnnotation],
    output_dir: str,
    padding: int = 50000,
):
    """
    Create detailed plot for a specific sample at a GD locus.

    Args:
        sample_id: Sample identifier
        locus: GDLocus object
        calls_df: DataFrame with CNV calls
        depth_df: DataFrame with depth data
        gtf: Optional GTFParser for gene annotations
        segdup: Optional SegDupAnnotation
        output_dir: Directory to save plots
        padding: Padding around locus boundaries
    """
    # Skip loci with no breakpoints
    if not locus.breakpoints:
        return
    
    chrom = locus.chrom
    # Filter strictly to this locus's bins (by cluster name).
    mask = (
        (depth_df["Cluster"] == locus.cluster) &
        (depth_df["Chr"] == chrom)
    )
    region_df = depth_df[mask].copy().sort_values("Start")

    if len(region_df) == 0 or sample_id not in region_df.columns:
        return

    # Bounds derived from actual data so flanks are always visible
    region_start = int(region_df["Start"].min())
    region_end = int(region_df["End"].max())

    # Get sample calls for this locus
    sample_calls = calls_df[
        (calls_df["cluster"] == locus.cluster) &
        (calls_df["sample"] == sample_id)
    ]

    # Create figure
    fig, axes = plt.subplots(3, 1, figsize=(14, 8),
                              gridspec_kw={"height_ratios": [1, 2, 1]})

    bin_starts = region_df["Start"].values
    bin_ends = region_df["End"].values
    bin_mids = (bin_starts + bin_ends) / 2
    sample_depth = region_df[sample_id].values

    # Panel 1: Annotations
    is_carrier = sample_calls["is_carrier"].any() if len(sample_calls) > 0 else False
    carrier_str = " [CARRIER]" if is_carrier else ""
    title = f"{sample_id} at {locus.cluster}{carrier_str}"
    draw_annotations_panel(axes[0], locus, region_start, region_end, chrom, title, gtf, segdup, show_gd_entries=True)
    axes[0].set_ylabel("Annotations")

    # Panel 2: Depth profile
    ax = axes[1]

    # Plot depth as bars
    bar_widths = bin_ends - bin_starts
    ax.bar(bin_mids, sample_depth, width=bar_widths * 0.9, alpha=0.7,
           color="steelblue", edgecolor="none")

    # Color bins by called CN if we have call data
    # Highlight breakpoint intervals
    for i, (start, end, name) in enumerate(locus.get_intervals()):
        interval_mask = (bin_mids >= start) & (bin_mids < end)
        ax.axvspan(start, end, alpha=0.1, color=plt.cm.Set3(i / len(locus.get_intervals())))

    # Draw reference lines
    ax.axhline(2.0, color="green", linestyle="-", alpha=0.5, linewidth=1, label="CN=2")
    ax.axhline(1.0, color="orange", linestyle="--", alpha=0.5, linewidth=1, label="CN=1")
    ax.axhline(3.0, color="purple", linestyle="--", alpha=0.5, linewidth=1, label="CN=3")

    # Draw breakpoint ranges
    for i, (bp_start, bp_end) in enumerate(locus.breakpoints):
        ax.axvspan(bp_start, bp_end, alpha=0.2, color="red", zorder=0)

    ax.set_xlim(region_start, region_end)
    ax.set_ylim(0, 5)
    ax.set_ylabel("Normalized Depth")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3, axis="y")

    # Panel 3: Call summary
    ax = axes[2]
    ax.set_xlim(region_start, region_end)
    ax.set_ylim(0, 1)

    # Show calls as colored regions
    y_pos = 0.5
    for _, call in sample_calls.iterrows():
        if call["is_carrier"]:
            color = "red" if call["svtype"] == "DEL" else "blue"
            alpha = 0.8 if call["is_best_match"] else 0.4
            ax.axvspan(call["start"], call["end"], alpha=alpha, color=color)
            label = f"{call['svtype']} CN={call['mean_cn']:.2f}"
            if call["is_best_match"]:
                label += " (best)"
            ax.text((call["start"] + call["end"]) / 2, y_pos, label,
                   ha="center", va="center", fontsize=9, fontweight="bold" if call["is_best_match"] else "normal")

    ax.set_yticks([])
    ax.set_xlabel(f"Position on {chrom}")
    ax.set_ylabel("Calls")

    plt.tight_layout()

    # Save plot
    sample_dir = os.path.join(output_dir, "sample_plots", locus.cluster.replace("/", "_"))
    os.makedirs(sample_dir, exist_ok=True)
    filename = f"{sample_id.replace('/', '_')}.png"
    plt.savefig(os.path.join(sample_dir, filename), dpi=150, bbox_inches="tight")
    plt.close()


def plot_carrier_summary(calls_df: pd.DataFrame, output_dir: str):
    """
    Create summary plots of carrier calls across all loci.

    Args:
        calls_df: DataFrame with all CNV calls
        output_dir: Directory to save plots
    """
    carriers = calls_df[calls_df["is_carrier"]].copy()

    if len(carriers) == 0:
        print("No carriers to plot.")
        return

    # Summary by locus
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))

    # Plot 1: Carriers per locus
    ax = axes[0]
    locus_counts = carriers.groupby(["cluster", "svtype"])["sample"].nunique().unstack(fill_value=0)
    locus_counts.plot(kind="barh", ax=ax, color={"DEL": "red", "DUP": "blue"}, alpha=0.7)
    ax.set_xlabel("Number of Carriers")
    ax.set_ylabel("Locus")
    ax.set_title("Carriers per GD Locus")
    ax.legend(title="SV Type")

    # Plot 2: Log probability score for carriers by SV type
    ax = axes[1]
    
    # Create common bins based on full data range
    all_scores = carriers["log_prob_score"].dropna()
    bins = np.linspace(all_scores.min(), all_scores.max(), 21)
    
    for svtype, color in [("DEL", "red"), ("DUP", "blue")]:
        sv_data = carriers[carriers["svtype"] == svtype]
        if len(sv_data) > 0:
            ax.hist(sv_data["log_prob_score"], bins=bins, alpha=0.6, color=color, 
                   label=f"{svtype} (n={len(sv_data)})", edgecolor="black")
    ax.set_xlabel("Log Probability Score")
    ax.set_ylabel("Count")
    ax.set_title("Log Probability Score Distribution in Carriers")
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "carrier_summary.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Created: carrier_summary.png")


def plot_confidence_distribution(calls_df: pd.DataFrame, output_dir: str):
    """
    Plot distribution of confidence scores.

    Args:
        calls_df: DataFrame with all CNV calls
        output_dir: Directory to save plots
    """
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))

    # Plot 1: Confidence by carrier status
    ax = axes[0]
    carriers = calls_df[calls_df["is_carrier"]]
    non_carriers = calls_df[~calls_df["is_carrier"]]
    
    # Create common bins based on full data range
    all_scores = calls_df["log_prob_score"].dropna()
    bins = np.linspace(all_scores.min(), all_scores.max(), 31)
    
    ax.hist(non_carriers["log_prob_score"], bins=bins, alpha=0.6, label="Non-carriers", 
            color="gray", edgecolor="black")
    ax.hist(carriers["log_prob_score"], bins=bins, alpha=0.6, label="Carriers",
            color="green", edgecolor="black")
    ax.set_xlabel("Log Probability Score")
    ax.set_ylabel("Count")
    ax.set_yscale("log")
    ax.legend()

    # Plot 2: Confidence vs Mean Depth
    ax = axes[1]
    
    # Plot non-carriers first (background)
    non_carriers_mask = ~calls_df["is_carrier"]
    ax.scatter(calls_df.loc[non_carriers_mask, "mean_depth"], 
               calls_df.loc[non_carriers_mask, "log_prob_score"], 
               c="gray", alpha=0.1, s=10, label="Non-carriers")
    
    # Plot carriers on top (more visible)
    carriers_mask = calls_df["is_carrier"]
    ax.scatter(calls_df.loc[carriers_mask, "mean_depth"], 
               calls_df.loc[carriers_mask, "log_prob_score"], 
               c="green", alpha=0.5, s=10, label="Carriers")
    
    ax.set_xlabel("Mean Depth (normalized)")
    ax.set_ylabel("Confidence")

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "confidence_distribution.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Created: confidence_distribution.png")


def create_carrier_pdf(
    calls_df: pd.DataFrame,
    depth_df: pd.DataFrame,
    gd_table: GDTable,
    gtf: Optional[GTFParser],
    segdup: Optional[SegDupAnnotation],
    output_dir: str,
    padding: int = 50000,
):
    """
    Create a PDF with plots for all carrier samples.

    Args:
        calls_df: DataFrame with CNV calls
        depth_df: DataFrame with depth data (from cn_posteriors)
        gd_table: GDTable with locus definitions
        gtf: Optional GTFParser
        segdup: Optional SegDupAnnotation
        output_dir: Directory to save PDF
        padding: Padding around locus boundaries
    """
    carriers = calls_df[calls_df["is_carrier"]]

    if len(carriers) == 0:
        print("No carriers to include in PDF.")
        return

    pdf_path = os.path.join(output_dir, "carrier_plots.pdf")
    print(f"Creating carrier PDF: {pdf_path}")

    with PdfPages(pdf_path) as pdf:
        # Group by locus and sample
        for cluster in sorted(carriers["cluster"].unique()):
            locus = gd_table.loci.get(cluster)
            if locus is None:
                continue
            
            # Skip loci with no breakpoints
            if not locus.breakpoints:
                print(f"  Warning: No breakpoints for {cluster}, skipping")
                continue

            cluster_carriers = carriers[carriers["cluster"] == cluster]["sample"].unique()
            print(f"  Adding {len(cluster_carriers)} carriers for {cluster}")

            for sample_id in sorted(cluster_carriers):
                # Create the plot
                chrom = locus.chrom
                # Filter strictly to this locus's bins (by cluster name).
                mask = (
                    (depth_df["Cluster"] == cluster) &
                    (depth_df["Chr"] == chrom)
                )
                region_df = depth_df[mask].sort_values("Start")

                if len(region_df) == 0 or sample_id not in region_df.columns:
                    continue

                # Bounds derived from actual data so flanks are always visible
                region_start = int(region_df["Start"].min())
                region_end = int(region_df["End"].max())

                # Get sample calls
                sample_calls = calls_df[
                    (calls_df["cluster"] == cluster) &
                    (calls_df["sample"] == sample_id)
                ]

                # Create figure
                fig, axes = plt.subplots(2, 1, figsize=(12, 4),
                                          gridspec_kw={"height_ratios": [1, 2]})

                bin_mids = (region_df["Start"].values + region_df["End"].values) / 2
                sample_depth = region_df[sample_id].values

                # Panel 1: Annotations
                carrier_call = sample_calls[sample_calls["is_carrier"]].iloc[0] if len(sample_calls[sample_calls["is_carrier"]]) > 0 else None
                call_info = f" - {carrier_call['svtype']} confidence={carrier_call['log_prob_score']:.2f}" if carrier_call is not None else ""
                title = f"{sample_id} at {cluster}{call_info}"
                draw_annotations_panel(axes[0], locus, region_start, region_end, chrom, title, gtf, segdup, show_gd_entries=False)

                # Panel 2: Depth
                ax = axes[1]
                
                # Only plot the most likely non-reference breakpoint combination
                # Filter for best match that is a carrier (non-reference)
                best_call = None
                if len(sample_calls) > 0:
                    carrier_calls = sample_calls[sample_calls["is_carrier"]]
                    if len(carrier_calls) > 0:
                        # Get the best match among carriers
                        best_matches = carrier_calls[carrier_calls["is_best_match"]]
                        if len(best_matches) > 0:
                            best_call = best_matches.iloc[0]
                        else:
                            # If no best_match flag, use highest log probability score carrier
                            best_call = carrier_calls.loc[carrier_calls["log_prob_score"].idxmax()]
                
                # Highlight the region covered by the best call
                if best_call is not None:
                    interval_start = best_call["start"]
                    interval_end = best_call["end"]
                    mean_depth = best_call["mean_depth"]
                    svtype = best_call["svtype"]
                    
                    # Color by SV type
                    color = "#FF6B6B" if svtype == "DEL" else "#6B9BD1"
                    
                    # Highlight the called region
                    ax.axvspan(interval_start, interval_end, alpha=0.2, color=color, zorder=1,
                              label=f"{svtype} region")
                    
                    # Draw line at exact mean depth
                    ax.hlines(mean_depth, interval_start, interval_end, 
                             colors="black", linewidth=2.5, alpha=0.8, zorder=2,
                             label=f"Mean depth={mean_depth:.2f}")
                
                # Plot depth as bars on top
                bar_widths = region_df["End"].values - region_df["Start"].values
                ax.bar(bin_mids, sample_depth, width=bar_widths * 0.9, alpha=0.6,
                       color="steelblue", edgecolor="none", zorder=3)

                # Reference lines
                ax.axhline(2.0, color="green", linestyle="-", alpha=0.4, linewidth=1.5, 
                          label="Reference CN=2", zorder=0)
                ax.axhline(1.0, color="orange", linestyle=":", alpha=0.4, linewidth=1, zorder=0)
                ax.axhline(3.0, color="purple", linestyle=":", alpha=0.4, linewidth=1, zorder=0)

                # Breakpoint lines
                # for i, (bp_start, bp_end) in enumerate(locus.breakpoints):
                #     ax.axvspan(bp_start, bp_end, alpha=0.2, color="red", zorder=0)

                ax.set_xlim(region_start, region_end)
                ax.set_ylim(0, 5)
                ax.set_xlabel(f"Position on {chrom}")
                ax.set_ylabel("Normalized Depth")
                ax.grid(True, alpha=0.3, axis="y", zorder=0)
                ax.legend(loc="upper right", fontsize=8)

                plt.tight_layout()
                pdf.savefig(fig)
                plt.close()

    print(f"  Saved carrier PDF: {pdf_path}")


# =============================================================================
# Main
# =============================================================================

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Call CNVs and plot GD CNV detection results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--calls", "-c",
        required=False,
        help="Pre-computed GD CNV calls file (gd_cnv_calls.tsv.gz). If not provided, calls will be made from posteriors.",
    )
    parser.add_argument(
        "--cn-posteriors",
        required=True,
        help="CN posteriors file (cn_posteriors.tsv.gz) with depth values",
    )
    parser.add_argument(
        "--bin-mappings",
        required=False,
        help="Bin mappings file (bin_mappings.tsv.gz) from gd_cnv_pyro.py. Required if --calls not provided.",
    )
    parser.add_argument(
        "--gd-table", "-g",
        required=True,
        help="GD locus definition table (TSV)",
    )
    parser.add_argument(
        "--gtf",
        required=False,
        help="Gene annotation GTF file (can be gzipped)",
    )
    parser.add_argument(
        "--segdup-bed",
        required=False,
        help="Segmental duplication BED file (can be gzipped)",
    )
    parser.add_argument(
        "--output-dir", "-o",
        required=True,
        help="Output directory for plots",
    )
    parser.add_argument(
        "--log-prob-threshold",
        type=float,
        default=-0.5,
        help="Minimum log probability score to call a CNV (only used if calling from posteriors)",
    )
    parser.add_argument(
        "--flanking-log-prob-threshold",
        type=float,
        default=-1.0,
        help="Minimum log probability score in flanking regions to classify a call as spanning (only used if calling from posteriors)",
    )
    parser.add_argument(
        "--padding",
        type=int,
        default=50000,
        help="Padding around locus boundaries (bp)",
    )
    parser.add_argument(
        "--plot-all-samples",
        action="store_true",
        help="Plot individual sample plots for all samples (not just carriers)",
    )
    parser.add_argument(
        "--sample",
        type=str,
        help="Plot only this specific sample",
    )

    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Output directory: {args.output_dir}")

    # Load data
    print("\nLoading data...")

    print(f"  Loading CN posteriors: {args.cn_posteriors}")
    cn_posteriors_df = pd.read_csv(args.cn_posteriors, sep="\t", compression="infer")
    print(f"    {len(cn_posteriors_df)} bin-sample records")
    
    print(f"  Loading GD table: {args.gd_table}")
    gd_table = GDTable(args.gd_table)
    print(f"    {len(gd_table.loci)} loci")
    
    # Call CNVs if not provided
    if args.calls:
        print(f"\n  Loading pre-computed calls: {args.calls}")
        calls_df = pd.read_csv(args.calls, sep="\t", compression="infer")
        print(f"    {len(calls_df)} call records")
    else:
        if not args.bin_mappings:
            raise ValueError("--bin-mappings required when --calls not provided")
        
        print(f"\n  Loading bin mappings: {args.bin_mappings}")
        bin_mappings_df = pd.read_csv(args.bin_mappings, sep="\t", compression="infer")
        print(f"    {len(bin_mappings_df)} bin mappings")
        
        # Call CNVs from posteriors
        calls_df = call_cnvs_from_posteriors(
            cn_posteriors_df,
            bin_mappings_df,
            gd_table,
            log_prob_threshold=args.log_prob_threshold,
            flanking_log_prob_threshold=args.flanking_log_prob_threshold,
        )
        
        # Save calls
        output_file = os.path.join(args.output_dir, "gd_cnv_calls.tsv.gz")
        calls_df.to_csv(output_file, sep="\t", index=False, compression="gzip")
        print(f"\n  Saved calls to: {output_file}")
        print(f"    {len(calls_df)} call records")
    
    # Convert cn_posteriors to depth_df format (bins x samples matrix)
    print("\n  Converting posteriors to depth matrix format...")
    depth_df = cn_posteriors_df.pivot_table(
        index=["cluster", "chr", "start", "end"],
        columns="sample",
        values="depth",
        aggfunc="first"
    ).reset_index()
    depth_df = depth_df.rename(columns={"cluster": "Cluster", "chr": "Chr", "start": "Start", "end": "End"})
    sample_cols_all = [c for c in depth_df.columns if c not in ["Cluster", "Chr", "Start", "End"]]
    print(f"    {len(depth_df)} bin-rows x {len(sample_cols_all)} samples (across all loci)")

    # Optional annotations
    gtf = None
    if args.gtf:
        print(f"\n  Loading GTF: {args.gtf}")
        gtf = GTFParser(args.gtf)

    segdup = None
    if args.segdup_bed:
        print(f"  Loading segdup BED: {args.segdup_bed}")
        segdup = SegDupAnnotation(args.segdup_bed)

    # Create summary plots
    print("\nCreating summary plots...")
    plot_carrier_summary(calls_df, args.output_dir)
    plot_confidence_distribution(calls_df, args.output_dir)

    # Create locus overview plots
    print("\nCreating locus overview plots...")
    for cluster, locus in gd_table.loci.items():
        print(f"  Processing {cluster}...")
        plot_locus_overview(
            locus, calls_df, depth_df, gtf, segdup,
            args.output_dir, padding=args.padding
        )

    # Create carrier PDF
    print("\nCreating carrier PDF...")
    create_carrier_pdf(
        calls_df, depth_df, gd_table, gtf, segdup,
        args.output_dir, padding=args.padding
    )

    # Create individual sample plots
    if args.plot_all_samples or args.sample:
        print("\nCreating individual sample plots...")
        sample_cols = [c for c in depth_df.columns if c not in ['Chr', 'Start', 'End']]
        carriers = calls_df[calls_df["is_carrier"]]
        
        for cluster, locus in gd_table.loci.items():
            if args.sample:
                samples_to_plot = [args.sample]
            elif args.plot_all_samples:
                samples_to_plot = sample_cols
            else:
                samples_to_plot = carriers[carriers["cluster"] == cluster]["sample"].unique()

            for sample_id in samples_to_plot:
                plot_sample_at_locus(
                    sample_id, locus, calls_df, depth_df, gtf, segdup,
                    args.output_dir, padding=args.padding
                )

    print("\n" + "=" * 80)
    print("Plotting complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
