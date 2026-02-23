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
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns


# =============================================================================
# Data Loading Classes
# =============================================================================

@dataclass
class GDLocus:
    """Represents a genomic disorder locus with its breakpoints."""
    cluster: str
    chrom: str
    breakpoints: List[Tuple[int, int]]  # List of breakpoint ranges (start, end)
    breakpoint_names: List[str]  # Names of breakpoints (e.g., ['1', '2', '3'] or ['A', 'B', 'C'])
    gd_entries: List[dict]
    is_nahr: bool
    is_terminal: bool

    @property
    def start(self) -> int:
        """Get overall start position (min of all breakpoint starts)."""
        return min(bp[0] for bp in self.breakpoints)

    @property
    def end(self) -> int:
        """Get overall end position (max of all breakpoint ends)."""
        return max(bp[1] for bp in self.breakpoints)

    def get_intervals(self) -> List[Tuple[int, int, str]]:
        """Get all intervals between adjacent breakpoints."""
        intervals = []
        for i in range(len(self.breakpoints) - 1):
            # Interval spans from end of BP_i to start of BP_{i+1}
            start = self.breakpoints[i][1]  # End of current breakpoint
            end = self.breakpoints[i + 1][0]  # Start of next breakpoint
            # Use actual breakpoint names from the table
            interval_name = f"{self.breakpoint_names[i]}-{self.breakpoint_names[i + 1]}"
            intervals.append((start, end, interval_name))
        return intervals


class GDTable:
    """Parser for GD locus definitions."""

    def __init__(self, filepath: str):
        self.df = pd.read_csv(filepath, sep="\t")
        self._validate_columns()
        self.loci = self._parse_loci()

    def _validate_columns(self):
        """Validate required columns exist."""
        required_cols = ["chr", "start_GRCh38", "end_GRCh38", "GD_ID", "svtype",
                         "NAHR", "terminal", "cluster", "BP1", "BP2"]
        missing = set(required_cols) - set(self.df.columns)
        if missing:
            raise ValueError(f"Missing required columns in GD table: {missing}")

    def _parse_loci(self) -> Dict[str, GDLocus]:
        """Parse GD table into loci with breakpoint ranges using BP1 and BP2 columns."""
        loci = {}
        for _, row in self.df.iterrows():
            cluster = row["cluster"]
            
            # Skip entries without a cluster (standalone regions)
            if pd.isna(cluster) or cluster == "":
                continue

            if cluster not in loci:
                loci[cluster] = {
                    "chrom": row["chr"],
                    "breakpoint_coords": {},  # Dict mapping BP name to list of coordinates
                    "entries": [],
                    "is_nahr": row["NAHR"] == "yes",
                    "is_terminal": row["terminal"] != "no" if pd.notna(row["terminal"]) else False,
                }

            # Get breakpoint names from BP1 and BP2 columns
            bp1 = row["BP1"]
            bp2 = row["BP2"]
            
            # Skip if BP1 or BP2 is missing
            if pd.isna(bp1) or pd.isna(bp2):
                continue
            
            # Convert to strings for consistent handling
            bp1 = str(bp1)
            bp2 = str(bp2)
            
            # Track coordinates for each breakpoint
            if bp1 not in loci[cluster]["breakpoint_coords"]:
                loci[cluster]["breakpoint_coords"][bp1] = []
            if bp2 not in loci[cluster]["breakpoint_coords"]:
                loci[cluster]["breakpoint_coords"][bp2] = []
            
            loci[cluster]["breakpoint_coords"][bp1].append(row["start_GRCh38"])
            loci[cluster]["breakpoint_coords"][bp2].append(row["end_GRCh38"])

            loci[cluster]["entries"].append({
                "GD_ID": row["GD_ID"],
                "start_GRCh38": row["start_GRCh38"],
                "end_GRCh38": row["end_GRCh38"],
                "svtype": row["svtype"],
                "BP1": bp1,
                "BP2": bp2,
            })

        result = {}
        for cluster, data in loci.items():
            # Convert breakpoint coordinates to ranges
            bp_ranges = []
            
            # Sort breakpoint names (numeric if possible, otherwise alphabetic)
            def sort_key(bp_name):
                try:
                    return (0, int(bp_name))
                except ValueError:
                    return (1, bp_name)
            
            bp_names_sorted = sorted(data["breakpoint_coords"].keys(), key=sort_key)
            
            for bp_name in bp_names_sorted:
                coords = data["breakpoint_coords"][bp_name]
                bp_range = (min(coords), max(coords))
                bp_ranges.append(bp_range)
            
            result[cluster] = GDLocus(
                cluster=cluster,
                chrom=data["chrom"],
                breakpoints=bp_ranges,
                breakpoint_names=bp_names_sorted,
                gd_entries=data["entries"],
                is_nahr=data["is_nahr"],
                is_terminal=data["is_terminal"],
            )
        return result


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
# Plotting Functions
# =============================================================================

def get_sample_columns(df: pd.DataFrame) -> list:
    """Extract sample column names from DataFrame."""
    metadata_cols = ["Chr", "Start", "End", "source_file", "Bin"]
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
        show_gd_entries: Whether to show GD entry regions
    """
    ax.set_xlim(region_start, region_end)
    ax.set_ylim(0, 1)
    ax.set_title(title)

    # Draw GD entry regions if requested
    if show_gd_entries:
        for entry in locus.gd_entries:
            color = "red" if entry["svtype"] == "DEL" else "blue"
            ax.axvspan(entry["start_GRCh38"], entry["end_GRCh38"],
                      ymin=0.7, ymax=0.9, alpha=0.3, color=color,
                      label=f"{entry['GD_ID']} ({entry['svtype']})")

    # Draw breakpoint ranges
    for i, (bp_start, bp_end) in enumerate(locus.breakpoints):
        y_min = 0.0 if show_gd_entries else 0.8
        y_max = 1.0
        ax.axvspan(bp_start, bp_end, ymin=y_min, ymax=y_max, alpha=0.2 if show_gd_entries else 0.3, color="red")
        ax.text((bp_start + bp_end) / 2, 0.95 if show_gd_entries else 0.9, f"BP{i+1}", ha="center", va="top" if show_gd_entries else "center", 
                fontsize=8, fontweight="bold", color="darkred")

    # Draw segdup regions
    if segdup:
        sd_regions = segdup.get_regions_in_range(chrom, region_start, region_end)
        for sd_start, sd_end in sd_regions:
            ax.axvspan(max(sd_start, region_start), min(sd_end, region_end),
                      ymin=0.4 if show_gd_entries else 0.5, ymax=0.6 if show_gd_entries else 0.7, 
                      alpha=0.5, color="orange")

    # Draw genes
    if gtf:
        genes = gtf.get_genes_in_region(chrom, region_start, region_end,
                                         gene_types=["protein_coding"])
        min_gene_size = 20000  # Minimum 20kb to show label
        max_genes = 8
        for gene in genes[:max_genes]:
            gene_start = max(gene["start"], region_start)
            gene_end = min(gene["end"], region_end)
            ax.hlines(0.2, gene_start, gene_end, colors="darkblue", linewidth=4)
            # Only label genes larger than minimum size
            if (gene_end - gene_start) >= min_gene_size:
                ax.text((gene_start + gene_end) / 2, 0.18,
                       gene["gene_name"], ha="center", va="top", fontsize=7, rotation=45)

    ax.set_yticks([])
    if not show_gd_entries:
        ax.set_ylabel("Annotations")


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
        depth_df: DataFrame with depth data
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
    region_start = locus.start - padding
    region_end = locus.end + padding

    # Get bins in this region
    mask = (
        (depth_df["Chr"] == chrom) &
        (depth_df["End"] > region_start) &
        (depth_df["Start"] < region_end)
    )
    region_df = depth_df[mask].copy()

    if len(region_df) == 0:
        print(f"  Warning: No bins found for locus {locus.cluster}")
        return

    sample_cols = get_sample_columns(region_df)
    n_samples = len(sample_cols)

    # Get carriers for this locus
    carriers = calls_df[
        (calls_df["cluster"] == locus.cluster) &
        (calls_df["is_carrier"])
    ]["sample"].unique()

    # Create figure
    fig_height = 4 + 0.4 * min(n_samples, 20)  # Scale with number of samples
    fig, axes = plt.subplots(3, 1, figsize=(14, fig_height), 
                              gridspec_kw={"height_ratios": [1, 4, 1]})

    # X-axis: genomic position
    bin_mids = (region_df["Start"].values + region_df["End"].values) / 2

    # Panel 1: Gene annotations and segdup regions
    ax = axes[0]
    ax.set_xlim(region_start, region_end)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Annotations")
    ax.set_title(f"{locus.cluster} ({chrom}:{locus.start:,}-{locus.end:,})")

    # Draw breakpoint ranges as shaded regions
    for i, (bp_start, bp_end) in enumerate(locus.breakpoints):
        ax.axvspan(bp_start, bp_end, ymin=0.85, ymax=1.0, alpha=0.4, color="red", zorder=10)
        ax.text((bp_start + bp_end) / 2, 0.925, f"BP{i+1}", ha="center", va="center", 
                fontsize=7, fontweight="bold", color="darkred")

    # Draw breakpoint intervals
    intervals = locus.get_intervals()
    colors = plt.cm.Set3(np.linspace(0, 1, len(intervals)))
    for i, (start, end, name) in enumerate(intervals):
        ax.axvspan(start, end, alpha=0.2, color=colors[i], label=name)
        ax.text((start + end) / 2, 0.75, name, ha="center", va="top", fontsize=8)

    # Draw segdup regions
    if segdup:
        sd_regions = segdup.get_regions_in_range(chrom, region_start, region_end)
        for sd_start, sd_end in sd_regions:
            ax.axvspan(max(sd_start, region_start), min(sd_end, region_end),
                      ymin=0.45, ymax=0.65, alpha=0.5, color="orange")
        if sd_regions:
            ax.text(region_start + 1000, 0.55, "SegDup", fontsize=7, va="center")

    # Draw genes
    if gtf:
        genes = gtf.get_genes_in_region(chrom, region_start, region_end,
                                         gene_types=["protein_coding"])
        y_pos = 0.3
        min_gene_size = 20000  # Minimum 20kb to show label
        for gene in genes[:10]:  # Limit to 10 genes
            gene_start = max(gene["start"], region_start)
            gene_end = min(gene["end"], region_end)
            ax.hlines(y_pos, gene_start, gene_end, colors="blue", linewidth=3)
            # Only label genes larger than minimum size
            if (gene_end - gene_start) >= min_gene_size:
                ax.text((gene_start + gene_end) / 2, y_pos - 0.05,
                       gene["gene_name"], ha="center", va="top", fontsize=7)

    ax.set_yticks([])

    # Panel 2: Depth heatmap across samples
    ax = axes[1]

    # Create depth matrix
    depth_matrix = region_df[sample_cols].values.T  # samples x bins

    # Sort samples: carriers first
    carrier_mask = np.array([s in carriers for s in sample_cols])
    sort_idx = np.argsort(~carrier_mask)  # Carriers first
    depth_matrix = depth_matrix[sort_idx]

    # Plot heatmap
    im = ax.imshow(depth_matrix, aspect="auto", cmap="RdBu_r",
                   vmin=0, vmax=4, interpolation="nearest",
                   extent=[region_start, region_end, n_samples, 0])

    # Mark carriers section with dividing line and bracket
    n_carriers = carrier_mask.sum()
    if n_carriers > 0 and n_carriers < n_samples:
        # Draw thin horizontal line separating carriers from non-carriers
        # ax.axhline(n_carriers, color="black", linewidth=1, linestyle="-", zorder=10)
        
        # Add bracket on the left to label carriers section
        # Use axes coordinates for positioning outside the plot
        # Since y=0 is at TOP and y=n_samples is at BOTTOM in data coordinates,
        # and carriers occupy rows 0 to n_carriers-1 (top of heatmap)
        bracket_y_top = 1.0  # Top of axes (y=0 in data)
        bracket_y_bottom = 1.0 - (n_carriers / n_samples)  # Bottom of carriers section
        
        # Draw bracket (left facing) in axes coordinates
        bracket_x = -0.02  # Position to the left of the y-axis
        bracket_width = 0.01
        
        # Top horizontal
        ax.plot([bracket_x, bracket_x - bracket_width], [bracket_y_top, bracket_y_top], 
                'k-', linewidth=1.5, transform=ax.transAxes, clip_on=False)
        # Vertical
        ax.plot([bracket_x - bracket_width, bracket_x - bracket_width], [bracket_y_top, bracket_y_bottom], 
                'k-', linewidth=1.5, transform=ax.transAxes, clip_on=False)
        # Bottom horizontal
        ax.plot([bracket_x - bracket_width, bracket_x], [bracket_y_bottom, bracket_y_bottom], 
                'k-', linewidth=1.5, transform=ax.transAxes, clip_on=False)
        
        # Add label in axes coordinates
        ax.text(bracket_x - bracket_width - 0.01, (bracket_y_top + bracket_y_bottom) / 2, 
                f"Carriers\n(n={n_carriers})", 
                rotation=90, va="center", ha="right", fontsize=8, fontweight="bold",
                transform=ax.transAxes)

    # Draw breakpoint ranges
    for i, (bp_start, bp_end) in enumerate(locus.breakpoints):
        ax.axvspan(bp_start, bp_end, alpha=0.1, color="red", zorder=5)

    ax.set_ylabel(f"Samples (n={n_samples})")
    ax.set_yticks([])

    # Add colorbar below the plot
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.08, aspect=40)
    cbar.set_label("Copy Number")

    # Panel 3: Mean depth profile with carrier/non-carrier separation
    ax = axes[2]

    # Calculate mean depth for carriers and non-carriers
    if n_carriers > 0 and n_carriers < n_samples:
        # Plot non-carriers
        non_carrier_depths = region_df[[s for s in sample_cols if s not in carriers]].mean(axis=1)
        ax.plot(bin_mids, non_carrier_depths.values, "b-", linewidth=1, alpha=0.7,
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
                group_depths = region_df[valid_samples].mean(axis=1)
                
                # Only plot if there's actual data
                if len(group_depths) > 0 and not group_depths.isna().all():
                    # Choose color based on svtype
                    if svtype == "DEL":
                        color = del_colors[color_idx % len(del_colors)]
                    else:  # DUP
                        color = dup_colors[color_idx % len(dup_colors)]
                    ax.plot(bin_mids, group_depths.values, "-", linewidth=1.5, alpha=0.8,
                           color=color, label=f"{svtype} {bp_interval} (n={len(valid_samples)})")
                    color_idx += 1
                    plotted_any = True
        
        # Only show legend if we actually plotted something
        if plotted_any:
            ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.35), 
                     ncol=6, fontsize=8, frameon=True, fancybox=True)
    else:
        mean_depth = region_df[sample_cols].mean(axis=1)
        ax.plot(bin_mids, mean_depth.values, "b-", linewidth=1.5)

    # Draw breakpoint ranges
    for i, (bp_start, bp_end) in enumerate(locus.breakpoints):
        ax.axvspan(bp_start, bp_end, alpha=0.1, color="red", zorder=0)

    ax.axhline(2.0, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlim(region_start, region_end)
    ax.set_ylim(0, 4)
    ax.set_xlabel(f"Position on {chrom}")
    ax.set_ylabel("Mean CN")
    ax.grid(True, alpha=0.3)

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
    region_start = locus.start - padding
    region_end = locus.end + padding

    # Get bins in this region
    mask = (
        (depth_df["Chr"] == chrom) &
        (depth_df["End"] > region_start) &
        (depth_df["Start"] < region_end)
    )
    region_df = depth_df[mask].copy()

    if len(region_df) == 0 or sample_id not in region_df.columns:
        return

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
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Plot 1: Carriers per locus
    ax = axes[0]
    locus_counts = carriers.groupby(["cluster", "svtype"])["sample"].nunique().unstack(fill_value=0)
    locus_counts.plot(kind="barh", ax=ax, color={"DEL": "red", "DUP": "blue"}, alpha=0.7)
    ax.set_xlabel("Number of Carriers")
    ax.set_ylabel("Locus")
    ax.set_title("Carriers per GD Locus")
    ax.legend(title="SV Type")

    # Plot 2: Mean CN distribution for carriers
    ax = axes[1]
    for svtype, color in [("DEL", "red"), ("DUP", "blue")]:
        sv_data = carriers[carriers["svtype"] == svtype]
        if len(sv_data) > 0:
            ax.hist(sv_data["mean_cn"], bins=20, alpha=0.6, color=color, 
                   label=f"{svtype} (n={len(sv_data)})", edgecolor="black")
    ax.axvline(2.0, color="gray", linestyle="--", alpha=0.5, label="CN=2")
    ax.set_xlabel("Mean Copy Number")
    ax.set_ylabel("Count")
    ax.set_title("Copy Number Distribution in Carriers")
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
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1: Confidence by carrier status
    ax = axes[0]
    carriers = calls_df[calls_df["is_carrier"]]
    non_carriers = calls_df[~calls_df["is_carrier"]]
    ax.hist(non_carriers["confidence"], bins=30, alpha=0.6, label="Non-carriers", 
            color="gray", edgecolor="black")
    ax.hist(carriers["confidence"], bins=30, alpha=0.6, label="Carriers",
            color="green", edgecolor="black")
    ax.set_xlabel("Confidence")
    ax.set_ylabel("Count")
    ax.set_title("Confidence Distribution")
    ax.legend()

    # Plot 2: Confidence vs Mean CN
    ax = axes[1]
    colors = calls_df["is_carrier"].map({True: "green", False: "gray"})
    ax.scatter(calls_df["mean_cn"], calls_df["confidence"], c=colors, alpha=0.5, s=20)
    ax.axvline(1.5, color="red", linestyle="--", alpha=0.5, label="DEL threshold")
    ax.axvline(2.5, color="blue", linestyle="--", alpha=0.5, label="DUP threshold")
    ax.set_xlabel("Mean Copy Number")
    ax.set_ylabel("Confidence")
    ax.set_title("Confidence vs Copy Number")
    ax.legend()

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
        depth_df: DataFrame with depth data
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
                region_start = locus.start - padding
                region_end = locus.end + padding

                mask = (
                    (depth_df["Chr"] == chrom) &
                    (depth_df["End"] > region_start) &
                    (depth_df["Start"] < region_end)
                )
                region_df = depth_df[mask]

                if len(region_df) == 0 or sample_id not in region_df.columns:
                    continue

                # Get sample calls
                sample_calls = calls_df[
                    (calls_df["cluster"] == cluster) &
                    (calls_df["sample"] == sample_id)
                ]

                # Create figure
                fig, axes = plt.subplots(2, 1, figsize=(14, 8),
                                          gridspec_kw={"height_ratios": [1, 2]})

                bin_mids = (region_df["Start"].values + region_df["End"].values) / 2
                sample_depth = region_df[sample_id].values

                # Panel 1: Annotations
                carrier_call = sample_calls[sample_calls["is_carrier"]].iloc[0] if len(sample_calls[sample_calls["is_carrier"]]) > 0 else None
                call_info = f" - {carrier_call['svtype']} CN={carrier_call['mean_cn']:.2f} confidence={carrier_call['confidence']:.2f}" if carrier_call is not None else ""
                title = f"{sample_id} at {cluster}{call_info}"
                draw_annotations_panel(axes[0], locus, region_start, region_end, chrom, title, gtf, segdup, show_gd_entries=False)

                # Panel 2: Depth
                ax = axes[1]
                
                # First, draw background shading for called copy states
                cn_colors = {
                    0: "#8B0000",  # Dark red for CN=0
                    1: "#FF6B6B",  # Light red for CN=1
                    2: "#90EE90",  # Light green for CN=2 (reference)
                    3: "#6B9BD1",  # Light blue for CN=3
                    4: "#4169E1",  # Royal blue for CN=4+
                }
                
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
                            # If no best_match flag, use highest confidence carrier
                            best_call = carrier_calls.loc[carrier_calls["confidence"].idxmax()]
                
                # Draw called CN as horizontal band for the best non-reference call only
                if best_call is not None:
                    interval_start = best_call["start"]
                    interval_end = best_call["end"]
                    mean_cn = best_call["mean_cn"]
                    
                    # Determine most likely CN state from probabilities
                    prob_cols = [col for col in best_call.index if col.startswith("prob_cn_")]
                    if prob_cols:
                        # Extract CN states and their probabilities
                        cn_probs = {}
                        for col in prob_cols:
                            cn = int(col.replace("prob_cn_", ""))
                            cn_probs[cn] = best_call[col]
                        # Find CN with highest probability
                        cn_state = max(cn_probs.items(), key=lambda x: x[1])[0]
                    else:
                        # Fallback: round mean CN
                        cn_state = int(round(mean_cn))
                    
                    cn_state = min(cn_state, 4)  # Cap at CN=4
                    
                    # Get color for this CN state
                    color = cn_colors.get(cn_state, "#CCCCCC")
                    
                    # Draw horizontal band at the called CN level
                    ax.axhspan(max(0, cn_state - 0.15), min(5, cn_state + 0.15),
                              xmin=(interval_start - region_start) / (region_end - region_start),
                              xmax=(interval_end - region_start) / (region_end - region_start),
                              alpha=0.3, color=color, zorder=1)
                    
                    # Draw line at exact mean CN
                    ax.hlines(mean_cn, interval_start, interval_end, 
                             colors="black", linewidth=2.5, alpha=0.8, zorder=2,
                             label=f"Called CN={cn_state} (mean={mean_cn:.2f})")
                
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
                for i, (bp_start, bp_end) in enumerate(locus.breakpoints):
                    ax.axvspan(bp_start, bp_end, alpha=0.2, color="red", zorder=0)

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
        description="Plot GD CNV detection results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--calls", "-c",
        required=True,
        help="GD CNV calls file (gd_cnv_calls.tsv.gz)",
    )
    parser.add_argument(
        "--depth-data", "-d",
        required=True,
        help="Binned read depth data (TSV)",
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

    print(f"  Loading calls: {args.calls}")
    calls_df = pd.read_csv(args.calls, sep="\t", compression="infer")
    print(f"    {len(calls_df)} call records")

    print(f"  Loading depth data: {args.depth_data}")
    depth_df = pd.read_csv(args.depth_data, sep="\t", compression="infer")
    if "#Chr" in depth_df.columns:
        depth_df["Chr"] = depth_df["#Chr"]
        depth_df = depth_df.drop("#Chr", axis=1)
    print(f"    {len(depth_df)} bins")

    print(f"  Loading GD table: {args.gd_table}")
    gd_table = GDTable(args.gd_table)
    print(f"    {len(gd_table.loci)} loci")

    # Optional annotations
    gtf = None
    if args.gtf:
        print(f"  Loading GTF: {args.gtf}")
        gtf = GTFParser(args.gtf)

    segdup = None
    if args.segdup_bed:
        print(f"  Loading segdup BED: {args.segdup_bed}")
        segdup = SegDupAnnotation(args.segdup_bed)

    # Normalize depth data
    print("\nNormalizing depth data...")
    sample_cols = get_sample_columns(depth_df)
    autosome_mask = ~depth_df["Chr"].isin(["chrX", "chrY"])
    if autosome_mask.any():
        column_medians = np.median(depth_df.loc[autosome_mask, sample_cols], axis=0)
    else:
        column_medians = np.median(depth_df[sample_cols], axis=0)
    depth_df[sample_cols] = 2.0 * depth_df[sample_cols] / column_medians[np.newaxis, :]

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
