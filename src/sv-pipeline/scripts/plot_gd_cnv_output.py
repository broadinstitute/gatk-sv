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
import sys
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
# Logging
# =============================================================================


class TeeStream:
    """Write to both the original stream and a log file."""

    def __init__(self, original_stream, log_file):
        self.original_stream = original_stream
        self.log_file = log_file

    def write(self, message: str):
        self.original_stream.write(message)
        self.log_file.write(message)

    def flush(self):
        self.original_stream.flush()
        self.log_file.flush()

    # Forward any other attribute lookups to the original stream so that
    # callers relying on e.g. ``fileno()`` or ``isatty()`` still work.
    def __getattr__(self, name):
        return getattr(self.original_stream, name)


def setup_logging(output_dir: str, filename: str = "log.txt"):
    """Redirect stdout and stderr so output goes to both the console and a log file.

    Call once at the beginning of ``main()``.  The log file is opened in
    write mode so each run produces a fresh log.

    Args:
        output_dir: Directory where the log file will be created.
        filename: Name of the log file (default ``log.txt``).

    Returns:
        The open file handle (kept so the caller can close it if desired).
    """
    log_path = os.path.join(output_dir, filename)
    log_fh = open(log_path, "w")
    sys.stdout = TeeStream(sys.__stdout__, log_fh)
    sys.stderr = TeeStream(sys.__stderr__, log_fh)
    return log_fh


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
        Dict mapping interval name to list of array_idx values
            (positional indices into the per-sample posterior arrays)
    """
    locus_bins = bin_mappings_df[bin_mappings_df["cluster"] == cluster]
    interval_bins = {}
    for interval_name, group in locus_bins.groupby("interval"):
        # Use the explicit array_idx column – this is the positional index
        # into the per-sample posterior arrays, independent of DataFrame
        # index state.
        interval_bins[interval_name] = group["array_idx"].tolist()
    return interval_bins


def compute_interval_cn_stats(
    cn_posteriors_df: pd.DataFrame,
    interval_bins: Dict[str, List[int]],
    sample_id: str,
) -> Dict[str, dict]:
    """
    Compute copy number statistics for each interval for a specific sample.
    
    Args:
        cn_posteriors_df: DataFrame with CN posterior probabilities
        interval_bins: Dict mapping interval name to list of array_idx values
            (positional indices into the per-sample posterior rows)
        sample_id: Sample identifier
    
    Returns:
        Dict mapping interval name to CN statistics
    """
    result = {}
    
    # Filter posteriors to this sample and reset index so that
    # .iloc[k] reliably gives the k-th bin (matching array_idx).
    sample_posteriors = (
        cn_posteriors_df[cn_posteriors_df["sample"] == sample_id]
        .reset_index(drop=True)
    )
    
    # Get probability columns (prob_cn_0, prob_cn_1, etc.)
    prob_cols = [c for c in sample_posteriors.columns if c.startswith("prob_cn_")]
    n_states = len(prob_cols)
    
    for interval_name, bin_indices in interval_bins.items():
        # bin_indices are array_idx values – use .iloc for positional access
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
    ploidy: int = 2,
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
        ploidy: Expected copy number for this sample on this chromosome
    
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
        
        # Determine which intervals this entry spans using BP1 and BP2.
        # Use the locus's authoritative breakpoint ordering directly
        # instead of reverse-engineering from interval name strings
        # (which breaks when BP names contain hyphens).
        covered_tuples = locus.get_intervals_between(bp1, bp2)
        covered_intervals = [name for _, _, name in covered_tuples]
        
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
                    p_del = max(probs[:ploidy].sum(), 1e-5)  # P(CN < ploidy)
                    if p_del > 0:
                        log_prob_score += weight * np.log(p_del)
                        total_weight += weight
        
        elif svtype == "DUP":
            for interval in covered_intervals:
                if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                    probs = interval_stats[interval]["cn_probs"]
                    weight = interval_stats[interval]["n_bins"]
                    p_dup = max(probs[ploidy + 1:].sum(), 1e-5)  # P(CN > ploidy)
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
                    p_flank = max(flank_probs[:ploidy].sum(), 1e-5)
                elif svtype == "DUP":
                    p_flank = max(flank_probs[ploidy + 1:].sum(), 1e-5)
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
            "BP1": bp1,
            "BP2": bp2,
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
    ploidy: int = 2,
) -> Dict[str, Optional[str]]:
    """
    Determine the most likely breakpoint pair for a GD CNV, separately for DEL and DUP.
    
    For loci with multiple possible breakpoint configurations (e.g., BP1-2, BP1-3),
    determine which best fits the observed data for each svtype.
    
    Args:
        locus: GDLocus object
        interval_stats: Dict mapping interval name to CN statistics
        calls: List of CNV call dictionaries
        ploidy: Expected copy number for this sample on this chromosome
    
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
        best_size = -1
        best_gd_id = None
        
        for call in carrier_calls:
            covered_intervals = set(call["intervals"])
            uncovered_intervals = all_intervals - covered_intervals
            
            score = 0.0
            total_weight = 0.0
            
            if svtype == "DEL":
                # Score affected intervals: P(CN < ploidy)
                for interval in covered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = 1.  # interval_stats[interval]["n_bins"]
                        p_del = probs[:ploidy].sum()
                        if p_del > 0:
                            score += weight * np.log(p_del)
                            total_weight += weight
                
                # Score unaffected intervals: P(CN >= ploidy)
                for interval in uncovered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = 1.  # interval_stats[interval]["n_bins"]
                        p_normal = probs[ploidy:].sum()
                        if p_normal > 0:
                            score += weight * np.log(p_normal)
                            total_weight += weight
            
            else:  # DUP
                # Score affected intervals: P(CN > ploidy)
                for interval in covered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = 1.  # interval_stats[interval]["n_bins"]
                        p_dup = probs[ploidy + 1:].sum()
                        if p_dup > 0:
                            score += weight * np.log(p_dup)
                            total_weight += weight
                
                # Score unaffected intervals: P(CN <= ploidy)
                for interval in uncovered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = 1.  # interval_stats[interval]["n_bins"]
                        p_normal = probs[:ploidy + 1].sum()
                        if p_normal > 0:
                            score += weight * np.log(p_normal)
                            total_weight += weight
            
            # Normalize
            if total_weight > 0:
                score = score / total_weight
            
            # Tiebreaker: when scores are equal, prefer the larger variant
            call_size = call["end"] - call["start"]
            if score > best_score or (score == best_score and call_size > best_size):
                best_score = score
                best_size = call_size
                best_gd_id = call["GD_ID"]
        
        best_by_svtype[svtype] = best_gd_id
    
    return best_by_svtype


def call_cnvs_from_posteriors(
    cn_posteriors_df: pd.DataFrame,
    bin_mappings_df: pd.DataFrame,
    gd_table: GDTable,
    log_prob_threshold: float = -0.5,
    flanking_log_prob_threshold: float = -1.0,
    ploidy_df: Optional[pd.DataFrame] = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Call CNVs from posterior probabilities.
    
    Args:
        cn_posteriors_df: DataFrame with CN posterior probabilities per bin/sample
        bin_mappings_df: DataFrame with bin-to-interval mappings
        gd_table: GDTable with locus definitions
        log_prob_threshold: Minimum log probability score to call a CNV
        flanking_log_prob_threshold: Minimum log probability score in flanking regions to classify as spanning
        ploidy_df: Optional DataFrame with columns (sample, contig, ploidy).
            If None, ploidy=2 is assumed for all sample/contig pairs.
        verbose: If True, print per-sample log probability scores for every
            GD entry at every locus, including flanking region scores.
    
    Returns:
        DataFrame with CNV calls for all samples and loci
    """
    print("\n" + "=" * 80)
    print("CALLING CNVs FROM POSTERIORS")
    print("=" * 80)
    
    all_results = []
    sample_ids = cn_posteriors_df["sample"].unique()
    
    # Build ploidy lookup: (sample, contig) -> ploidy, default 2
    ploidy_lookup: Dict[Tuple[str, str], int] = {}
    if ploidy_df is not None:
        for _, row in ploidy_df.iterrows():
            ploidy_lookup[(str(row["sample"]), str(row["contig"]))] = int(row["ploidy"])
        print(f"  Loaded ploidy for {len(ploidy_lookup)} sample/contig pairs")
    else:
        print("  No ploidy table provided; assuming diploid (ploidy=2) everywhere")
    
    # ---- Validate alignment between cn_posteriors and bin_mappings ----
    n_bins = len(bin_mappings_df)
    n_samples = len(sample_ids)
    expected_rows = n_bins * n_samples
    if len(cn_posteriors_df) != expected_rows:
        print(f"  WARNING: cn_posteriors has {len(cn_posteriors_df)} rows, "
              f"expected {expected_rows} ({n_bins} bins × {n_samples} samples)")

    # Spot-check that the first sample's bin coordinates match bin_mappings
    first_sample = sample_ids[0]
    first_sample_rows = cn_posteriors_df[cn_posteriors_df["sample"] == first_sample]
    if len(first_sample_rows) == n_bins:
        post_coords = list(zip(
            first_sample_rows["chr"].values,
            first_sample_rows["start"].values,
            first_sample_rows["end"].values,
        ))
        map_coords = list(zip(
            bin_mappings_df["chr"].values,
            bin_mappings_df["start"].values,
            bin_mappings_df["end"].values,
        ))
        if post_coords != map_coords:
            raise ValueError(
                "Bin coordinates in cn_posteriors do not match bin_mappings. "
                "The files may have been generated by different runs or with "
                "a version that had a bin-ordering bug. Please re-run "
                "gd_cnv_pyro.py to regenerate both files."
            )
        print(f"  Validated: bin coordinates match between posteriors and mappings ({n_bins} bins)")
    else:
        print(f"  WARNING: sample {first_sample} has {len(first_sample_rows)} bins, "
              f"expected {n_bins}; skipping alignment check")

    # Pre-group posteriors by sample for faster lookups.
    # Reset each sample's index to 0..N-1 so that .iloc[k] reliably
    # selects the k-th bin, matching the array_idx values from
    # bin_mappings_df.
    print("  Organizing data for fast access...")
    posteriors_by_sample = {}
    for sample_id in sample_ids:
        sample_data = (
            cn_posteriors_df[cn_posteriors_df["sample"] == sample_id]
            .reset_index(drop=True)
        )
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
            
            # Look up ploidy for this sample/chromosome
            sample_ploidy = ploidy_lookup.get((str(sample_id), locus.chrom), 2)
            
            # Call GD CNVs
            calls = call_gd_cnv(locus, interval_stats, log_prob_threshold, flanking_log_prob_threshold, ploidy=sample_ploidy)
            
            # Determine best breakpoints
            best_by_svtype = determine_best_breakpoints(locus, interval_stats, calls, ploidy=sample_ploidy)

            # ---- verbose per-sample logging ----
            if verbose:
                for call in calls:
                    tag_parts = []
                    if call["is_carrier"]:
                        tag_parts.append("CARRIER")
                    if call.get("is_spanning"):
                        tag_parts.append("SPANNING")
                    tag = " [" + ",".join(tag_parts) + "]" if tag_parts else ""
                    flank_str = (
                        f"{call['flanking_log_prob_score']:.4f}"
                        if not np.isnan(call.get("flanking_log_prob_score", np.nan))
                        else "NA"
                    )
                    print(
                        f"    {sample_id:30s}  {call['GD_ID']:25s}  "
                        f"{call['svtype']:4s}  ploidy={sample_ploidy}  "
                        f"logP={call['log_prob_score']:+.4f}  "
                        f"flank={flank_str}  "
                        f"intervals={','.join(call['intervals'])}{tag}"
                    )

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
                    "BP1": call["BP1"],
                    "BP2": call["BP2"],
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
    min_gene_label_spacing: float = 0.05,
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
        min_gene_label_spacing: Minimum distance between gene label centres as a
            fraction of the plot width.  Labels closer than this threshold are
            suppressed to avoid overlap (lines are always drawn).
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
        ax.text((bp_start + bp_end) / 2, 0.02, f"BP {bp_name}", ha="center", va="center", 
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
        plot_width = region_end - region_start
        min_spacing = min_gene_label_spacing * plot_width
        last_labeled_center = -np.inf
        for gene in genes:
            gene_start = max(gene["start"], region_start)
            gene_end = min(gene["end"], region_end)
            gene_center = (gene_start + gene_end) / 2
            ax.hlines(y_pos, gene_start, gene_end, colors="blue", linewidth=4, zorder=3)
            # Only label genes whose centre is far enough from the last labelled gene
            if gene_center - last_labeled_center >= min_spacing:
                ax.text(gene_center, y_pos - 0.05,
                       gene["gene_name"], ha="center", va="bottom", fontsize=7, zorder=3)
                last_labeled_center = gene_center

    ax.set_yticks([])


def _build_ploidy_lookup(
    ploidy_df: Optional[pd.DataFrame],
) -> Dict[Tuple[str, str], int]:
    """Build (sample, contig) -> ploidy mapping from a ploidy DataFrame."""
    lookup: Dict[Tuple[str, str], int] = {}
    if ploidy_df is not None:
        for _, row in ploidy_df.iterrows():
            lookup[(str(row["sample"]), str(row["contig"]))] = int(row["ploidy"])
    return lookup


def _sort_samples_by_ploidy(
    sample_cols: List[str],
    chrom: str,
    ploidy_lookup: Dict[Tuple[str, str], int],
) -> List[str]:
    """Return *sample_cols* sorted by ploidy (ascending), then by name."""
    return sorted(sample_cols, key=lambda s: (ploidy_lookup.get((s, chrom), 2), s))


def _build_gap_positions(
    bin_mids: np.ndarray,
    bin_starts: np.ndarray,
    bin_ends: np.ndarray,
) -> np.ndarray:
    """Build a positions array with NaN sentinels at bin gaps."""
    bin_spacing = np.median(bin_starts[1:] - bin_ends[:-1]) if len(bin_starts) > 1 else 0
    gap_threshold = 3 * bin_spacing

    positions: List[float] = []
    for i in range(len(bin_starts)):
        positions.append(bin_mids[i])
        if i < len(bin_starts) - 1:
            if bin_starts[i + 1] - bin_ends[i] > gap_threshold:
                positions.append(np.nan)
    return np.array(positions)


def _depths_with_gaps(
    raw_depths: np.ndarray,
    plot_positions: np.ndarray,
) -> List[float]:
    """Map per-bin depths onto the gap-aware position array."""
    result: List[float] = []
    j = 0
    for pos in plot_positions:
        if np.isnan(pos):
            result.append(np.nan)
        else:
            result.append(raw_depths[j])
            j += 1
    return result


def _build_heatmap_matrix(
    region_df: pd.DataFrame,
    sample_cols: List[str],
    region_start: int,
    region_end: int,
    n_viz_bins: int = 1000,
) -> np.ndarray:
    """Build a dense (n_samples, n_viz_bins) matrix for heatmap display.

    The visualisation grid has a fixed number of pixels (*n_viz_bins*)
    spanning ``[region_start, region_end)``.  Each real bin is painted
    only at the pixel columns it truly covers, so gaps between bins
    (e.g. breakpoint regions) naturally remain NaN (white).

    Args:
        region_df: DataFrame with ``Start``, ``End`` and sample depth columns.
        sample_cols: Sample column names to include.
        region_start: Left edge of the visualisation region.
        region_end: Right edge of the visualisation region.
        n_viz_bins: Number of pixel columns in the output matrix.

    Returns:
        2-D array of shape ``(len(sample_cols), n_viz_bins)``.
    """
    n_samples = len(sample_cols)
    region_span = max(region_end - region_start, 1)
    viz_matrix = np.full((n_samples, n_viz_bins), np.nan)

    bin_starts = region_df["Start"].values
    bin_ends = region_df["End"].values
    data = region_df[sample_cols].values  # shape (n_bins, n_samples)

    for i, (bs, be) in enumerate(zip(bin_starts, bin_ends)):
        # Map genomic coords to pixel indices
        idx_lo = int((bs - region_start) / region_span * n_viz_bins)
        idx_hi = int((be - region_start) / region_span * n_viz_bins)
        idx_lo = max(idx_lo, 0)
        idx_hi = min(idx_hi, n_viz_bins)
        if idx_hi <= idx_lo:
            idx_hi = min(idx_lo + 1, n_viz_bins)
        viz_matrix[:, idx_lo:idx_hi] = data[i, :, np.newaxis]

    return viz_matrix


def _draw_overview_column(
    axes_col: List,
    region_df: pd.DataFrame,
    locus: GDLocus,
    calls_df: pd.DataFrame,
    region_start: int,
    region_end: int,
    carrier_cols: List[str],
    non_carrier_cols: List[str],
    all_ploidies: List[int],
    ploidy_lookup: Dict[Tuple[str, str], int],
    sample_cols: List[str],
    carriers,
    gtf: Optional[GTFParser],
    segdup: Optional[SegDupAnnotation],
    col_title: str,
    min_gene_label_spacing: float = 0.05,
    show_colorbar: bool = True,
):
    """Draw annotation + carrier heatmap + non-carrier heatmap + mean-depth
    panels into a list of axes (one column of the overview figure).

    Args:
        axes_col: List of axes, length = 3 + len(all_ploidies).
        region_df: DataFrame with ``Start``, ``End`` and per-sample depth columns.
        locus: GDLocus.
        calls_df: Full calls DataFrame.
        region_start / region_end: Genomic extent for x-axis.
        carrier_cols / non_carrier_cols: Sample lists (already ploidy-sorted).
        all_ploidies: Sorted list of distinct ploidy values.
        ploidy_lookup: (sample, contig) -> ploidy.
        sample_cols: All sample columns in *region_df*.
        carriers: Array/set of carrier sample IDs.
        gtf / segdup: Optional annotation objects.
        col_title: Title placed on the annotation panel.
        min_gene_label_spacing: Forwarded to draw_annotations_panel.
        show_colorbar: Whether to add a colorbar below the non-carrier heatmap.
    """
    chrom = locus.chrom
    n_carriers = len(carrier_cols)
    n_non_carriers = len(non_carrier_cols)
    n_ploidy_panels = len(all_ploidies)

    # Panel 1: Gene annotations and segdup regions
    ax = axes_col[0]
    draw_annotations_panel(
        ax, locus, region_start, region_end, chrom, col_title,
        gtf=gtf, segdup=segdup, show_gd_entries=True,
        min_gene_label_spacing=min_gene_label_spacing,
    )

    # Panel 2: Carriers heatmap
    ax = axes_col[1]
    if n_carriers > 0:
        viz_matrix = _build_heatmap_matrix(
            region_df, carrier_cols, region_start, region_end)
        cmap = plt.cm.RdBu_r.copy()
        cmap.set_bad(color='white')
        ax.imshow(viz_matrix, aspect="auto", cmap=cmap,
                  vmin=0, vmax=4, interpolation="nearest",
                  extent=[region_start, region_end, n_carriers, 0])
        ax.set_ylabel(f"Carriers (n={n_carriers})")
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_xlim(region_start, region_end)
        ax.set_ylim(0, n_carriers)
    else:
        ax.axis('off')
        ax.text(0.5, 0.5, 'No carriers', ha='center', va='center',
                transform=ax.transAxes)

    # Panel 3: Non-carriers heatmap
    ax = axes_col[2]
    if n_non_carriers > 0:
        viz_matrix = _build_heatmap_matrix(
            region_df, non_carrier_cols, region_start, region_end)
        cmap = plt.cm.RdBu_r.copy()
        cmap.set_bad(color='white')
        im = ax.imshow(viz_matrix, aspect="auto", cmap=cmap,
                       vmin=0, vmax=4, interpolation="nearest",
                       extent=[region_start, region_end, n_non_carriers, 0])
        ax.set_ylabel(f"Non-carriers (n={n_non_carriers})")
        ax.set_yticks([])
        ax.set_xlim(region_start, region_end)
        ax.set_ylim(0, n_non_carriers)
        ax.set_xlabel(f"Position on {chrom}")
        if show_colorbar:
            cbar = plt.colorbar(im, ax=ax, orientation='horizontal',
                                pad=0.08, aspect=40)
            cbar.set_label("Normalized read depth")
    else:
        ax.axis('off')
        ax.text(0.5, 0.5, 'No non-carriers', ha='center', va='center',
                transform=ax.transAxes)

    # ---- Mean depth panels: one per ploidy state ----
    bin_mids = (region_df["Start"].values + region_df["End"].values) / 2
    bin_starts = region_df["Start"].values
    bin_ends = region_df["End"].values
    plot_positions = _build_gap_positions(bin_mids, bin_starts, bin_ends)

    for panel_idx, ploidy_val in enumerate(all_ploidies):
        ax = axes_col[3 + panel_idx]

        ploidy_samples = [s for s in sample_cols
                          if ploidy_lookup.get((s, chrom), 2) == ploidy_val]
        ploidy_carriers = [s for s in ploidy_samples if s in carriers]
        ploidy_non_carriers = [s for s in ploidy_samples if s not in carriers]

        if len(ploidy_carriers) > 0 and len(ploidy_non_carriers) > 0:
            nc_depths = _depths_with_gaps(
                region_df[ploidy_non_carriers].mean(axis=1).values, plot_positions)
            ax.plot(plot_positions, nc_depths, "b-", linewidth=1, alpha=0.7,
                    label=f"Non-carriers (n={len(ploidy_non_carriers)})")

            carrier_calls = calls_df[
                (calls_df["cluster"] == locus.cluster) &
                (calls_df["is_carrier"]) &
                (calls_df["is_best_match"]) &
                (calls_df["sample"].isin(ploidy_carriers))
            ].copy()

            def get_bp_labels(row):
                for entry in locus.gd_entries:
                    if (entry["start_GRCh38"] == row["start"] and
                            entry["end_GRCh38"] == row["end"] and
                            entry["svtype"] == row["svtype"]):
                        return f"{entry['BP1']}-{entry['BP2']}"
                for entry in locus.gd_entries:
                    if entry["start_GRCh38"] == row["start"] and entry["end_GRCh38"] == row["end"]:
                        return f"{entry['BP1']}-{entry['BP2']}"
                return "unknown"

            if len(carrier_calls) > 0:
                carrier_calls["bp_interval"] = carrier_calls.apply(
                    get_bp_labels, axis=1)
                del_colors = ['#FF6B6B', '#FF4757', '#EE5A6F', '#C23B4E']
                dup_colors = ['#6B9BD1', '#4169E1', '#5080D0', '#3A5BB8']
                color_idx = 0
                plotted_any = False
                for (svtype, bp_interval), group in carrier_calls.groupby(
                        ["svtype", "bp_interval"]):
                    valid = [s for s in group["sample"].unique()
                             if s in ploidy_carriers]
                    if len(valid) > 0:
                        gd = _depths_with_gaps(
                            region_df[valid].mean(axis=1).values,
                            plot_positions)
                        if not all(np.isnan(gd)):
                            color = (del_colors if svtype == "DEL"
                                     else dup_colors)[color_idx % 4]
                            ax.plot(plot_positions, gd, "-", linewidth=1.5,
                                    alpha=0.8, color=color,
                                    label=f"{svtype} {bp_interval} (n={len(valid)})")
                            color_idx += 1
                            plotted_any = True
                if plotted_any:
                    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.35),
                              ncol=6, fontsize=7, frameon=True, fancybox=True)
        else:
            if len(ploidy_samples) > 0:
                md = _depths_with_gaps(
                    region_df[ploidy_samples].mean(axis=1).values,
                    plot_positions)
                ax.plot(plot_positions, md, "b-", linewidth=1.5,
                        label=f"All (n={len(ploidy_samples)})")
                ax.legend(fontsize=7)

        for bp_start, bp_end in locus.breakpoints:
            ax.axvline(bp_start, color="red", linestyle="-", alpha=1, zorder=5)
            ax.axvline(bp_end, color="red", linestyle="-", alpha=1, zorder=5)
            ax.axvspan(bp_start, bp_end, alpha=0.1, color="red", zorder=0)

        ax.axhline(float(ploidy_val), color="gray", linestyle=":", alpha=0.5)
        ax.set_xlim(region_start, region_end)
        ax.set_ylim(0, max(4, ploidy_val + 2))
        ax.set_ylabel(f"Depth (P={ploidy_val})")
        ax.grid(True, alpha=0.3, axis="y")

        if panel_idx == n_ploidy_panels - 1:
            ax.set_xlabel(f"Position on {chrom}")
        else:
            ax.set_xticks([])


def plot_locus_overview(
    locus: GDLocus,
    calls_df: pd.DataFrame,
    depth_df: pd.DataFrame,
    gtf: Optional[GTFParser],
    segdup: Optional[SegDupAnnotation],
    output_dir: str,
    padding: int = 50000,
    ploidy_df: Optional[pd.DataFrame] = None,
    min_gene_label_spacing: float = 0.05,
    raw_counts_df: Optional[pd.DataFrame] = None,
    raw_sample_medians: Optional[Dict[str, float]] = None,
):
    """
    Create overview plot for a GD locus showing all samples.

    When *raw_counts_df* is provided the figure has two columns:
      - **Left**: normalised raw depth (before filtering/rebinning)
      - **Right**: processed rebinned counts from the model

    Without raw data the figure is a single column (backward-compatible).

    Args:
        locus: GDLocus object
        calls_df: DataFrame with CNV calls for this locus
        depth_df: DataFrame with depth data (from cn_posteriors)
        gtf: Optional GTFParser for gene annotations
        segdup: Optional SegDupAnnotation for segdup regions
        output_dir: Directory to save plots
        padding: Padding around locus boundaries
        ploidy_df: Optional ploidy estimates (sample, contig, ploidy).
            Used to sort heatmaps and split mean-depth panel by ploidy.
        min_gene_label_spacing: Forwarded to draw_annotations_panel.
        raw_counts_df: Optional DataFrame with raw (un-normalised) read
            counts (Chr, Start, End + sample columns).  When provided the
            plot gains a left-hand column showing the raw normalised depth.
        raw_sample_medians: Dict mapping sample ID to its autosomal median
            raw count.  Required when *raw_counts_df* is given.
    """
    # Skip loci with no breakpoints
    if not locus.breakpoints:
        print(f"  Warning: No breakpoints defined for locus {locus.cluster}, skipping")
        return

    chrom = locus.chrom
    # Filter strictly to this locus's bins (by cluster name).
    mask = (
        (depth_df["Cluster"] == locus.cluster) &
        (depth_df["Chr"] == chrom)
    )
    region_df = depth_df[mask].copy().sort_values("Start")

    if len(region_df) == 0:
        print(f"  Warning: No depth data found for locus {locus.cluster} "
              f"({chrom}:{locus.start:,}-{locus.end:,}), skipping plot")
        return

    # Bounds derived from actual data so flanks are always visible
    region_start = int(region_df["Start"].min())
    region_end = int(region_df["End"].max())

    sample_cols = get_sample_columns(region_df)

    # ---- ploidy helpers ----
    ploidy_lookup = _build_ploidy_lookup(ploidy_df)

    # Get carriers for this locus
    carriers = calls_df[
        (calls_df["cluster"] == locus.cluster) &
        (calls_df["is_carrier"])
    ]["sample"].unique()

    # Separate carriers and non-carriers, sorted by ploidy
    carrier_cols = _sort_samples_by_ploidy(
        [s for s in sample_cols if s in carriers], chrom, ploidy_lookup)
    non_carrier_cols = _sort_samples_by_ploidy(
        [s for s in sample_cols if s not in carriers], chrom, ploidy_lookup)
    n_carriers = len(carrier_cols)
    n_non_carriers = len(non_carrier_cols)

    # Determine distinct ploidies present on this chromosome
    all_ploidies = sorted(set(
        ploidy_lookup.get((s, chrom), 2) for s in sample_cols
    ))
    n_ploidy_panels = len(all_ploidies)

    # ---- Build raw-normalised DataFrame for the left column (if available) ----
    have_raw = (raw_counts_df is not None and raw_sample_medians is not None)
    raw_region_df: Optional[pd.DataFrame] = None
    if have_raw:
        raw_mask = (
            (raw_counts_df["Chr"] == chrom)
            & (raw_counts_df["End"] > region_start)
            & (raw_counts_df["Start"] < region_end)
        )
        raw_region = raw_counts_df[raw_mask].copy().sort_values("Start")
        if len(raw_region) > 0:
            # Normalise each sample: 2 * raw / autosomal_median → CN-2 ≈ 2.0
            raw_norm_data = {}
            usable_samples = []
            for s in sample_cols:
                if s in raw_region.columns and s in raw_sample_medians:
                    med = raw_sample_medians[s]
                    if med > 0:
                        raw_norm_data[s] = 2.0 * raw_region[s].values.astype(float) / med
                        usable_samples.append(s)
            if usable_samples:
                raw_region_df = pd.concat(
                    [raw_region[["Chr", "Start", "End"]].reset_index(drop=True),
                     pd.DataFrame(raw_norm_data, columns=usable_samples)],
                    axis=1,
                )
        if raw_region_df is None or len(raw_region_df) == 0:
            have_raw = False
            raw_region_df = None

    # ---- figure layout ----
    n_cols = 2 if have_raw else 1
    carrier_height = min(6, max(1, n_carriers * 0.15)) if n_carriers > 0 else 0.5
    non_carrier_height = min(6, max(1, n_non_carriers * 0.15)) if n_non_carriers > 0 else 0.5
    depth_panel_height = 1.0
    fig_height = 3 + carrier_height + non_carrier_height + n_ploidy_panels * depth_panel_height

    height_ratios = (
        [1, carrier_height, non_carrier_height]
        + [depth_panel_height] * n_ploidy_panels
    )
    n_rows = 3 + n_ploidy_panels
    fig_width = 16 * n_cols

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(fig_width, fig_height),
        gridspec_kw={"height_ratios": height_ratios},
        squeeze=False,
    )

    locus_label = f"{locus.cluster} ({chrom}:{locus.start:,}-{locus.end:,})"

    if have_raw:
        # Left column: raw normalised depth
        _draw_overview_column(
            [axes[r, 0] for r in range(n_rows)],
            raw_region_df, locus, calls_df,
            region_start, region_end,
            # Use the same carrier/non-carrier split but only samples present
            # in the raw data
            [s for s in carrier_cols if s in raw_region_df.columns],
            [s for s in non_carrier_cols if s in raw_region_df.columns],
            all_ploidies, ploidy_lookup,
            [s for s in sample_cols if s in raw_region_df.columns],
            carriers, gtf, segdup,
            col_title=f"Raw normalised — {locus_label}",
            min_gene_label_spacing=min_gene_label_spacing,
            show_colorbar=False,
        )
        # Right column: processed rebinned counts
        _draw_overview_column(
            [axes[r, 1] for r in range(n_rows)],
            region_df, locus, calls_df,
            region_start, region_end,
            carrier_cols, non_carrier_cols,
            all_ploidies, ploidy_lookup,
            sample_cols, carriers, gtf, segdup,
            col_title=f"Processed — {locus_label}",
            min_gene_label_spacing=min_gene_label_spacing,
            show_colorbar=True,
        )
    else:
        # Single column (backward-compatible)
        _draw_overview_column(
            [axes[r, 0] for r in range(n_rows)],
            region_df, locus, calls_df,
            region_start, region_end,
            carrier_cols, non_carrier_cols,
            all_ploidies, ploidy_lookup,
            sample_cols, carriers, gtf, segdup,
            col_title=locus_label,
            min_gene_label_spacing=min_gene_label_spacing,
            show_colorbar=True,
        )

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
    min_gene_label_spacing: float = 0.05,
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
        min_gene_label_spacing: Forwarded to draw_annotations_panel.
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
    draw_annotations_panel(axes[0], locus, region_start, region_end, chrom, title, gtf, segdup,
                           show_gd_entries=True, min_gene_label_spacing=min_gene_label_spacing)
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
    min_gene_label_spacing: float = 0.05,
    raw_counts_df: Optional[pd.DataFrame] = None,
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
        min_gene_label_spacing: Forwarded to draw_annotations_panel.
        raw_counts_df: Optional DataFrame with raw (un-normalised) read
            counts.  Columns: ``Chr``, ``Start``, ``End``, and one per
            sample.  When provided, the raw depth is normalised per-sample
            and overlaid on each carrier depth panel.
    """
    carriers = calls_df[calls_df["is_carrier"]]

    if len(carriers) == 0:
        print("No carriers to include in PDF.")
        return

    # Precompute per-sample genome-wide autosomal median raw counts so we
    # can normalise the raw depth trace to a CN-2 baseline.
    raw_sample_medians: Dict[str, float] = {}
    if raw_counts_df is not None:
        autosomal = [str(c) for c in range(1, 23)] + [f"chr{c}" for c in range(1, 23)]
        raw_auto = raw_counts_df[raw_counts_df["Chr"].isin(autosomal)]
        raw_sample_cols = [c for c in raw_counts_df.columns
                          if c not in ("Chr", "Start", "End", "source_file", "Bin")]
        for s in raw_sample_cols:
            med = raw_auto[s].median()
            if med > 0:
                raw_sample_medians[s] = med
        print(f"  Computed raw autosomal medians for {len(raw_sample_medians)} samples")

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
                draw_annotations_panel(axes[0], locus, region_start, region_end, chrom, title, gtf, segdup,
                                       show_gd_entries=False, min_gene_label_spacing=min_gene_label_spacing)

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

                # Overlay raw depth trace if available
                if raw_counts_df is not None and sample_id in raw_sample_medians:
                    raw_region = raw_counts_df[
                        (raw_counts_df["Chr"] == chrom)
                        & (raw_counts_df["End"] > region_start)
                        & (raw_counts_df["Start"] < region_end)
                    ].sort_values("Start")
                    if len(raw_region) > 0 and sample_id in raw_region.columns:
                        raw_mids = (raw_region["Start"].values + raw_region["End"].values) / 2
                        raw_depth = raw_region[sample_id].values.astype(float)
                        # Normalise: depth / genome_median * 2  →  CN-2 ≈ 2.0
                        norm_raw = 2.0 * raw_depth / raw_sample_medians[sample_id]
                        ax.plot(raw_mids, norm_raw, color="darkorange", linewidth=0.6,
                                alpha=0.7, zorder=4, label="Raw depth")

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
# Truth Table Evaluation
# =============================================================================


def load_truth_table(filepath: str) -> pd.DataFrame:
    """
    Load a truth table of manually labelled GD carrier calls.

    Expected columns:
        chr, start, end, GD_ID, cluster_ID, SVTYPE, carriers, non_carriers

    The *carriers* column is a comma-separated list of sample IDs.
    *cluster_ID* may be blank.  *non_carriers* is ignored.

    Args:
        filepath: Path to TSV truth table.

    Returns:
        DataFrame with one row per GD site, carriers parsed into a set.
    """
    df = pd.read_csv(filepath, sep="\t", dtype=str)
    # Normalise column names to lowercase for robustness
    df.columns = [c.strip() for c in df.columns]
    required = {"GD_ID", "carriers"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Truth table missing required columns: {missing}")

    # Parse carriers into sets; treat empty / NaN as no carriers
    def _parse_carriers(val):
        if pd.isna(val) or str(val).strip() == "":
            return set()
        return set(s.strip() for s in str(val).split(",") if s.strip())

    df["carrier_set"] = df["carriers"].apply(_parse_carriers)
    return df


def evaluate_against_truth(
    calls_df: pd.DataFrame,
    truth_df: pd.DataFrame,
    output_dir: str,
    batch_samples: Optional[set] = None,
) -> pd.DataFrame:
    """
    Cross-reference predicted GD calls against a truth table and report
    sensitivity and precision for every site with at least one truth or
    predicted carrier.

    Matching is performed by *GD_ID*.  For each GD_ID the set of predicted
    carriers (samples where ``is_carrier == True``) is compared to the truth
    carrier set.

    Truth carrier sets are intersected with *batch_samples* (when provided)
    so that sensitivity is not penalised for labelled carriers absent from
    the current batch.

    Args:
        calls_df: Predicted calls DataFrame (from call_cnvs_from_posteriors
            or loaded from ``--calls``).
        truth_df: Truth table DataFrame produced by :func:`load_truth_table`.
        output_dir: Directory to write the report TSV.
        batch_samples: Optional set of sample IDs present in the current
            batch.  If provided, truth carriers not in this set are removed
            before scoring.

    Returns:
        Per-site report DataFrame.
    """
    print("\n" + "=" * 80)
    print("EVALUATING PREDICTIONS AGAINST TRUTH TABLE")
    print("=" * 80)

    if batch_samples is not None:
        print(f"  Batch contains {len(batch_samples)} samples; "
              "truth carriers will be restricted to this set.")

    # Build predicted carrier sets keyed by GD_ID.
    # Only count samples whose call is both a carrier AND the best-match
    # breakpoint configuration for its svtype.  Without the is_best_match
    # filter, a sample carrying a BP1-3 DUP would also appear as a
    # (spurious) carrier of BP1-2 and BP2-3, inflating false-positive
    # counts for those sub-intervals.
    pred_by_gd: Dict[str, set] = {}
    pred_meta_by_gd: Dict[str, dict] = {}  # fallback metadata from calls
    for gd_id, grp in calls_df.groupby("GD_ID"):
        gd_id_str = str(gd_id)
        carrier_mask = grp["is_carrier"] == True  # noqa: E712
        if "is_best_match" in grp.columns:
            carrier_mask = carrier_mask & (grp["is_best_match"] == True)  # noqa: E712
        pred_by_gd[gd_id_str] = set(
            grp.loc[carrier_mask, "sample"].unique()
        )
        # Grab metadata from the first row of this GD_ID group so we can
        # fill in chr/start/end/cluster/svtype when the truth table has no
        # entry for this GD_ID.
        first = grp.iloc[0]
        pred_meta_by_gd[gd_id_str] = {
            "chr": first.get("chrom", ""),
            "start": first.get("start", ""),
            "end": first.get("end", ""),
            "cluster_ID": first.get("cluster", ""),
            "SVTYPE": first.get("svtype", ""),
        }

    # Build truth carrier sets keyed by GD_ID
    truth_by_gd: Dict[str, dict] = {}
    for _, row in truth_df.iterrows():
        gd_id = str(row["GD_ID"])
        carrier_set = row["carrier_set"]
        # Restrict truth carriers to samples in the current batch
        if batch_samples is not None:
            carrier_set = carrier_set & batch_samples
        truth_by_gd[gd_id] = {
            "carrier_set": carrier_set,
            "chr": row.get("chr", ""),
            "start": row.get("start", ""),
            "end": row.get("end", ""),
            "cluster_ID": row.get("cluster_ID", ""),
            "SVTYPE": row.get("SVTYPE", ""),
        }

    # Union of all GD_IDs with at least one truth or predicted carrier
    all_gd_ids = sorted(
        {gd for gd, s in pred_by_gd.items() if len(s) > 0}
        | {gd for gd, d in truth_by_gd.items() if len(d["carrier_set"]) > 0}
    )

    rows = []
    total_tp = total_fp = total_fn = 0
    for gd_id in all_gd_ids:
        truth_set = truth_by_gd.get(gd_id, {}).get("carrier_set", set())
        pred_set = pred_by_gd.get(gd_id, set())

        tp = len(truth_set & pred_set)
        fp = len(pred_set - truth_set)
        fn = len(truth_set - pred_set)
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
        precision = tp / (tp + fp) if (tp + fp) > 0 else float("nan")

        total_tp += tp
        total_fp += fp
        total_fn += fn

        meta = truth_by_gd.get(gd_id, {})
        # Fall back to prediction metadata when the truth table has no
        # entry for this GD_ID (e.g. a novel/unannotated site).
        fallback = pred_meta_by_gd.get(gd_id, {})
        rows.append({
            "GD_ID": gd_id,
            "chr": meta.get("chr", "") or fallback.get("chr", ""),
            "start": meta.get("start", "") or fallback.get("start", ""),
            "end": meta.get("end", "") or fallback.get("end", ""),
            "cluster_ID": meta.get("cluster_ID", "") or fallback.get("cluster_ID", ""),
            "SVTYPE": meta.get("SVTYPE", "") or fallback.get("SVTYPE", ""),
            "n_truth_carriers": len(truth_set),
            "n_pred_carriers": len(pred_set),
            "TP": tp,
            "FP": fp,
            "FN": fn,
            "sensitivity": round(sensitivity, 4) if not np.isnan(sensitivity) else "NA",
            "precision": round(precision, 4) if not np.isnan(precision) else "NA",
            "FP_samples": ",".join(sorted(pred_set - truth_set)) if fp > 0 else "",
            "FN_samples": ",".join(sorted(truth_set - pred_set)) if fn > 0 else "",
        })

    report_df = pd.DataFrame(rows)

    # Print summary
    overall_sens = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else float("nan")
    overall_prec = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else float("nan")
    print(f"\n  Sites evaluated: {len(all_gd_ids)}")
    print(f"  Overall TP={total_tp}  FP={total_fp}  FN={total_fn}")
    print(f"  Overall sensitivity: {overall_sens:.4f}")
    print(f"  Overall precision:   {overall_prec:.4f}")

    # Per-site summary to stdout
    print(f"\n  {'GD_ID':30s}  {'SVTYPE':6s}  truth  pred    TP    FP    FN   sens   prec")
    print("  " + "-" * 100)
    for r in rows:
        sens_str = f"{r['sensitivity']:.2f}" if r["sensitivity"] != "NA" else "  NA"
        prec_str = f"{r['precision']:.2f}" if r["precision"] != "NA" else "  NA"
        print(f"  {r['GD_ID']:30s}  {str(r['SVTYPE']):6s}  "
              f"{r['n_truth_carriers']:5d}  {r['n_pred_carriers']:4d}  "
              f"{r['TP']:4d}  {r['FP']:4d}  {r['FN']:4d}  "
              f"{sens_str:>5s}  {prec_str:>5s}")

    # Write report
    output_path = os.path.join(output_dir, "truth_evaluation_report.tsv")
    report_df.to_csv(output_path, sep="\t", index=False)
    print(f"\n  Saved report: {output_path}")
    print("=" * 80 + "\n")

    return report_df


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
        "--ploidy-table",
        required=False,
        help="Ploidy estimates TSV (ploidy_estimates.tsv) from gd_cnv_pyro.py. "
             "Columns: sample, contig, median_depth, ploidy. "
             "If not provided, ploidy=2 is assumed for all sample/contig pairs.",
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
        "--raw-counts",
        required=False,
        help="Raw read-count matrix TSV (can be gzipped).  Same format as "
             "the input to gd_cnv_pyro.py (#Chr / Start / End + sample columns "
             "with un-normalised integer counts).  When provided, the raw depth "
             "profile is superimposed on each carrier plot, normalised by the "
             "sample's genome-wide autosomal median and scaled to a CN-2 baseline.",
    )
    parser.add_argument(
        "--output-dir", "-o",
        required=True,
        help="Output directory for plots",
    )
    parser.add_argument(
        "--log-prob-threshold",
        type=float,
        default=-0.3,
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
    parser.add_argument(
        "--truth-table",
        required=False,
        help="Optional truth table TSV with manually labelled GD carriers. "
             "Columns: chr, start, end, GD_ID, cluster_ID, SVTYPE, carriers, non_carriers. "
             "If provided, a sensitivity/precision report is produced.",
    )
    parser.add_argument(
        "--sample-posteriors",
        required=False,
        help="Sample posteriors table (TSV) with columns: sample, sample_var_map. "
             "Used to identify the set of samples in the current batch. "
             "When combined with --truth-table, truth carriers absent from "
             "this sample set are excluded from sensitivity scoring.",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print detailed per-sample log probability scores for all GD "
             "entries at every locus, including flanking region scores.",
    )
    parser.add_argument(
        "--min-gene-label-spacing",
        type=float,
        default=0.05,
        help="Minimum distance between adjacent gene label centres as a fraction "
             "of the plot width.  Gene lines are always drawn; labels are "
             "suppressed when centres are closer than this threshold (default: 0.05 = 5%%).",
    )

    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Mirror all stdout/stderr to a log file in the output directory
    log_fh = setup_logging(args.output_dir)

    print(f"Output directory: {args.output_dir}")

    # Load data
    print("\nLoading data...")

    print(f"  Loading CN posteriors: {args.cn_posteriors}")
    cn_posteriors_df = pd.read_csv(args.cn_posteriors, sep="\t", compression="infer")
    print(f"    {len(cn_posteriors_df)} bin-sample records")
    
    print(f"  Loading GD table: {args.gd_table}")
    gd_table = GDTable(args.gd_table)
    print(f"    {len(gd_table.loci)} loci")
    
    # Load ploidy table if provided
    ploidy_df = None
    if args.ploidy_table:
        print(f"  Loading ploidy table: {args.ploidy_table}")
        ploidy_df = pd.read_csv(args.ploidy_table, sep="\t")
        print(f"    {len(ploidy_df)} sample/contig ploidy records")
    
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
            ploidy_df=ploidy_df,
            verbose=args.verbose,
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

    # Load raw counts matrix if provided
    raw_counts_df = None
    raw_sample_medians: Dict[str, float] = {}
    if args.raw_counts:
        print(f"\n  Loading raw counts: {args.raw_counts}")
        raw_counts_df = pd.read_csv(args.raw_counts, sep="\t", compression="infer")
        if "#Chr" in raw_counts_df.columns:
            raw_counts_df.rename(columns={"#Chr": "Chr"}, inplace=True)
        raw_counts_df["Start"] = raw_counts_df["Start"].astype(int)
        raw_counts_df["End"] = raw_counts_df["End"].astype(int)
        print(f"    {len(raw_counts_df)} bins")
        # Precompute per-sample autosomal median raw counts for normalisation
        autosomal = [str(c) for c in range(1, 23)] + [f"chr{c}" for c in range(1, 23)]
        raw_auto = raw_counts_df[raw_counts_df["Chr"].isin(autosomal)]
        raw_sample_cols = [c for c in raw_counts_df.columns
                          if c not in ("Chr", "Start", "End", "source_file", "Bin")]
        for s in raw_sample_cols:
            med = raw_auto[s].median()
            if med > 0:
                raw_sample_medians[s] = med
        print(f"    Computed raw autosomal medians for {len(raw_sample_medians)} samples")

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
            args.output_dir, padding=args.padding,
            ploidy_df=ploidy_df,
            min_gene_label_spacing=args.min_gene_label_spacing,
            raw_counts_df=raw_counts_df,
            raw_sample_medians=raw_sample_medians if raw_sample_medians else None,
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
                    args.output_dir, padding=args.padding,
                    min_gene_label_spacing=args.min_gene_label_spacing,
                )

    # Load sample posteriors to identify batch samples (if provided)
    batch_samples = None
    if args.sample_posteriors:
        print(f"\n  Loading sample posteriors: {args.sample_posteriors}")
        sample_post_df = pd.read_csv(args.sample_posteriors, sep="\t")
        batch_samples = set(sample_post_df["sample"].astype(str).unique())
        print(f"    {len(batch_samples)} samples in batch")

    # Evaluate against truth table if provided
    if args.truth_table:
        print("\nLoading truth table...")
        truth_df = load_truth_table(args.truth_table)
        print(f"  {len(truth_df)} truth entries from {args.truth_table}")
        evaluate_against_truth(calls_df, truth_df, args.output_dir,
                               batch_samples=batch_samples)

    # Create carrier PDF
    print("\nCreating carrier PDF...")
    create_carrier_pdf(
        calls_df, depth_df, gd_table, gtf, segdup,
        args.output_dir, padding=args.padding,
        min_gene_label_spacing=args.min_gene_label_spacing,
        raw_counts_df=raw_counts_df,
    )

    print("\n" + "=" * 80)
    print("Plotting complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
