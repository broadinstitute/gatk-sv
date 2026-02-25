#!/bin/python

"""
Genomic Disorder (GD) CNV Detection from Binned Read Counts

This script detects copy number variants at known genomic disorder loci using
a hierarchical Bayesian model. Each locus (cluster) is defined by one or more
breakpoints, and the script determines which breakpoint set best explains the
observed depth signal for each sample.

Input:
    - Binned read count file (TSV)
    - GD table with locus definitions (TSV)
    - Segmental duplication regions (BED.gz) to mask

Output:
    - Per-locus CNV calls with breakpoint assignments
    - Copy number posteriors for each sample at each locus
"""

import argparse
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import torch

from pathlib import Path
from tqdm import tqdm

import pyro
import pyro.distributions as dist
from pyro import poutine
from pyro.ops.indexing import Vindex
from pyro.infer import config_enumerate, infer_discrete
from pyro.infer.autoguide import AutoDiagonalNormal, AutoDelta
from pyro.infer import JitTraceEnum_ELBO, TraceEnum_ELBO
from pyro.infer.svi import SVI


def get_sample_columns(df: pd.DataFrame) -> list:
    """
    Extract sample column names from a DataFrame by excluding metadata columns.

    Args:
        df: DataFrame with bins as rows and samples as columns

    Returns:
        List of sample column names
    """
    metadata_cols = ["Chr", "Start", "End", "source_file"]
    sample_cols = [col for col in df.columns if col not in metadata_cols]
    return sample_cols


@dataclass
class GDLocus:
    """Represents a genomic disorder locus with its breakpoints."""
    cluster: str
    chrom: str
    breakpoints: List[Tuple[int, int]]  # Sorted list of breakpoint ranges (start, end)
    breakpoint_names: List[str]  # Names of breakpoints (e.g., ['1', '2', '3'] or ['A', 'B', 'C'])
    gd_entries: List[dict]  # Original GD table entries for this locus
    is_nahr: bool
    is_terminal: bool

    @property
    def n_breakpoints(self) -> int:
        return len(self.breakpoints)

    @property
    def start(self) -> int:
        """Get overall start position (min of all breakpoint starts)."""
        return min(bp[0] for bp in self.breakpoints)

    @property
    def end(self) -> int:
        """Get overall end position (max of all breakpoint ends)."""
        return max(bp[1] for bp in self.breakpoints)

    def get_intervals(self) -> List[Tuple[int, int, str]]:
        """
        Get all intervals between adjacent breakpoints.
        Intervals span from the end of one breakpoint region to the start of the next.

        Returns:
            List of (start, end, interval_name) tuples
        """
        intervals = []
        for i in range(len(self.breakpoints) - 1):
            # Interval starts at end of BP_i and ends at start of BP_{i+1}
            start = self.breakpoints[i][1]  # End of current breakpoint
            end = self.breakpoints[i + 1][0]  # Start of next breakpoint
            # Use actual breakpoint names from the table
            interval_name = f"{self.breakpoint_names[i]}-{self.breakpoint_names[i + 1]}"
            intervals.append((start, end, interval_name))
        return intervals

    def get_gd_intervals(self) -> Dict[str, Tuple[int, int]]:
        """
        Get intervals for each GD entry (each DEL/DUP defined in the table).

        Returns:
            Dict mapping GD_ID to (start, end) tuple
        """
        return {
            entry["GD_ID"]: (entry["start_GRCh38"], entry["end_GRCh38"])
            for entry in self.gd_entries
        }

    def get_flanking_regions(self) -> List[Tuple[int, int, str]]:
        """
        Get flanking regions on either side of the locus (100% of locus size each side).

        These regions are used to detect large spanning CNVs that overlap the locus
        incidentally but are not true GD events. A true GD only affects the region
        between its breakpoints; a spanning variant would also show copy number change
        in these flanking regions.

        Note: The left flank is clipped to position 0 if the locus is near the
        chromosome start. Flanks at the chromosome end will simply have no bins.

        Returns:
            List of (start, end, name) tuples. Contains 0, 1, or 2 elements
            depending on chromosome proximity.
        """
        locus_size = self.end - self.start
        flanks = []

        # Left flank: extend left by 100% of locus size
        left_start = max(0, self.start - locus_size)
        left_end = self.start
        if left_end > left_start:
            flanks.append((left_start, left_end, "left_flank"))

        # Right flank: extend right by 100% of locus size
        right_start = self.end
        right_end = self.end + locus_size
        flanks.append((right_start, right_end, "right_flank"))

        return flanks

    @property
    def del_entries(self) -> List[dict]:
        """Get all deletion entries for this locus."""
        return [e for e in self.gd_entries if e["svtype"] == "DEL"]

    @property
    def dup_entries(self) -> List[dict]:
        """Get all duplication entries for this locus."""
        return [e for e in self.gd_entries if e["svtype"] == "DUP"]

    @property
    def svtypes(self) -> List[str]:
        """Get unique svtypes at this locus."""
        return list(set(e["svtype"] for e in self.gd_entries))


class GDTable:
    """Parser and container for GD locus definitions."""

    def __init__(self, filepath: str):
        """
        Load GD table from TSV file.

        Expected columns:
            chr, start_GRCh38, end_GRCh38, GD_ID, svtype, NAHR, terminal, cluster

        Args:
            filepath: Path to GD table TSV file
        """
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
        """
        Parse GD table into loci grouped by cluster.
        Breakpoints are extracted from BP1 and BP2 columns.

        Returns:
            Dict mapping cluster name to GDLocus object
        """
        loci = {}

        # Group by cluster, only processing entries with a defined cluster
        for _, row in self.df.iterrows():
            cluster = row["cluster"]
            
            # For entries without a cluster (standalone regions), use GD_ID as cluster name
            if pd.isna(cluster) or cluster == "":
                cluster = row["GD_ID"]

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
            
            # If BP1 or BP2 is missing, use default values "1" and "2"
            if pd.isna(bp1) or pd.isna(bp2):
                bp1 = "1"
                bp2 = "2"
            else:
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

            # Add entry
            loci[cluster]["entries"].append({
                "GD_ID": row["GD_ID"],
                "start_GRCh38": row["start_GRCh38"],
                "end_GRCh38": row["end_GRCh38"],
                "svtype": row["svtype"],
                "NAHR": row["NAHR"],
                "terminal": row["terminal"],
                "BP1": bp1,
                "BP2": bp2,
            })

        # Convert to GDLocus objects
        result = {}
        for cluster, data in loci.items():
            # Convert breakpoint coordinates to ranges
            bp_ranges = []
            
            if data["breakpoint_coords"]:
                # Sort breakpoint names (numeric if possible, otherwise alphabetic)
                def sort_key(bp_name):
                    # Try to convert to int for numeric breakpoints
                    try:
                        return (0, int(bp_name))
                    except ValueError:
                        # Alphabetic breakpoints
                        return (1, bp_name)
                
                bp_names_sorted = sorted(data["breakpoint_coords"].keys(), key=sort_key)
                
                for bp_name in bp_names_sorted:
                    coords = data["breakpoint_coords"][bp_name]
                    # Breakpoint range spans from min to max observed coordinate
                    bp_range = (min(coords), max(coords))
                    bp_ranges.append(bp_range)
            else:
                # No breakpoints found at all - this shouldn't happen, but handle it
                print(f"  ERROR: No breakpoints found for {cluster}, skipping")
                continue
            
            result[cluster] = GDLocus(
                cluster=cluster,
                chrom=data["chrom"],
                breakpoints=bp_ranges,
                breakpoint_names=bp_names_sorted,
                gd_entries=data["entries"],
                is_nahr=data["is_nahr"],
                is_terminal=data["is_terminal"],
            )
            
            # Print breakpoint ranges for debugging
            print(f"  {cluster}: {len(bp_ranges)} breakpoints")
            for i, (start, end) in enumerate(bp_ranges):
                width = end - start
                print(f"    BP{i+1}: {start:,} - {end:,} (width: {width:,} bp)")

        return result

    def get_locus(self, cluster: str) -> Optional[GDLocus]:
        """Get a specific locus by cluster name."""
        return self.loci.get(cluster)

    def get_all_loci(self) -> Dict[str, GDLocus]:
        """Get all loci."""
        return self.loci

    def get_loci_by_chrom(self, chrom: str) -> Dict[str, GDLocus]:
        """Get all loci on a specific chromosome."""
        return {k: v for k, v in self.loci.items() if v.chrom == chrom}


class SegDupMask:
    """
    Handles segmental duplication regions for masking unreliable depth signal.

    Segmental duplications have variable copy number in the reference and
    unreliable depth signal, so we mask them from CNV analysis.
    """

    def __init__(self, filepath: str):
        """
        Load segmental duplication regions from BED file.

        Expected format (BED6):
            chr, start, end, name, score, strand

        Args:
            filepath: Path to compressed BED file (.bed.gz)
        """
        self.df = pd.read_csv(
            filepath,
            sep="\t",
            header=None,
            names=["chr", "start", "end", "name", "score", "strand"],
            compression="gzip" if filepath.endswith(".gz") else None,
        )
        self._build_interval_index()
        print(f"Loaded {len(self.df)} segmental duplication regions")

    def _build_interval_index(self):
        """Build efficient interval index for overlap queries."""
        # Group by chromosome for efficient lookup
        self.intervals_by_chrom = {}
        for chrom, group in self.df.groupby("chr"):
            # Store as numpy arrays for fast operations
            self.intervals_by_chrom[chrom] = {
                "starts": group["start"].values,
                "ends": group["end"].values,
            }

    def get_overlap_fraction(self, chrom: str, start: int, end: int) -> float:
        """
        Calculate the fraction of a region that overlaps with segmental duplications.

        Args:
            chrom: Chromosome name
            start: Region start position
            end: Region end position

        Returns:
            Fraction of region overlapping with segdups (0.0 to 1.0)
        """
        if chrom not in self.intervals_by_chrom:
            return 0.0

        intervals = self.intervals_by_chrom[chrom]
        region_length = end - start

        if region_length <= 0:
            return 0.0

        # Find overlapping intervals
        sd_starts = intervals["starts"]
        sd_ends = intervals["ends"]

        # Check for any overlap: sd_start < end AND sd_end > start
        overlaps = (sd_starts < end) & (sd_ends > start)

        if not overlaps.any():
            return 0.0

        # Calculate total overlap length (handling overlapping segdups)
        overlap_starts = np.maximum(sd_starts[overlaps], start)
        overlap_ends = np.minimum(sd_ends[overlaps], end)

        # Merge overlapping intervals to avoid double-counting
        # Sort by start position
        sorted_idx = np.argsort(overlap_starts)
        overlap_starts = overlap_starts[sorted_idx]
        overlap_ends = overlap_ends[sorted_idx]

        # Merge overlapping intervals
        merged_length = 0
        current_start = overlap_starts[0]
        current_end = overlap_ends[0]

        for i in range(1, len(overlap_starts)):
            if overlap_starts[i] <= current_end:
                # Overlapping or adjacent, extend current interval
                current_end = max(current_end, overlap_ends[i])
            else:
                # Non-overlapping, add current interval and start new one
                merged_length += current_end - current_start
                current_start = overlap_starts[i]
                current_end = overlap_ends[i]

        # Add final interval
        merged_length += current_end - current_start

        return min(merged_length / region_length, 1.0)

    def is_masked(self, chrom: str, start: int, end: int,
                  threshold: float = 0.5) -> bool:
        """
        Check if a region should be masked due to segdup overlap.

        Args:
            chrom: Chromosome name
            start: Region start position
            end: Region end position
            threshold: Minimum overlap fraction to trigger masking

        Returns:
            True if region should be masked
        """
        return self.get_overlap_fraction(chrom, start, end) >= threshold


class DepthData:
    """Data container for CNV detection model"""

    def __init__(
        self,
        df: pd.DataFrame,
        device: str = "cpu",
        dtype: torch.dtype = torch.float32,
        subsample_bins: int = None,
        subsample_samples: int = None,
        seed: int = 42,
        clamp_threshold: float = 5.0,
    ):
        """
        Args:
            df: DataFrame with bins as rows and samples as columns
            device: torch device
            dtype: torch data type
            subsample_bins: If specified, randomly subsample this many bins
            subsample_samples: If specified, randomly subsample this many samples
            seed: Random seed for subsampling
            clamp_threshold: Maximum value for depth; values above this are clamped
        """
        # Store original dataframe
        self.original_df = df

        # Get sample columns (exclude metadata)
        sample_cols = get_sample_columns(df)

        # Subsample if requested
        if subsample_bins is not None or subsample_samples is not None:
            np.random.seed(seed)

            # Subsample bins (rows)
            if subsample_bins is not None and subsample_bins < len(df):
                print(f"Subsampling {subsample_bins} bins from {len(df)} total bins")
                bin_indices = np.random.choice(
                    len(df), size=subsample_bins, replace=False
                )
                bin_indices = np.sort(bin_indices)  # Keep sorted for interpretability
                df = df.iloc[bin_indices].copy()

            # Subsample samples (columns)
            if subsample_samples is not None and subsample_samples < len(sample_cols):
                print(
                    f"Subsampling {subsample_samples} samples from {len(sample_cols)} total samples"
                )
                sample_indices = np.random.choice(
                    len(sample_cols), size=subsample_samples, replace=False
                )
                selected_samples = [sample_cols[i] for i in sample_indices]
                sample_cols = selected_samples

        # Sort chromosomes in proper order (chr1-chr22, chrX, chrY)
        chr_order = {f"chr{i}": i for i in range(1, 23)}
        chr_order["chrX"] = 23
        chr_order["chrY"] = 24
        df["_chr_order"] = df["Chr"].map(chr_order)
        df = df.sort_values(["_chr_order", "Start"]).drop("_chr_order", axis=1)

        # Extract metadata
        self.chr = df["Chr"].values
        self.start = df["Start"].values
        self.end = df["End"].values
        self.sample_ids = sample_cols

        # Extract normalized read depth matrix (bins x samples)
        depth_matrix = df[sample_cols].values

        # Clamp values above threshold
        if clamp_threshold is not None:
            n_clamped = np.sum(depth_matrix > clamp_threshold)
            if n_clamped > 0:
                print(f"Clamping {n_clamped} values above threshold {clamp_threshold}")
                depth_matrix = np.clip(depth_matrix, None, clamp_threshold)

        # Convert to torch tensors
        self.depth = torch.tensor(depth_matrix, dtype=dtype, device=device)
        self.n_bins = self.depth.shape[0]
        self.n_samples = self.depth.shape[1]

        print(f"Loaded data: {self.n_bins} bins x {self.n_samples} samples")
        print(f"Depth range: [{self.depth.min():.3f}, {self.depth.max():.3f}]")
        print(f"Depth mean: {self.depth.mean():.3f}, median: {self.depth.median():.3f}")


class CNVModel:
    """
    Hierarchical Bayesian model for CNV detection from normalized read depth.

    Model structure:
    - Copy number (CN): discrete latent variable [0,1,2,3,4,5] with prior favoring CN=2
    - Per-bin mean bias: modulates expected depth at each bin
    - Per-sample variance: controls noise in each sample
    - Per-bin variance: controls noise at each bin
    - Observed depth ~ Normal(CN * sample_mean * bin_bias, variance)
    """

    def __init__(
        self,
        n_states: int = 6,
        alpha_ref: float = 50.0,
        alpha_non_ref: float = 1.0,
        var_bias_bin: float = 0.1,
        var_sample: float = 0.2,
        var_bin: float = 0.2,
        device: str = "cpu",
        dtype: torch.dtype = torch.float32,
        debug: bool = False,
        guide_type: str = "diagonal",
    ):
        """
        Args:
            n_states: Number of copy number states (default: 6 for [0,1,2,3,4,5])
            alpha_ref: Dirichlet concentration for reference state (CN=2)
            alpha_non_ref: Dirichlet concentration for non-reference states
            var_bias_bin: Variance for per-bin mean bias (log-normal)
            var_sample: Variance for per-sample variance factor (log-normal)
            var_bin: Variance for per-bin variance factor (log-normal)
            device: torch device
            dtype: torch data type
            debug: Whether to print debug statements in model()
            guide_type: Type of variational guide ('diagonal' for AutoDiagonalNormal, 'delta' for AutoDelta)
        """
        self.n_states = n_states
        self.alpha_ref = alpha_ref
        self.alpha_non_ref = alpha_non_ref
        self.var_bias_bin = var_bias_bin
        self.var_sample = var_sample
        self.var_bin = var_bin
        self.device = device
        self.dtype = dtype
        self.debug = debug
        self.guide_type = guide_type

        # Training history
        self.loss_history = {"epoch": [], "elbo": []}

        # Define which sites to expose to the guide (continuous latent variables)
        self.latent_sites = ["bin_bias", "sample_var", "bin_var", "cn_probs"]

        # Initialize guide based on type
        blocked_model = poutine.block(self.model, expose=self.latent_sites)
        if guide_type == "delta":
            self.guide = AutoDelta(blocked_model)
        elif guide_type == "diagonal":
            self.guide = AutoDiagonalNormal(blocked_model)
        else:
            raise ValueError(
                f"Unknown guide_type: {guide_type}. Choose 'diagonal' or 'delta'."
            )

    @config_enumerate(default="parallel")
    def model(self, depth: torch.Tensor, n_bins: int = None, n_samples: int = None):
        """
        Probabilistic model for CNV detection.

        Args:
            depth: Observed normalized read depth (n_bins x n_samples)
        """

        if self.debug:
            print("\n=== MODEL DEBUG ===")
            print(f"depth.shape: {depth.shape}")

        zero_t = torch.zeros(1, device=self.device, dtype=self.dtype)
        one_t = torch.ones(1, device=self.device, dtype=self.dtype)
        cn_states = torch.arange(0, self.n_states, device=self.device, dtype=self.dtype)

        # Plates for bins and samples
        plate_bins = pyro.plate("bins", n_bins, dim=-2, device=self.device)
        plate_samples = pyro.plate("samples", n_samples, dim=-1, device=self.device)

        # Per-sample variance factor (log-normal prior)
        with plate_samples:
            sample_var = pyro.sample(
                "sample_var", dist.Exponential(1.0 / self.var_sample)
            )
        if self.debug:
            print(f"sample_var.shape: {sample_var.shape}")

        # Per-bin latent variables
        with plate_bins:
            # Per-bin mean bias factor (log-normal prior, centered at 1.0)
            bin_bias = pyro.sample(
                "bin_bias", dist.LogNormal(zero_t, self.var_bias_bin)
            )
            if self.debug:
                print(f"bin_bias.shape: {bin_bias.shape}")

            # Per-bin variance factor (log-normal prior)
            bin_var = pyro.sample("bin_var", dist.Exponential(1.0 / self.var_bin))
            if self.debug:
                print(f"bin_var.shape: {bin_var.shape}")

            # Copy number prior (Dirichlet-Categorical)
            # Heavily weight CN=2 (diploid)
            alpha_cn = self.alpha_non_ref * one_t.expand(self.n_states)
            alpha_cn[2] = self.alpha_ref  # Favor CN=2
            cn_probs = pyro.sample("cn_probs", dist.Dirichlet(alpha_cn))
        if self.debug:
            print(f"cn_probs.shape: {cn_probs.shape}")

        # Per-bin, per-sample copy number and observation
        with plate_bins, plate_samples:
            # Sample copy number (discrete latent variable)
            # Shape: (n_bins, n_samples)
            cn = pyro.sample("cn", dist.Categorical(cn_probs))
            if self.debug:
                print(f"cn.shape: {cn.shape}")

            # Expected depth: CN * bin_bias
            # For diploid (CN=2), expected depth should be ~2.0
            # bin_bias modulates this expectation per bin
            # CN=0 -> 0, CN=1 -> bin_bias, CN=2 -> 2*bin_bias, etc.
            if self.debug:
                print(f"bin_bias.shape: {bin_bias.shape}")
            expected_depth = Vindex(cn_states)[cn] * bin_bias
            if self.debug:
                print(f"expected_depth.shape: {expected_depth.shape}")

            # Variance: combination of sample and bin variance
            # Use summation to combine the two variance factors
            if self.debug:
                print(f"bin_var).shape: {bin_var.shape}")
            variance = sample_var + bin_var
            if self.debug:
                print(f"variance.shape: {variance.shape}")
            std = torch.sqrt(variance)
            if self.debug:
                print(f"std.shape: {std.shape}")

            # Observed depth
            if self.debug:
                print(
                    f"About to sample obs with expected_depth.shape={expected_depth.shape}, std.shape={std.shape}, depth.shape={depth.shape}"
                )
            pyro.sample("obs", dist.Normal(expected_depth, std), obs=depth)
        if self.debug:
            print("=== END MODEL DEBUG ===\n")

    def train(
        self,
        data,
        max_iter: int = 1000,
        lr_init: float = 0.01,
        lr_min: float = 0.001,
        lr_decay: float = 500,
        adam_beta1: float = 0.9,
        adam_beta2: float = 0.999,
        log_freq: int = 50,
        jit: bool = False,
        early_stopping: bool = True,
        patience: int = 50,
        min_delta: float = 1000.0,
    ):
        """
        Train the model using stochastic variational inference (SVI).

        Args:
            data: DepthData object
            max_iter: Maximum number of training iterations
            lr_init: Initial learning rate
            lr_min: Minimum learning rate
            lr_decay: Learning rate decay constant
            adam_beta1: Adam optimizer beta1 parameter
            adam_beta2: Adam optimizer beta2 parameter
            log_freq: Frequency of logging (iterations)
            jit: Whether to use JIT compilation
            early_stopping: Whether to use early stopping
            patience: Number of epochs to wait for improvement before stopping
            min_delta: Minimum change in loss to qualify as improvement
        """
        print("Initializing training...")
        pyro.clear_param_store()

        # Create optimizer with learning rate schedule
        scheduler = pyro.optim.LambdaLR(
            {
                "optimizer": torch.optim.Adam,
                "optim_args": {"lr": 1.0, "betas": (adam_beta1, adam_beta2)},
                "lr_lambda": lambda k: (
                    lr_min + (lr_init - lr_min) * np.exp(-k / lr_decay)
                ),
            }
        )

        # Create ELBO loss
        if jit:
            loss = JitTraceEnum_ELBO()
        else:
            loss = TraceEnum_ELBO()

        # Create SVI object
        svi = SVI(self.model, self.guide, optim=scheduler, loss=loss)

        print(f"Training for up to {max_iter} iterations...")
        if early_stopping:
            print(f"Early stopping enabled: patience={patience}, min_delta={min_delta}")

        # Early stopping tracking
        best_loss = float("inf")
        patience_counter = 0
        stopped_early = False

        try:
            # Configure tqdm with mininterval to force updates
            with tqdm(
                range(max_iter),
                desc="Training",
                unit="epoch",
                mininterval=0.1,
                dynamic_ncols=True,
            ) as pbar:
                for epoch in pbar:
                    # Train step
                    epoch_loss = svi.step(
                        depth=data.depth, n_bins=data.n_bins, n_samples=data.n_samples
                    )
                    scheduler.step()

                    # Record loss
                    self.loss_history["epoch"].append(epoch)
                    self.loss_history["elbo"].append(epoch_loss)

                    # Update progress bar with loss
                    pbar.set_postfix({"loss": f"{epoch_loss:.4f}"})

                    # Log progress
                    if (epoch + 1) % log_freq == 0:
                        tqdm.write(f"[epoch {epoch + 1:04d}]  loss: {epoch_loss:.4f}")

                    # Early stopping check
                    if early_stopping:
                        if epoch_loss < best_loss - min_delta:
                            # Improvement detected
                            best_loss = epoch_loss
                            patience_counter = 0
                        else:
                            # No improvement
                            patience_counter += 1

                        if patience_counter >= patience:
                            tqdm.write(
                                f"\nEarly stopping triggered at epoch {epoch + 1}"
                            )
                            tqdm.write(
                                f"No improvement for {patience} epochs (best loss: {best_loss:.4f})"
                            )
                            stopped_early = True
                            break

            if stopped_early:
                print(f"\nTraining stopped early after {epoch + 1} epochs")
                print(f"Best loss: {best_loss:.4f}, Final loss: {epoch_loss:.4f}")
            else:
                print(
                    f"\nTraining complete after {max_iter} epochs, final loss: {epoch_loss:.4f}"
                )

        except KeyboardInterrupt:
            print("\nTraining interrupted by user.")

    def get_map_estimates(self, data):
        """
        Get MAP (maximum a posteriori) estimates of all latent variables.

        Returns:
            Dictionary with MAP estimates
        """
        print("Computing MAP estimates...")

        # SVI mode: use guide parameters
        # Get guide trace (contains MAP estimates of continuous variables)
        guide_trace = poutine.trace(self.guide).get_trace(
            depth=data.depth, n_bins=data.n_bins, n_samples=data.n_samples
        )

        # Get model trace conditioned on guide
        trained_model = poutine.replay(self.model, trace=guide_trace)

        # Get discrete MAP using infer_discrete
        inferred_model = infer_discrete(
            trained_model, temperature=0, first_available_dim=-3
        )
        trace = poutine.trace(inferred_model).get_trace(
            depth=data.depth, n_bins=data.n_bins, n_samples=data.n_samples
        )

        # Extract all latent variables
        map_estimates = {}

        # Continuous variables from guide
        map_estimates["bin_bias"] = (
            guide_trace.nodes["bin_bias"]["value"].detach().cpu().numpy()
        )
        map_estimates["sample_var"] = (
            guide_trace.nodes["sample_var"]["value"].detach().cpu().numpy()
        )
        map_estimates["bin_var"] = (
            guide_trace.nodes["bin_var"]["value"].detach().cpu().numpy()
        )
        map_estimates["cn_probs"] = (
            guide_trace.nodes["cn_probs"]["value"].detach().cpu().numpy()
        )

        # Discrete variable (copy number)
        map_estimates["cn"] = trace.nodes["cn"]["value"].detach().cpu().numpy()

        return map_estimates

    def run_discrete_inference(self, data, n_samples: int = 1000, log_freq: int = 100):
        """
        Run discrete inference to get posterior distribution over copy numbers.

        Args:
            data: DepthData object
            n_samples: Number of samples for discrete inference
            log_freq: Logging frequency

        Returns:
            Dictionary with copy number posterior probabilities
        """
        print(f"Running discrete inference with {n_samples} samples...")

        # Get guide trace
        guide_trace = poutine.trace(self.guide).get_trace(
            depth=data.depth, n_bins=data.n_bins, n_samples=data.n_samples
        )
        trained_model = poutine.replay(self.model, trace=guide_trace)

        # Accumulate samples
        cn_samples = []

        with torch.no_grad():
            with tqdm(
                range(n_samples),
                desc="Discrete inference",
                unit="sample",
                mininterval=0.1,
                dynamic_ncols=True,
            ) as pbar:
                for i in pbar:
                    inferred_model = infer_discrete(
                        trained_model, temperature=1, first_available_dim=-3
                    )
                    trace = poutine.trace(inferred_model).get_trace(
                        depth=data.depth, n_bins=data.n_bins, n_samples=data.n_samples
                    )
                    cn = trace.nodes["cn"]["value"].detach().cpu().numpy()
                    cn_samples.append(cn)

        # Stack samples and compute frequencies
        cn_samples = np.array(cn_samples)  # (n_samples, n_bins, n_samples)

        # Convert to one-hot encoding and compute frequencies
        cn_freq = np.zeros((data.n_bins, data.n_samples, self.n_states))
        for state in range(self.n_states):
            cn_freq[..., state] = (cn_samples == state).mean(axis=0)

        print("Discrete inference complete.")
        return {"cn_posterior": cn_freq}


def extract_locus_bins(
    df: pd.DataFrame,
    locus: GDLocus,
    segdup_mask: Optional[SegDupMask] = None,
    segdup_threshold: float = 0.5,
    padding: int = 0,
) -> pd.DataFrame:
    """
    Extract bins overlapping a GD locus, optionally masking segmental duplications.

    Args:
        df: DataFrame with all bins
        locus: GDLocus object defining the region
        segdup_mask: Optional SegDupMask for filtering out segdup regions
        segdup_threshold: Minimum overlap fraction with segdups to mask a bin
        padding: Extend the region by this many base pairs on each side

    Returns:
        DataFrame with bins for this locus
    """
    chrom = locus.chrom
    # Get overall start and end from breakpoint ranges
    start = min(bp[0] for bp in locus.breakpoints) - padding
    end = max(bp[1] for bp in locus.breakpoints) + padding

    # Filter to bins in this region
    mask = (
        (df["Chr"] == chrom) &
        (df["End"] > start) &
        (df["Start"] < end)
    )
    locus_df = df[mask].copy()

    if len(locus_df) == 0:
        return locus_df

    # Apply segdup masking if provided
    if segdup_mask is not None:
        keep_mask = []
        for _, row in locus_df.iterrows():
            overlap = segdup_mask.get_overlap_fraction(
                row["Chr"], row["Start"], row["End"]
            )
            keep_mask.append(overlap < segdup_threshold)

        n_masked = len(keep_mask) - sum(keep_mask)
        if n_masked > 0:
            print(f"  Masked {n_masked}/{len(locus_df)} bins due to segdup overlap")

        locus_df = locus_df[keep_mask]

    return locus_df


def compute_flank_regions_from_bins(
    locus_df: pd.DataFrame,
    locus: GDLocus,
    target_size: int,
) -> List[Tuple[int, int, str]]:
    """
    Compute flanking regions based on actual available (unfiltered) bins.

    Accumulates bins outward from each side of the locus until the total covered
    base pairs of unfiltered bins reaches target_size. This avoids the problem of
    large segmental duplication deserts making a fixed-size geometric window empty.

    Args:
        locus_df: DataFrame of unfiltered bins available for this locus
        locus: GDLocus object
        target_size: Target cumulative bin coverage in bp for each flank

    Returns:
        List of (start, end, name) tuples for "left_flank" and/or "right_flank".
        A flank is omitted if no bins exist on that side.
    """
    flanks = []
    locus_start = locus.start
    locus_end = locus.end

    print(f"    [flanks] locus body: {locus_start:,}-{locus_end:,}  target_size: {target_size:,} bp")
    print(f"    [flanks] input bins: {len(locus_df)}  "
          f"range: {int(locus_df['Start'].min()):,}-{int(locus_df['End'].max()):,}" if len(locus_df) > 0
          else f"    [flanks] input bins: 0")

    # Use bin midpoints for candidacy — consistent with how bins are assigned
    # throughout the rest of the codebase. Using strict End/Start comparisons
    # misses bins that straddle the locus boundary (common with zero-width breakpoints).
    bin_mids = (locus_df["Start"] + locus_df["End"]) / 2

    # Left flank: bins whose midpoint is left of locus start, accumulated right-to-left
    left_bins = locus_df[bin_mids < locus_start].sort_values("Start", ascending=False)
    print(f"    [flanks] left candidate bins (mid < {locus_start:,}): {len(left_bins)}", end="")
    if len(left_bins) > 0:
        print(f"  span {int(left_bins['Start'].min()):,}-{int(left_bins['End'].max()):,}  "
              f"total bp: {int((left_bins['End'] - left_bins['Start']).sum()):,}")
        cumulative = 0
        leftmost_start = locus_start  # will be updated as we walk left
        for _, row in left_bins.iterrows():
            cumulative += int(row["End"]) - int(row["Start"])
            leftmost_start = int(row["Start"])
            if cumulative >= target_size:
                break
        print(f"    [flanks] left_flank: {leftmost_start:,}-{locus_start:,}  "
              f"accumulated {cumulative:,} bp (target {target_size:,})")
        flanks.append((leftmost_start, locus_start, "left_flank"))
    else:
        print()  # newline after the count

    # Right flank: bins whose midpoint is right of locus end, accumulated left-to-right
    right_bins = locus_df[bin_mids >= locus_end].sort_values("Start")
    print(f"    [flanks] right candidate bins (mid >= {locus_end:,}): {len(right_bins)}", end="")
    if len(right_bins) > 0:
        print(f"  span {int(right_bins['Start'].min()):,}-{int(right_bins['End'].max()):,}  "
              f"total bp: {int((right_bins['End'] - right_bins['Start']).sum()):,}")
        cumulative = 0
        rightmost_end = locus_end  # will be updated as we walk right
        for _, row in right_bins.iterrows():
            cumulative += int(row["End"]) - int(row["Start"])
            rightmost_end = int(row["End"])
            if cumulative >= target_size:
                break
        print(f"    [flanks] right_flank: {locus_end:,}-{rightmost_end:,}  "
              f"accumulated {cumulative:,} bp (target {target_size:,})")
        flanks.append((locus_end, rightmost_end, "right_flank"))
    else:
        print()  # newline after the count

    return flanks


def rebin_locus_intervals(
    df: pd.DataFrame,
    locus: GDLocus,
    max_bins_per_interval: int = 10,
    flank_regions: Optional[List[Tuple[int, int, str]]] = None,
) -> pd.DataFrame:
    """
    Reduce the number of bins by rebinning each interval and flanking region to have
    at most max_bins_per_interval bins. Uses weighted averaging based on bin size.

    Applies to:
      - Intervals between adjacent breakpoints
      - Left and right flanking regions (bin-derived if flank_regions is provided,
        otherwise falls back to geometric 100%-of-locus-size windows)
      - Each breakpoint range (SD blocks) independently

    Args:
        df: DataFrame with bins for this locus
        locus: GDLocus object
        max_bins_per_interval: Maximum number of bins per interval (default: 10)

    Returns:
        DataFrame with rebinned data
    """
    if len(df) == 0:
        return df

    # Get sample columns
    metadata_cols = ["Chr", "Start", "End", "source_file", "Bin"]
    sample_cols = [col for col in df.columns if col not in metadata_cols]

    def rebin_region(region_df: pd.DataFrame) -> List[dict]:
        """Rebin a set of bins to at most max_bins_per_interval using weighted averaging."""
        if len(region_df) <= max_bins_per_interval:
            return region_df.to_dict("records")

        n_new_bins = max_bins_per_interval
        bins_per_group = len(region_df) / n_new_bins
        region_df = region_df.sort_values("Start")
        rows = []

        for i in range(n_new_bins):
            start_idx = int(i * bins_per_group)
            end_idx = int((i + 1) * bins_per_group) if i < n_new_bins - 1 else len(region_df)
            group = region_df.iloc[start_idx:end_idx]
            if len(group) == 0:
                continue

            bin_sizes = (group["End"] - group["Start"]).values
            weights = bin_sizes / bin_sizes.sum()
            new_bin = {
                "Chr": group.iloc[0]["Chr"],
                "Start": int(group["Start"].min()),
                "End": int(group["End"].max()),
            }
            for col in sample_cols:
                new_bin[col] = (group[col].values * weights).sum()
            if "source_file" in group.columns:
                new_bin["source_file"] = group.iloc[0]["source_file"]
            rows.append(new_bin)

        return rows

    # All rebinnable regions: between-breakpoint intervals + flanking regions
    all_rebin_regions = locus.get_intervals() + (
        flank_regions if flank_regions is not None else locus.get_flanking_regions()
    )

    rebinned_rows = []
    processed_indices = set()

    for region_start, region_end, _ in all_rebin_regions:
        region_indices = []
        for idx, row in df.iterrows():
            bin_mid = (row["Start"] + row["End"]) / 2
            if region_start <= bin_mid < region_end:
                region_indices.append(idx)
                processed_indices.add(idx)

        if not region_indices:
            continue

        region_df = df.loc[region_indices].copy()
        rebinned_rows.extend(rebin_region(region_df))

    # Rebin bins inside each breakpoint range (SD blocks) independently
    for bp_start, bp_end in locus.breakpoints:
        bp_indices = []
        for idx, row in df.iterrows():
            if idx in processed_indices:
                continue
            bin_mid = (row["Start"] + row["End"]) / 2
            if bp_start <= bin_mid < bp_end:
                bp_indices.append(idx)
                processed_indices.add(idx)

        if bp_indices:
            rebinned_rows.extend(rebin_region(df.loc[bp_indices].copy()))

    # Pass through any remaining bins not covered by any region or breakpoint
    for idx, row in df.iterrows():
        if idx not in processed_indices:
            rebinned_rows.append(row.to_dict())

    if len(rebinned_rows) == 0:
        return pd.DataFrame(columns=df.columns)

    return pd.DataFrame(rebinned_rows)


def assign_bins_to_intervals(
    df: pd.DataFrame,
    locus: GDLocus,
    flank_regions: Optional[List[Tuple[int, int, str]]] = None,
) -> Dict[str, List[int]]:
    """
    Assign bins to breakpoint intervals and flanking regions within a locus.

    Args:
        df: DataFrame with bins for this locus
        locus: GDLocus object
        flank_regions: Pre-computed bin-derived flank regions. If None, falls back
            to geometric windows from locus.get_flanking_regions().

    Returns:
        Dict mapping region name to list of bin indices. Includes between-breakpoint
        interval names (e.g. "1-2"), flanking names ("left_flank", "right_flank"),
        and "breakpoint_ranges" for bins that fall inside the breakpoint SD blocks.
    """
    all_named_regions = locus.get_intervals() + (
        flank_regions if flank_regions is not None else locus.get_flanking_regions()
    )
    interval_bins: Dict[str, List[int]] = {name: [] for _, _, name in all_named_regions}
    interval_bins["breakpoint_ranges"] = []

    for idx, row in df.iterrows():
        bin_mid = (row["Start"] + row["End"]) / 2

        matched = False
        for start, end, name in all_named_regions:
            if start <= bin_mid < end:
                interval_bins[name].append(idx)
                matched = True
                break

        if not matched:
            interval_bins["breakpoint_ranges"].append(idx)

    return interval_bins


def compute_interval_cn_stats(
    cn_posterior: np.ndarray,
    interval_bins: Dict[str, List[int]],
    sample_idx: int,
) -> Dict[str, dict]:
    """
    Compute copy number statistics for each interval.

    Args:
        cn_posterior: Array of shape (n_bins, n_samples, n_states)
        cn_map: Array of shape (n_bins, n_samples) with MAP CN estimates
        interval_bins: Dict mapping interval name to bin indices
        sample_idx: Sample index

    Returns:
        Dict mapping interval name to CN statistics
    """
    result = {}

    for interval_name, bin_indices in interval_bins.items():
        if len(bin_indices) == 0:
            result[interval_name] = {
                "n_bins": 0,
                "cn_probs": np.zeros(6),
            }
            continue

        # Get CN values for this interval
        interval_cn_probs = cn_posterior[bin_indices, sample_idx, :]
        
        # Average posterior probabilities across bins
        mean_probs = interval_cn_probs.mean(axis=0)

        result[interval_name] = {
            "n_bins": len(bin_indices),
            "cn_probs": mean_probs,
        }

    return result


def call_gd_cnv(
    locus: GDLocus,
    interval_stats: Dict[str, dict],
    log_prob_threshold: float = -0.5) -> List[dict]:
    """
    Call GD CNVs based on interval copy number statistics.

    For each GD entry in the locus (each DEL/DUP definition), check if the
    copy number in the corresponding region supports a CNV call.

    Args:
        locus: GDLocus object
        interval_stats: Dict mapping interval name to CN statistics
        log_prob_threshold: Minimum log probability score to call a CNV

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
        # The breakpoint names in BP1 and BP2 define the start and end breakpoints
        # We need to find all intervals between these two breakpoints
        
        # Get all interval names for this locus
        all_intervals = locus.get_intervals()  # List of (start, end, name) tuples
        
        # Build a mapping of breakpoint names to their positions
        # Interval names are like "BP1-2", "BP2-3", etc. or "A-B", "B-C", etc.
        # Extract the breakpoint ordering from the interval names
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
            # BP1 or BP2 not found in ordering - skip this entry
            print(f"  Warning: Could not find {bp1} or {bp2} in breakpoint order for {gd_id}")
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
                    # Include interval if it falls between pos1 and pos2
                    if start_pos >= pos1 and end_pos <= pos2:
                        covered_intervals.append(interval_name)
                except ValueError:
                    continue

        if len(covered_intervals) == 0:
            continue

        # Determine if this is a CNV using weighted average of log probabilities
        # Weight by interval size (number of bins)
        log_prob_score = 0.0
        total_weight = 0.0
        
        if svtype == "DEL":
            # For affected intervals: weighted average of log P(CN < 2)
            for interval in covered_intervals:
                if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                    probs = interval_stats[interval]["cn_probs"]
                    weight = interval_stats[interval]["n_bins"]
                    p_del = max(probs[0] + probs[1], 1e-3)  # P(CN=0) + P(CN=1)
                    if p_del > 0:
                        log_prob_score += weight * np.log(p_del)
                        total_weight += weight
        
        elif svtype == "DUP":
            # For affected intervals: weighted average of log P(CN > 2)
            for interval in covered_intervals:
                if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                    probs = interval_stats[interval]["cn_probs"]
                    weight = interval_stats[interval]["n_bins"]
                    p_dup = max(probs[3:].sum(), 1e-3)  # P(CN >= 3)
                    if p_dup > 0:
                        log_prob_score += weight * np.log(p_dup)
                        total_weight += weight
        
        # Normalize by total weight to get weighted average
        if total_weight > 0:
            log_prob_score = log_prob_score / total_weight
        else:
            log_prob_score = np.nan  # No data to support call
        
        # Determine if this is a CNV based on log probability score
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
                    p_flank = max(flank_probs[0] + flank_probs[1], 1e-3)
                elif svtype == "DUP":
                    p_flank = max(flank_probs[3:].sum(), 1e-3)
                else:
                    continue
                flank_score_sum += flank_weight * np.log(p_flank)
                flank_weight_total += flank_weight

            if flank_weight_total > 0:
                flanking_log_prob_score = flank_score_sum / flank_weight_total
                # If flanking regions also show a CN change, this is a spanning variant
                if flanking_log_prob_score > log_prob_threshold:
                    is_spanning = True
                    is_cnv = False

        # Count total bins across covered intervals
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
            "n_bins": n_bins
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

    Scoring approach:
    - For DEL: log-sum P(CN < 2) in affected intervals + log-sum P(CN >= 2) in unaffected intervals
    - For DUP: log-sum P(CN > 2) in affected intervals + log-sum P(CN <= 2) in unaffected intervals

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

        # Multiple carrier calls of same svtype - score each based on CN probabilities
        best_score = -np.inf
        best_gd_id = None

        for call in carrier_calls:
            covered_intervals = set(call["intervals"])
            uncovered_intervals = all_intervals - covered_intervals

            # Compute weighted average score based on svtype
            score = 0.0
            total_weight = 0.0

            if svtype == "DEL":
                # For affected intervals: weighted average of log P(CN < 2)
                for interval in covered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = interval_stats[interval]["n_bins"]
                        p_del = probs[0] + probs[1]  # P(CN=0) + P(CN=1)
                        if p_del > 0:
                            score += weight * np.log(p_del)
                            total_weight += weight

                # For unaffected intervals: weighted average of log P(CN >= 2)
                for interval in uncovered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = interval_stats[interval]["n_bins"]
                        p_normal = probs[2:].sum()  # P(CN >= 2)
                        if p_normal > 0:
                            score += weight * np.log(p_normal)
                            total_weight += weight

            else:  # DUP
                # For affected intervals: weighted average of log P(CN > 2)
                for interval in covered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = interval_stats[interval]["n_bins"]
                        p_dup = probs[3:].sum()  # P(CN >= 3)
                        if p_dup > 0:
                            score += weight * np.log(p_dup)
                            total_weight += weight

                # For unaffected intervals: weighted average of log P(CN <= 2)
                for interval in uncovered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = interval_stats[interval]["n_bins"]
                        p_normal = probs[:3].sum()  # P(CN <= 2)
                        if p_normal > 0:
                            score += weight * np.log(p_normal)
                            total_weight += weight
            
            # Normalize by total weight
            if total_weight > 0:
                score = score / total_weight

            # Update best if this score is higher
            if score > best_score:
                best_score = score
                best_gd_id = call["GD_ID"]

        best_by_svtype[svtype] = best_gd_id

    return best_by_svtype


def read_data(file_path: str) -> pd.DataFrame:
    """Load binned read count data from TSV file."""
    print(f"Loading: {file_path}")
    try:
        # Read gzipped TSV file
        df_file = pd.read_csv(file_path, sep="\t", compression="infer")

        # Get source file identifier
        source_file = str(Path(file_path).parent.parent.name)

        # Rename chromosome column if it exists
        if "#Chr" in df_file.columns:
            df_file["Chr"] = df_file["#Chr"]
            df_file = df_file.drop("#Chr", axis=1)

        # Create bin identifier
        df_file["Bin"] = (
            df_file["Chr"].astype(str) +
            ":" +
            df_file["Start"].astype(str) +
            "-" +
            df_file["End"].astype(str)
        )
        df_file = df_file.set_index("Bin")
        df_file["source_file"] = source_file
        return df_file
    except Exception as e:
        print(f"Error loading {file_path}")
        raise e


def filter_low_quality_bins(
    df: pd.DataFrame,
    median_min: float = 1.5,
    median_max: float = 2.5,
    mad_max: float = 0.5,
):
    """
    Filter out low quality bins based on median and MAD thresholds.

    Args:
        df: DataFrame with bins as rows and samples as columns
        median_min: Minimum median depth for bins
        median_max: Maximum median depth for bins
        mad_max: Maximum MAD for bins

    Returns:
        Filtered DataFrame
    """
    # Get sample columns (exclude metadata)
    sample_cols = get_sample_columns(df)

    # Compute median and MAD for each bin across samples
    depths = df[sample_cols].values
    medians = np.median(depths, axis=1)
    mads = np.median(np.abs(depths - medians[:, np.newaxis]), axis=1)

    print(f"\n{'=' * 80}")
    print("FILTERING LOW QUALITY BINS")
    print(f"{'=' * 80}")
    print(f"Starting bins: {len(df)}")

    # Filter based on thresholds
    keep_mask = (
        (medians >= median_min) &
        (medians <= median_max) &
        (mads <= mad_max)
    )

    n_filtered = (~keep_mask).sum()
    print(f"Thresholds: median [{median_min}, {median_max}], MAD <= {mad_max}")
    print(f"Bins filtered: {n_filtered}")
    print(f"Bins remaining: {keep_mask.sum()}")
    print(f"{'=' * 80}\n")

    return df[keep_mask].copy()


@dataclass
class LocusBinMapping:
    """Tracks mapping between array indices and locus/interval assignments."""
    cluster: str
    locus: GDLocus
    interval_name: str
    array_idx: int  # Index in the combined depth tensor
    chrom: str  # Chromosome
    start: int  # Bin start position
    end: int  # Bin end position


def collect_all_locus_bins(
    df: pd.DataFrame,
    gd_table: GDTable,
    segdup_mask: Optional[SegDupMask],
    segdup_threshold: float = 0.5,
    locus_padding: int = 0,
    min_bins_per_locus: int = 5,
    max_bins_per_interval: int = 10,
) -> Tuple[pd.DataFrame, List[LocusBinMapping], Dict[str, GDLocus]]:
    """
    Collect all bins across all GD loci into a single DataFrame.

    Args:
        df: DataFrame with all bins
        gd_table: GDTable with locus definitions
        segdup_mask: Optional SegDupMask for filtering
        segdup_threshold: Minimum overlap fraction with segdups to mask a bin
        locus_padding: Padding around locus boundaries
        min_bins_per_locus: Minimum bins required to include a locus
        max_bins_per_interval: Maximum bins per interval after rebinning (0 = no rebinning)

    Returns:
        Tuple of:
        - Combined DataFrame with all locus bins
        - List of LocusBinMapping objects tracking bin assignments
        - Dict of included loci (cluster -> GDLocus)
    """
    all_locus_dfs = []
    all_mappings = []
    included_loci = {}
    current_idx = 0

    print(f"\n{'=' * 80}")
    print("COLLECTING BINS ACROSS ALL GD LOCI")
    print(f"{'=' * 80}")

    # Pre-build a per-chromosome cache of filtered (segdup-masked) bins.
    # Flank computation walks outward through these; using the full chromosome
    # guarantees filtered bins are found even when multi-megabase segdup deserts
    # surround the locus.
    print("\nBuilding per-chromosome filtered bin cache...")
    chrom_filtered: Dict[str, pd.DataFrame] = {}
    for chrom, chrom_df in df.groupby("Chr"):
        if segdup_mask is not None:
            keep = []
            for _, row in chrom_df.iterrows():
                overlap = segdup_mask.get_overlap_fraction(chrom, row["Start"], row["End"])
                keep.append(overlap < segdup_threshold)
            chrom_filtered[chrom] = chrom_df[keep].copy()
        else:
            chrom_filtered[chrom] = chrom_df.copy()
        print(f"  {chrom}: {len(chrom_filtered[chrom])} filtered bins "
              f"(of {len(chrom_df)} total)")

    for cluster, locus in gd_table.get_all_loci().items():
        print(f"\nProcessing locus: {cluster}")
        print(f"  Chromosome: {locus.chrom}")
        print(f"  Breakpoints: {locus.breakpoints}")
        print(f"  GD entries: {len(locus.gd_entries)} ({', '.join(locus.svtypes)})")

        locus_size = locus.end - locus.start

        # Compute flank coordinates from the full chromosome's filtered bins so
        # that even multi-megabase segdup deserts don't prevent flank discovery.
        chrom_bins = chrom_filtered.get(locus.chrom, pd.DataFrame())
        flank_regions = compute_flank_regions_from_bins(chrom_bins, locus, locus_size)
        if flank_regions:
            for fs, fe, fn in flank_regions:
                print(f"  {fn}: {fs:,}-{fe:,} (bin-derived)")
        else:
            print(f"  Warning: no flanking bins found for locus {cluster}")

        # Determine extraction bounds: locus body + computed flank extents.
        left_bound = locus.start
        right_bound = locus.end
        for fs, fe, fn in flank_regions:
            if fn == "left_flank":
                left_bound = fs
            elif fn == "right_flank":
                right_bound = fe

        # Extract only the bins within the active region (locus + flanks),
        # applying segdup masking. locus_padding still applies as a minimum.
        active_padding = max(locus_padding, locus.start - left_bound, right_bound - locus.end)
        locus_df = extract_locus_bins(
            df, locus, segdup_mask,
            segdup_threshold=segdup_threshold,
            padding=active_padding,
        )

        if len(locus_df) < min_bins_per_locus:
            print(f"  Skipping: only {len(locus_df)} bins (minimum {min_bins_per_locus})")
            continue

        print(f"  Bins after filtering: {len(locus_df)}")

        # Trim to active region: drop any bins outside [left_bound, right_bound)
        bin_mids = (locus_df["Start"] + locus_df["End"]) / 2
        locus_df = locus_df[(bin_mids >= left_bound) & (bin_mids < right_bound)].copy()
        print(f"  Bins after trimming to active region [{left_bound:,}, {right_bound:,}): {len(locus_df)}")

        # Rebin to reduce number of bins per interval/flank if requested
        if max_bins_per_interval > 0:
            locus_df_orig = locus_df
            locus_df = rebin_locus_intervals(locus_df, locus, max_bins_per_interval, flank_regions)
            if len(locus_df) < len(locus_df_orig):
                print(f"  Bins after rebinning: {len(locus_df)} (reduced from {len(locus_df_orig)})")

        # Assign bins to intervals and flanking regions
        interval_bins = assign_bins_to_intervals(locus_df, locus, flank_regions)
        total_assigned = sum(len(v) for v in interval_bins.values())
        for region_name, bins in interval_bins.items():
            print(f"    {region_name}: {len(bins)} bins")
        print(f"    total: {total_assigned} bins")

        # Build a fast lookup: bin midpoint -> assigned region name
        all_named_regions = locus.get_intervals() + flank_regions

        # Create mappings for each bin
        for idx in locus_df.index:
            bin_row = locus_df.loc[idx]
            bin_mid = (bin_row["Start"] + bin_row["End"]) / 2

            # Check named regions (intervals + flanks) in order
            assigned_interval = None
            for start, end, name in all_named_regions:
                if start <= bin_mid < end:
                    assigned_interval = name
                    break

            # Bins inside breakpoint ranges: assign to the nearest between-breakpoint interval
            if assigned_interval is None:
                intervals = locus.get_intervals()
                if intervals:
                    if bin_mid < intervals[0][0]:
                        assigned_interval = intervals[0][2]
                    elif bin_mid >= intervals[-1][1]:
                        assigned_interval = intervals[-1][2]
                    else:
                        for i in range(len(intervals) - 1):
                            if intervals[i][1] <= bin_mid < intervals[i + 1][0]:
                                dist_left = bin_mid - intervals[i][1]
                                dist_right = intervals[i + 1][0] - bin_mid
                                assigned_interval = intervals[i][2] if dist_left <= dist_right else intervals[i + 1][2]
                                break

            all_mappings.append(LocusBinMapping(
                cluster=cluster,
                locus=locus,
                interval_name=assigned_interval or (all_named_regions[0][2] if all_named_regions else "unknown"),
                array_idx=current_idx,
                chrom=bin_row["Chr"],
                start=int(bin_row["Start"]),
                end=int(bin_row["End"]),
            ))
            current_idx += 1

        all_locus_dfs.append(locus_df)
        included_loci[cluster] = locus

    if len(all_locus_dfs) == 0:
        print("\nNo loci with sufficient bins found!")
        return pd.DataFrame(), [], {}

    # Combine all DataFrames
    combined_df = pd.concat(all_locus_dfs, axis=0)

    # Reset index to ensure contiguous integer indices
    combined_df = combined_df.reset_index(drop=True)

    print(f"\n{'=' * 80}")
    print(f"TOTAL: {len(combined_df)} bins across {len(included_loci)} loci")
    print(f"{'=' * 80}\n")

    return combined_df, all_mappings, included_loci


def get_locus_interval_bins(
    mappings: List[LocusBinMapping],
    cluster: str,
) -> Dict[str, List[int]]:
    """
    Get array indices for each interval in a specific locus.

    Args:
        mappings: List of LocusBinMapping objects
        cluster: Cluster name to filter by

    Returns:
        Dict mapping interval name to list of array indices
    """
    interval_bins = {}
    for mapping in mappings:
        if mapping.cluster != cluster:
            continue
        if mapping.interval_name not in interval_bins:
            interval_bins[mapping.interval_name] = []
        interval_bins[mapping.interval_name].append(mapping.array_idx)
    return interval_bins


def write_posterior_tables(
    combined_data: "DepthData",
    map_estimates: dict,
    cn_posterior: dict,
    mappings: List[LocusBinMapping],
    output_dir: str,
):
    """
    Write comprehensive posterior tables to disk.
    
    Args:
        combined_data: DepthData object with all bins
        map_estimates: Dictionary with MAP estimates from model
        cn_posterior: Dictionary with CN posterior probabilities
        mappings: List of LocusBinMapping objects
        output_dir: Output directory for tables
    """
    print("\n" + "=" * 80)
    print("WRITING POSTERIOR TABLES")
    print("=" * 80)
    
    # 1. Copy state probabilities for all bins and samples
    print("\nWriting copy state posteriors...")
    cn_post = np.asarray(cn_posterior["cn_posterior"]).squeeze()  # shape: (n_bins, n_samples, n_states)
    cn_map = np.asarray(map_estimates["cn"]).squeeze()  # shape: (n_bins, n_samples)
    depth = np.asarray(combined_data.depth.cpu().numpy())  # shape: (n_bins, n_samples)
    
    # Ensure proper dimensions
    if cn_post.ndim == 2:
        cn_post = cn_post.reshape(cn_post.shape[0], cn_post.shape[1], 1)
    if cn_map.ndim == 1:
        cn_map = cn_map.reshape(-1, 1)
    if depth.ndim == 1:
        depth = depth.reshape(-1, 1)
    
    cn_rows = []
    for bin_idx in range(combined_data.n_bins):
        mapping = mappings[bin_idx]
        
        for sample_idx, sample_id in enumerate(combined_data.sample_ids):
            row = {
                "cluster": mapping.cluster,
                "interval": mapping.interval_name,
                "chr": mapping.chrom,
                "start": mapping.start,
                "end": mapping.end,
                "sample": sample_id,
                "depth": depth[bin_idx, sample_idx].tolist() if isinstance(depth[bin_idx, sample_idx], np.ndarray) else float(depth[bin_idx, sample_idx]),
            }
            
            # Add probability for each CN state
            for cn_state in range(cn_post.shape[2]):
                prob_val = cn_post[bin_idx, sample_idx, cn_state]
                row[f"prob_cn_{cn_state}"] = prob_val.tolist() if isinstance(prob_val, np.ndarray) else float(prob_val)
            
            # Add MAP estimate
            map_val = cn_map[bin_idx, sample_idx]
            row["cn_map"] = int(map_val.tolist() if isinstance(map_val, np.ndarray) else map_val)
            
            cn_rows.append(row)
    
    cn_df = pd.DataFrame(cn_rows)
    cn_output = os.path.join(output_dir, "cn_posteriors.tsv.gz")
    cn_df.to_csv(cn_output, sep="\t", index=False, compression="gzip")
    print(f"  Saved: {cn_output}")
    print(f"  Rows: {len(cn_df):,} ({combined_data.n_bins:,} bins × {combined_data.n_samples} samples)")
    
    # 2. Sample-specific variable posteriors
    print("\nWriting sample-specific variable posteriors...")
    sample_rows = []
    
    # Convert to numpy array and squeeze extra dimensions
    sample_var = np.asarray(map_estimates["sample_var"]).squeeze()
    
    # Ensure it's at least 1D
    if sample_var.ndim == 0:
        sample_var = sample_var.reshape(1)
    
    for sample_idx, sample_id in enumerate(combined_data.sample_ids):
        var_val = sample_var[sample_idx]
        row = {
            "sample": sample_id,
            "sample_var_map": var_val.tolist() if isinstance(var_val, np.ndarray) else float(var_val),
        }
        sample_rows.append(row)
    
    sample_df = pd.DataFrame(sample_rows)
    sample_output = os.path.join(output_dir, "sample_posteriors.tsv.gz")
    sample_df.to_csv(sample_output, sep="\t", index=False, compression="gzip")
    print(f"  Saved: {sample_output}")
    print(f"  Rows: {len(sample_df):,} ({combined_data.n_samples} samples)")
    
    # 3. Bin-specific variable posteriors
    print("\nWriting bin-specific variable posteriors...")
    bin_rows = []
    
    # Convert to numpy arrays and ensure proper shape
    bin_bias = np.asarray(map_estimates["bin_bias"]).squeeze()
    bin_var = np.asarray(map_estimates["bin_var"]).squeeze()
    cn_probs = np.asarray(map_estimates["cn_probs"]).squeeze()
    
    # Ensure we have the right number of dimensions
    if bin_bias.ndim == 0:
        bin_bias = bin_bias.reshape(1)
    if bin_var.ndim == 0:
        bin_var = bin_var.reshape(1)
    if cn_probs.ndim == 1:
        cn_probs = cn_probs.reshape(-1, 1)
    
    for bin_idx in range(combined_data.n_bins):
        mapping = mappings[bin_idx]
        
        row = {
            "cluster": mapping.cluster,
            "interval": mapping.interval_name,
            "chr": mapping.chrom,
            "start": mapping.start,
            "end": mapping.end,
            "bin_bias_map": bin_bias[bin_idx].tolist() if isinstance(bin_bias[bin_idx], np.ndarray) else float(bin_bias[bin_idx]),
            "bin_var_map": bin_var[bin_idx].tolist() if isinstance(bin_var[bin_idx], np.ndarray) else float(bin_var[bin_idx]),
        }
        
        # Add CN probability priors (per-bin learned from data)
        for cn_state in range(cn_probs.shape[1]):
            prob_val = cn_probs[bin_idx, cn_state]
            row[f"cn_prior_{cn_state}"] = prob_val.tolist() if isinstance(prob_val, np.ndarray) else float(prob_val)
        
        bin_rows.append(row)
    
    bin_df = pd.DataFrame(bin_rows)
    bin_output = os.path.join(output_dir, "bin_posteriors.tsv.gz")
    bin_df.to_csv(bin_output, sep="\t", index=False, compression="gzip")
    print(f"  Saved: {bin_output}")
    print(f"  Rows: {len(bin_df):,} ({combined_data.n_bins:,} bins)")
    
    print("\n" + "=" * 80)


def write_locus_metadata(
    included_loci: Dict[str, GDLocus],
    mappings: List[LocusBinMapping],
    output_dir: str,
):
    """
    Write locus metadata and bin mappings for use by downstream scripts.
    
    Args:
        included_loci: Dict of cluster -> GDLocus objects
        mappings: List of LocusBinMapping objects
        output_dir: Output directory
    """
    print("\nWriting locus metadata...")
    
    # 1. Write bin-to-interval mappings
    bin_mapping_rows = []
    for mapping in mappings:
        bin_mapping_rows.append({
            "cluster": mapping.cluster,
            "interval": mapping.interval_name,
            "chr": mapping.chrom,
            "start": mapping.start,
            "end": mapping.end,
            "array_idx": mapping.array_idx,
        })
    
    bin_mapping_df = pd.DataFrame(bin_mapping_rows)
    bin_mapping_output = os.path.join(output_dir, "bin_mappings.tsv.gz")
    bin_mapping_df.to_csv(bin_mapping_output, sep="\t", index=False, compression="gzip")
    print(f"  Saved: {bin_mapping_output}")
    print(f"  Rows: {len(bin_mapping_df):,} bins")
    
    # 2. Write locus definitions with interval coordinates
    locus_rows = []
    for cluster, locus in included_loci.items():
        # Get intervals for this locus
        for start, end, name in locus.get_intervals():
            locus_rows.append({
                "cluster": cluster,
                "interval": name,
                "chr": locus.chrom,
                "start": start,
                "end": end,
            })
    
    locus_df = pd.DataFrame(locus_rows)
    locus_output = os.path.join(output_dir, "locus_intervals.tsv.gz")
    locus_df.to_csv(locus_output, sep="\t", index=False, compression="gzip")
    print(f"  Saved: {locus_output}")
    print(f"  Rows: {len(locus_df):,} intervals")


def run_gd_analysis(
    df: pd.DataFrame,
    gd_table: GDTable,
    segdup_mask: Optional[SegDupMask],
    args: argparse.Namespace,
    device: str = "cpu",
):
    """
    Run GD CNV analysis on all loci using a single unified model.
    
    This function performs model training and inference only.
    CNV calling is handled by downstream scripts (plot_gd_cnv_output.py).
    
    This function performs model training and inference only.
    CNV calling is handled by downstream scripts (plot_gd_cnv_output.py).

    This function collects bins from all GD loci, trains a single model on all
    bins together, and writes out posterior probabilities and metadata.

    Args:
        df: DataFrame with normalized read depth
        gd_table: GDTable with locus definitions
        segdup_mask: Optional SegDupMask for filtering
        args: Command line arguments
        device: Torch device
    """
    # Collect all bins across all loci
    combined_df, mappings, included_loci = collect_all_locus_bins(
        df, gd_table, segdup_mask,
        segdup_threshold=args.segdup_threshold,
        locus_padding=args.locus_padding,
        min_bins_per_locus=args.min_bins_per_locus,
        max_bins_per_interval=args.max_bins_per_interval,
    )

    if len(combined_df) == 0:
        print("No bins to analyze!")
        return pd.DataFrame()

    # Create data object for all loci combined
    print("\nCreating combined depth data...")
    combined_data = DepthData(
        combined_df,
        device=device,
        dtype=torch.float32,
        clamp_threshold=args.clamp_threshold,
    )

    # Initialize and train a single model on all bins
    print("\nInitializing unified CNV model...")
    model = CNVModel(
        n_states=6,
        alpha_ref=args.alpha_ref,
        alpha_non_ref=args.alpha_non_ref,
        var_bias_bin=args.var_bias_bin,
        var_sample=args.var_sample,
        var_bin=args.var_bin,
        device=device,
        dtype=torch.float32,
        guide_type=args.guide_type,
    )

    print("\nTraining unified model on all GD loci bins...")
    model.train(
        data=combined_data,
        max_iter=args.max_iter,
        lr_init=args.lr_init,
        lr_min=args.lr_min,
        lr_decay=args.lr_decay,
        log_freq=args.log_freq,
        jit=args.jit,
        early_stopping=args.early_stopping,
        patience=args.patience,
        min_delta=args.min_delta,
    )

    # Get MAP estimates and posterior for all bins
    print("\nComputing MAP estimates...")
    map_estimates = model.get_map_estimates(combined_data)

    print("\nRunning discrete inference...")
    cn_posterior = model.run_discrete_inference(
        combined_data,
        n_samples=args.n_discrete_samples,
        log_freq=args.log_freq,
    )

    # Write comprehensive posterior tables
    write_posterior_tables(
        combined_data,
        map_estimates,
        cn_posterior,
        mappings,
        args.output_dir,
    )
    
    # Write locus metadata for downstream calling/plotting
    write_locus_metadata(
        included_loci,
        mappings,
        args.output_dir,
    )
    
    print("\n" + "=" * 80)
    print("Model training and inference complete!")
    print("Use plot_gd_cnv_output.py to call CNVs and generate plots.")
    print("=" * 80)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Bayesian genomic disorder CNV detection from binned read depth",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Input/Output
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input TSV file with normalized read depth (bins x samples)",
    )
    parser.add_argument(
        "-g", "--gd-table",
        required=True,
        help="GD locus definition table (TSV)",
    )
    parser.add_argument(
        "-s", "--segdup-bed",
        required=False,
        help="Segmental duplication regions (BED.gz) for masking",
    )
    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Output directory for results",
    )

    # Locus processing
    parser.add_argument(
        "--locus-padding",
        type=int,
        default=10000,
        help="Padding around locus boundaries (bp)",
    )
    parser.add_argument(
        "--segdup-threshold",
        type=float,
        default=0.5,
        help="Minimum segdup overlap fraction to mask a bin",
    )
    parser.add_argument(
        "--min-bins-per-locus",
        type=int,
        default=5,
        help="Minimum bins required to analyze a locus",
    )
    parser.add_argument(
        "--max-bins-per-interval",
        type=int,
        default=10,
        help="Maximum bins per interval after rebinning (0 = no rebinning)",
    )
    parser.add_argument(
        "--min-bins-per-call",
        type=int,
        default=3,
        help="Minimum bins required for a CNV call",
    )

    # CNV calling is now handled by plot_gd_cnv_output.py
    # Removed: --log-prob-threshold, --del-threshold, --dup-threshold, --min-confidence

    # Model parameters
    parser.add_argument(
        "--alpha-ref",
        type=float,
        default=1.0,
        help="Dirichlet concentration for reference CN state (CN=2)",
    )
    parser.add_argument(
        "--alpha-non-ref",
        type=float,
        default=1.0,
        help="Dirichlet concentration for non-reference CN states",
    )
    parser.add_argument(
        "--var-bias-bin",
        type=float,
        default=0.01,
        help="Variance for per-bin mean bias",
    )
    parser.add_argument(
        "--var-sample",
        type=float,
        default=0.001,
        help="Variance for per-sample variance factor",
    )
    parser.add_argument(
        "--var-bin",
        type=float,
        default=0.001,
        help="Variance for per-bin variance factor",
    )
    parser.add_argument(
        "--guide-type",
        type=str,
        default="delta",
        choices=["delta", "diagonal"],
        help="Type of variational guide",
    )
    parser.add_argument(
        "--clamp-threshold",
        type=float,
        default=5.0,
        help="Maximum value for depth clamping",
    )

    # Training parameters
    parser.add_argument(
        "--max-iter",
        type=int,
        default=2000,
        help="Maximum training iterations per locus",
    )
    parser.add_argument(
        "--lr-init",
        type=float,
        default=0.02,
        help="Initial learning rate",
    )
    parser.add_argument(
        "--lr-min",
        type=float,
        default=0.01,
        help="Minimum learning rate",
    )
    parser.add_argument(
        "--lr-decay",
        type=float,
        default=500,
        help="Learning rate decay constant",
    )
    parser.add_argument(
        "--log-freq",
        type=int,
        default=100,
        help="Logging frequency (iterations)",
    )
    parser.add_argument(
        "--jit",
        action="store_true",
        default=False,
        help="Enable JIT compilation",
    )
    parser.add_argument(
        "--early-stopping",
        action="store_true",
        default=True,
        help="Enable early stopping",
    )
    parser.add_argument(
        "--no-early-stopping",
        action="store_false",
        dest="early_stopping",
        help="Disable early stopping",
    )
    parser.add_argument(
        "--patience",
        type=int,
        default=50,
        help="Early stopping patience",
    )
    parser.add_argument(
        "--min-delta",
        type=float,
        default=100.0,
        help="Minimum improvement for early stopping",
    )

    # Inference parameters
    parser.add_argument(
        "--n-discrete-samples",
        type=int,
        default=500,
        help="Number of samples for discrete inference",
    )

    # Data filtering
    parser.add_argument(
        "--skip-bin-filter",
        action="store_true",
        default=False,
        help="Skip bin quality filtering",
    )
    parser.add_argument(
        "--median-min",
        type=float,
        default=1.0,
        help="Minimum median depth for bins",
    )
    parser.add_argument(
        "--median-max",
        type=float,
        default=3.0,
        help="Maximum median depth for bins",
    )
    parser.add_argument(
        "--mad-max",
        type=float,
        default=2.0,
        help="Maximum MAD for bins",
    )

    # Device
    parser.add_argument(
        "--device",
        type=str,
        default="cpu",
        choices=["cpu", "cuda"],
        help="Device for computation",
    )

    return parser.parse_args()


def main():
    """Main function to run GD CNV detection pipeline."""
    args = parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Output directory: {args.output_dir}")

    # Load GD table
    print(f"\nLoading GD table: {args.gd_table}")
    gd_table = GDTable(args.gd_table)
    print(f"Loaded {len(gd_table.loci)} loci")
    for cluster, locus in gd_table.loci.items():
        # Get overall start and end from breakpoint ranges
        if locus.breakpoints:
            overall_start = min(bp[0] for bp in locus.breakpoints)
            overall_end = max(bp[1] for bp in locus.breakpoints)
            print(f"  {cluster}: {locus.chrom}:{overall_start}-{overall_end} "
                  f"({len(locus.gd_entries)} entries, {locus.n_breakpoints} breakpoints)")
        else:
            print(f"  {cluster}: {locus.chrom} - NO BREAKPOINTS DEFINED")

    # Load segdup mask if provided
    segdup_mask = None
    if args.segdup_bed:
        print(f"\nLoading segdup mask: {args.segdup_bed}")
        segdup_mask = SegDupMask(args.segdup_bed)

    # Load read depth data
    df = read_data(args.input)

    # Normalize by sample median over autosomal bins
    sample_cols = get_sample_columns(df)
    autosome_mask = ~df["Chr"].isin(["chrX", "chrY"])
    if autosome_mask.any():
        column_medians = np.median(df.loc[autosome_mask, sample_cols], axis=0)
    else:
        column_medians = np.median(df[sample_cols], axis=0)

    print(f"Column medians: min={column_medians.min():.3f}, "
          f"max={column_medians.max():.3f}, mean={column_medians.mean():.3f}")

    # Normalize such that CN=2 corresponds to depth of 2.0
    df[sample_cols] = 2.0 * df[sample_cols] / column_medians[np.newaxis, :]

    # Filter low quality bins
    if not args.skip_bin_filter:
        df = filter_low_quality_bins(
            df,
            median_min=args.median_min,
            median_max=args.median_max,
            mad_max=args.mad_max,
        )

    # Set up Pyro
    pyro.enable_validation(True)
    pyro.distributions.enable_validation(True)
    pyro.set_rng_seed(42)
    torch.manual_seed(42)
    np.random.seed(42)

    # Run GD analysis (training and inference only)
    run_gd_analysis(
        df, gd_table, segdup_mask, args, device=args.device
    )

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print("\nOutput files:")
    print(f"  - cn_posteriors.tsv.gz")
    print(f"  - sample_posteriors.tsv.gz")
    print(f"  - bin_posteriors.tsv.gz")
    print(f"  - bin_mappings.tsv.gz")
    print(f"  - locus_intervals.tsv.gz")
    print("\nNext step: Run plot_gd_cnv_output.py to call CNVs and generate plots.")
    print("=" * 80)

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
