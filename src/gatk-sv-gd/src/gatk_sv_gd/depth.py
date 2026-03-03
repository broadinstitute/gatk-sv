"""
Depth data containers and Bayesian CNV model.

Contains:
    - ExclusionMask: genomic exclusion region handling
    - DepthData: depth matrix container for the Pyro model
    - CNVModel: hierarchical Bayesian CNV detection model
"""

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm

import pyro
import pyro.distributions as dist
from pyro import poutine
from pyro.ops.indexing import Vindex
from pyro.infer import config_enumerate, infer_discrete
from pyro.infer.autoguide import AutoDiagonalNormal, AutoDelta
from pyro.infer import JitTraceEnum_ELBO, TraceEnum_ELBO
from pyro.infer.svi import SVI

from gatk_sv_gd._util import get_sample_columns

class ExclusionMask:
    """
    Handles genomic exclusion regions for masking unreliable depth signal.

    Exclusion regions can include segmental duplications, centromeres,
    satellite repeats, or any other intervals the user wishes to mask.
    Multiple BED files can be loaded and merged via :meth:`merge`.

    Any BED file with at least three columns (chr, start, end) is accepted.
    Additional columns are ignored.
    """

    def __init__(self, filepath: str, label: str = "exclusion regions"):
        """
        Load genomic regions from a BED file.

        Only the first three columns (chr, start, end) are required.  Extra
        BED columns (name, score, strand, …) are automatically ignored.

        Args:
            filepath: Path to BED file (plain or bgzipped .bed.gz).
            label: Human-readable label used in log messages, e.g.
                the basename of the input file.
        """
        self.label = label
        self.df = pd.read_csv(
            filepath,
            sep="\t",
            header=None,
            usecols=[0, 1, 2],
            names=["chr", "start", "end"],
            compression="gzip" if filepath.endswith(".gz") else None,
        )
        self._build_interval_index()
        print(f"Loaded {len(self.df)} {label} regions")

    def _build_interval_index(self):
        """Build efficient interval index for overlap queries."""
        # Group by chromosome for efficient lookup
        self.intervals_by_chrom = {}
        self.merged_by_chrom = {}  # Pre-merged non-overlapping intervals
        for chrom, group in self.df.groupby("chr"):
            starts = group["start"].values
            ends = group["end"].values
            self.intervals_by_chrom[chrom] = {"starts": starts, "ends": ends}
            # Pre-merge overlapping intervals for O(log N) binary-search queries
            idx = np.argsort(starts)
            s, e = starts[idx], ends[idx]
            ms, me = [int(s[0])], [int(e[0])]
            for i in range(1, len(s)):
                if int(s[i]) <= me[-1]:
                    me[-1] = max(me[-1], int(e[i]))
                else:
                    ms.append(int(s[i]))
                    me.append(int(e[i]))
            self.merged_by_chrom[chrom] = {
                "starts": np.array(ms),
                "ends": np.array(me),
            }

    def get_overlap_fraction(self, chrom: str, start: int, end: int) -> float:
        """
        Calculate the fraction of a region that overlaps with exclusion intervals.

        Args:
            chrom: Chromosome name
            start: Region start position
            end: Region end position

        Returns:
            Fraction of region overlapping with exclusion intervals (0.0 to 1.0)
        """
        if chrom not in self.merged_by_chrom:
            return 0.0
        region_length = end - start
        if region_length <= 0:
            return 0.0
        merged = self.merged_by_chrom[chrom]
        ms, me = merged["starts"], merged["ends"]
        # Binary search: find merged intervals overlapping [start, end)
        lo = int(np.searchsorted(me, start, side="right"))
        hi = int(np.searchsorted(ms, end, side="left"))
        if lo >= hi:
            return 0.0
        ov_s = np.maximum(ms[lo:hi], start)
        ov_e = np.minimum(me[lo:hi], end)
        return min(float(np.sum(np.maximum(0, ov_e - ov_s))) / region_length, 1.0)

    def get_overlap_fractions_batch(
        self, chrom: str, starts: np.ndarray, ends: np.ndarray
    ) -> np.ndarray:
        """
        Compute overlap fractions for multiple bins at once.

        Uses vectorized searchsorted to identify bins that overlap any
        exclusion interval, then only iterates over those bins (typically
        a small fraction).

        Args:
            chrom: Chromosome name
            starts: Array of bin start positions
            ends: Array of bin end positions

        Returns:
            Array of overlap fractions (0.0 to 1.0) for each bin
        """
        n = len(starts)
        if n == 0 or chrom not in self.merged_by_chrom:
            return np.zeros(n)
        merged = self.merged_by_chrom[chrom]
        ms, me = merged["starts"], merged["ends"]
        starts = np.asarray(starts)
        ends = np.asarray(ends)
        # Vectorized: find lo/hi bracket for each bin in one searchsorted call
        lo_arr = np.searchsorted(me, starts, side="right")
        hi_arr = np.searchsorted(ms, ends, side="left")
        has_overlap = lo_arr < hi_arr
        result = np.zeros(n)
        # Only iterate over bins that actually overlap an exclusion interval
        for i in np.where(has_overlap)[0]:
            s, e = starts[i], ends[i]
            lo, hi = int(lo_arr[i]), int(hi_arr[i])
            ov_s = np.maximum(ms[lo:hi], s)
            ov_e = np.minimum(me[lo:hi], e)
            result[i] = min(float(np.sum(np.maximum(0, ov_e - ov_s))) / (e - s), 1.0)
        return result

    def is_masked(self, chrom: str, start: int, end: int,
                  threshold: float = 0.5) -> bool:
        """
        Check if a region should be masked due to exclusion overlap.

        Args:
            chrom: Chromosome name
            start: Region start position
            end: Region end position
            threshold: Minimum overlap fraction to trigger masking

        Returns:
            True if region should be masked
        """
        return self.get_overlap_fraction(chrom, start, end) >= threshold

    def merge(self, other_bed: str, label: str = "additional regions") -> None:
        """Merge intervals from another BED file into this mask.

        The intervals are appended to the existing set and the internal
        index is rebuilt so that subsequent overlap queries cover both the
        original and the newly added regions.

        Args:
            other_bed: Path to a BED file (plain or bgzipped).
            label: Human-readable label for log messages describing the
                regions being merged (e.g. the basename of the BED file).
        """
        other_df = pd.read_csv(
            other_bed,
            sep="\t",
            header=None,
            usecols=[0, 1, 2],
            names=["chr", "start", "end"],
            compression="gzip" if other_bed.endswith(".gz") else None,
        )
        n_before = len(self.df)
        self.df = pd.concat([self.df, other_df], ignore_index=True)
        self._build_interval_index()
        print(f"Merged {len(other_df)} {label} regions into mask "
              f"({n_before} → {len(self.df)} total)")


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

        # NOTE: Do NOT sort the DataFrame here.  The caller
        # (collect_all_locus_bins) builds a parallel `mappings` list whose
        # order must stay aligned with the rows of this DataFrame.  Sorting
        # would silently rearrange depth values relative to those mappings,
        # causing every downstream table (cn_posteriors, bin_posteriors, etc.)
        # to attribute values to the wrong genomic bins.

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
                "sample_var", dist.Exponential(torch.tensor(1.0 / self.var_sample, device=self.device, dtype=self.dtype))
            )
        if self.debug:
            print(f"sample_var.shape: {sample_var.shape}")

        # Per-bin latent variables
        with plate_bins:
            # Per-bin mean bias factor (log-normal prior, centered at 1.0)
            bin_bias = pyro.sample(
                "bin_bias", dist.LogNormal(zero_t, torch.tensor(self.var_bias_bin, device=self.device, dtype=self.dtype))
            )
            if self.debug:
                print(f"bin_bias.shape: {bin_bias.shape}")

            # Per-bin variance factor (log-normal prior)
            bin_var = pyro.sample("bin_var", dist.Exponential(torch.tensor(1.0 / self.var_bin, device=self.device, dtype=self.dtype)))
            if self.debug:
                print(f"bin_var.shape: {bin_var.shape}")

            # Copy number prior (Dirichlet-Categorical)
            # Heavily weight CN=2 (diploid)
            alpha_cn = self.alpha_non_ref * one_t.expand(self.n_states)
            alpha_cn = alpha_cn.clone()
            alpha_cn[2] = torch.tensor(self.alpha_ref, device=self.device, dtype=self.dtype)
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

        # Accumulate counts directly instead of storing all samples, to avoid
        # building an (n_samples x n_bins x n_samples) array in memory.
        cn_counts = np.zeros((data.n_bins, data.n_samples, self.n_states), dtype=np.int32)

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
                    cn = trace.nodes["cn"]["value"].detach().cpu().numpy()  # (n_bins, n_samples)
                    # Accumulate one-hot counts in-place — no sample list needed
                    for state in range(self.n_states):
                        cn_counts[..., state] += (cn == state).astype(np.int32)

                    # Periodically free GPU cache to prevent fragmentation
                    if (i + 1) % 100 == 0 and self.device == "cuda":
                        torch.cuda.empty_cache()

        cn_freq = cn_counts / n_samples  # normalise to probabilities

        print("Discrete inference complete.")
        return {"cn_posterior": cn_freq}


