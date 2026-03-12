"""
Depth data containers and Bayesian CNV model.

Contains:
    - ExclusionMask: genomic exclusion region handling
    - DepthData: depth matrix container for the Pyro model
    - CNVModel: hierarchical Bayesian CNV detection model
"""

from typing import Dict, List, Tuple, Union

import numpy as np
import pandas as pd
import torch
from intervaltree import IntervalTree
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


def build_diploid_pair_states(max_hap_cn: int = 2) -> List[Tuple[int, int]]:
    """Return canonical unordered diploid pair states with h1 <= h2."""
    return [
        (h1, h2)
        for h1 in range(max_hap_cn + 1)
        for h2 in range(h1, max_hap_cn + 1)
    ]


def pair_state_minor_baf(pair_states: List[Tuple[int, int]]) -> np.ndarray:
    """Expected minor-allele BAF for each diploid pair state."""
    values = []
    for h1, h2 in pair_states:
        total = h1 + h2
        if total <= 0:
            values.append(0.0)
        else:
            values.append(min(h1, h2) / total)
    return np.asarray(values, dtype=np.float32)


def pair_state_total_cn(pair_states: List[Tuple[int, int]]) -> np.ndarray:
    """Total CN implied by each diploid pair state."""
    return np.asarray([h1 + h2 for h1, h2 in pair_states], dtype=np.int64)


class ExclusionMask:
    """
    Handles genomic exclusion regions for masking unreliable depth signal.

    Uses :class:`intervaltree.IntervalTree` per chromosome for exact
    interval-overlap queries (no midpoint heuristics).

    Exclusion regions can include segmental duplications, centromeres,
    satellite repeats, or any other intervals the user wishes to mask.
    Multiple BED files are accepted; all intervals are concatenated and
    merged *once* before the index is built, ensuring cross-file overlaps
    are collapsed properly.

    Any BED file with at least three columns (chr, start, end) is accepted.
    Additional columns are ignored.
    """

    def __init__(self, filepaths: Union[str, List[str]], label: str = "exclusion regions"):
        """
        Load genomic regions from one or more BED files.

        All files are read and concatenated first, then overlapping
        intervals are merged across the combined set before the lookup
        index is built.  This guarantees that cross-file overlaps are
        collapsed correctly and the index is only constructed once.

        Only the first three columns (chr, start, end) are required.  Extra
        BED columns (name, score, strand, …) are automatically ignored.

        Args:
            filepaths: Path to a BED file (plain or bgzipped .bed.gz),
                or a list of such paths.  All files are merged before
                the index is built.
            label: Human-readable label used in log messages.
        """
        self.label = label
        if isinstance(filepaths, str):
            filepaths = [filepaths]
        dfs = []
        for fp in filepaths:
            df = pd.read_csv(
                fp,
                sep="\t",
                header=None,
                usecols=[0, 1, 2],
                names=["chr", "start", "end"],
                compression="gzip" if fp.endswith(".gz") else None,
            )
            print(f"  Loaded {len(df):,} regions from {fp}")
            dfs.append(df)
        self.df = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame(
            columns=["chr", "start", "end"]
        )
        self._build_interval_index()
        print(f"Exclusion mask: {len(self.df):,} total {label} regions "
              f"({len(filepaths)} file(s)) → "
              f"{sum(len(t) for t in self.trees.values())} merged intervals")

    def _build_interval_index(self):
        """Build per-chromosome :class:`IntervalTree` instances.

        Intervals are merged (via ``tree.merge_overlaps()``) so that
        overlap-fraction queries can sum non-overlapping pieces directly.
        """
        self.trees: Dict[str, IntervalTree] = {}
        for chrom, group in self.df.groupby("chr"):
            tree = IntervalTree()
            for s, e in zip(group["start"].values, group["end"].values):
                s, e = int(s), int(e)
                if e > s:
                    tree.addi(s, e)
            tree.merge_overlaps()
            self.trees[chrom] = tree

    def has_any_overlap(self, chrom: str, start: int, end: int) -> bool:
        """Return *True* if ``[start, end)`` overlaps any exclusion interval.

        This is a fast O(log N + k) query with no fraction computation.
        """
        if chrom not in self.trees:
            return False
        return bool(self.trees[chrom].overlap(start, end))

    def get_overlap_fraction(self, chrom: str, start: int, end: int) -> float:
        """
        Calculate the fraction of ``[start, end)`` overlapping exclusion intervals.

        The tree is pre-merged so every hit is non-overlapping; the total
        overlap is the sum of per-interval intersections.

        Args:
            chrom: Chromosome name
            start: Region start position
            end: Region end position

        Returns:
            Fraction of region overlapping with exclusion intervals (0.0 to 1.0)
        """
        region_length = end - start
        if region_length <= 0 or chrom not in self.trees:
            return 0.0
        hits = self.trees[chrom].overlap(start, end)
        if not hits:
            return 0.0
        total = sum(min(iv.end, end) - max(iv.begin, start) for iv in hits)
        return min(total / region_length, 1.0)

    def get_overlap_fractions_batch(
        self, chrom: str, starts: np.ndarray, ends: np.ndarray
    ) -> np.ndarray:
        """
        Compute overlap fractions for multiple bins at once.

        Args:
            chrom: Chromosome name
            starts: Array of bin start positions
            ends: Array of bin end positions

        Returns:
            Array of overlap fractions (0.0 to 1.0) for each bin
        """
        n = len(starts)
        if n == 0 or chrom not in self.trees:
            return np.zeros(n)
        tree = self.trees[chrom]
        starts = np.asarray(starts)
        ends = np.asarray(ends)
        result = np.zeros(n)
        for i in range(n):
            s, e = int(starts[i]), int(ends[i])
            length = e - s
            if length <= 0:
                continue
            hits = tree.overlap(s, e)
            if hits:
                total = sum(min(iv.end, e) - max(iv.begin, s) for iv in hits)
                result[i] = min(total / length, 1.0)
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

    def add_beds(self, filepaths: Union[str, List[str]]) -> None:
        """Add intervals from one or more additional BED files.

        All new intervals are appended to the existing set and the
        internal index is rebuilt once so that cross-file overlaps are
        merged before any query is made.

        Args:
            filepaths: A single BED file path or a list of paths.
        """
        if isinstance(filepaths, str):
            filepaths = [filepaths]
        dfs = [self.df]
        for fp in filepaths:
            df = pd.read_csv(
                fp,
                sep="\t",
                header=None,
                usecols=[0, 1, 2],
                names=["chr", "start", "end"],
                compression="gzip" if fp.endswith(".gz") else None,
            )
            print(f"  Adding {len(df):,} regions from {fp}")
            dfs.append(df)
        self.df = pd.concat(dfs, ignore_index=True)
        self._build_interval_index()
        print(f"Exclusion mask updated: {len(self.df):,} total regions → "
              f"{sum(len(t) for t in self.trees.values())} merged intervals")


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

        # Compute interval sizes (bp) for each bin
        interval_sizes = (self.end - self.start).astype(float)
        self.interval_sizes = torch.tensor(
            interval_sizes, dtype=dtype, device=device
        ).unsqueeze(-1)  # (n_bins, 1) for broadcasting over samples

        print(f"Loaded data: {self.n_bins} bins x {self.n_samples} samples")
        print(f"Depth range: [{self.depth.min():.3f}, {self.depth.max():.3f}]")
        print(f"Depth mean: {self.depth.mean():.3f}, median: {self.depth.median():.3f}")
        print(f"Interval sizes: [{interval_sizes.min():.0f}, {interval_sizes.max():.0f}] bp, "
              f"median={np.median(interval_sizes):.0f} bp")

        # Optional SNP/BAF summaries can be attached later once the
        # preprocess bin mappings are available.
        self.baf_median = None
        self.minor_baf_median = None
        self.baf_variance = None
        self.baf_n_sites = None
        self.has_baf = False

    def attach_baf_summary(self, baf_summary_df: pd.DataFrame, mappings) -> None:
        """Attach per-bin, per-sample BAF summaries to this data object.

        The summary table is expected to contain rows keyed by preprocess
        ``array_idx`` and sample identifier, with columns such as
        ``baf_median``, ``minor_baf_median``, ``baf_variance``, and
        ``baf_n_sites``.  Only exact sample-id matches are attached.
        """
        if baf_summary_df is None or len(baf_summary_df) == 0:
            return

        sample_to_idx = {str(sample_id): idx for idx, sample_id in enumerate(self.sample_ids)}
        n_bins = self.n_bins
        n_samples = self.n_samples

        baf_median = np.full((n_bins, n_samples), np.nan, dtype=np.float32)
        minor_baf_median = np.full((n_bins, n_samples), np.nan, dtype=np.float32)
        baf_variance = np.full((n_bins, n_samples), np.nan, dtype=np.float32)
        baf_n_sites = np.zeros((n_bins, n_samples), dtype=np.int32)

        n_attached = 0
        for row in baf_summary_df.itertuples(index=False):
            try:
                bin_idx = int(row.array_idx)
            except (TypeError, ValueError):
                continue
            sample_idx = sample_to_idx.get(str(row.sample))
            if sample_idx is None or not (0 <= bin_idx < n_bins):
                continue

            baf_median[bin_idx, sample_idx] = float(row.baf_median)
            minor_baf_median[bin_idx, sample_idx] = float(row.minor_baf_median)
            baf_variance[bin_idx, sample_idx] = float(row.baf_variance)
            baf_n_sites[bin_idx, sample_idx] = int(row.baf_n_sites)
            n_attached += 1

        self.baf_median = torch.tensor(baf_median, dtype=self.depth.dtype, device=self.depth.device)
        self.minor_baf_median = torch.tensor(minor_baf_median, dtype=self.depth.dtype, device=self.depth.device)
        self.baf_variance = torch.tensor(baf_variance, dtype=self.depth.dtype, device=self.depth.device)
        self.baf_n_sites = torch.tensor(baf_n_sites, dtype=torch.int32, device=self.depth.device)
        self.has_baf = n_attached > 0

        print(f"Attached BAF summaries: {n_attached:,} matched bin × sample rows")


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
        state_prior_weight: float = 1.0,
        var_bias_bin: float = 0.1,
        var_sample: float = 0.2,
        var_bin: float = 0.2,
        bin_size_factor: float = 10000.0,
        device: str = "cpu",
        dtype: torch.dtype = torch.float32,
        debug: bool = False,
        guide_type: str = "diagonal",
    ):
        """
        Args:
            n_states: Number of latent states. Default 6 corresponds to the
                canonical unordered diploid pair states over per-haplotype CN
                values [0, 1, 2]: (0,0), (0,1), (0,2), (1,1), (1,2), (2,2).
            alpha_ref: Dirichlet concentration for reference state ((1,1))
            alpha_non_ref: Dirichlet concentration for non-reference states
            state_prior_weight: Multiplicative weight applied to the learned
                per-bin pair-state log-prior when reconstructing analytical
                discrete posteriors. Values < 1.0 temper overly sharp priors.
            var_bias_bin: Variance for per-bin mean bias (log-normal)
            var_sample: Variance for per-sample variance factor (log-normal)
            var_bin: Variance for per-bin variance factor (log-normal)
            bin_size_factor: Reference bin size (bp) for variance scaling.
                The total variance is multiplied by
                ``bin_size_factor / interval_size`` so that smaller bins
                have proportionally higher variance.
            device: torch device
            dtype: torch data type
            debug: Whether to print debug statements in model()
            guide_type: Type of variational guide ('diagonal' for AutoDiagonalNormal, 'delta' for AutoDelta)
        """
        self.n_states = n_states
        self.alpha_ref = alpha_ref
        self.alpha_non_ref = alpha_non_ref
        self.state_prior_weight = state_prior_weight
        self.var_bias_bin = var_bias_bin
        self.var_sample = var_sample
        self.var_bin = var_bin
        self.bin_size_factor = bin_size_factor
        self.device = device
        self.dtype = dtype
        self.debug = debug
        self.guide_type = guide_type
        self.pair_states = build_diploid_pair_states(max_hap_cn=2)
        if n_states != len(self.pair_states):
            raise ValueError(
                f"n_states={n_states} does not match diploid pair-state count "
                f"{len(self.pair_states)}"
            )
        self.total_cn_by_state = torch.tensor(
            pair_state_total_cn(self.pair_states),
            device=self.device,
            dtype=self.dtype,
        )
        self.minor_baf_by_state = torch.tensor(
            pair_state_minor_baf(self.pair_states),
            device=self.device,
            dtype=self.dtype,
        )
        self.ref_state_idx = self.pair_states.index((1, 1))
        self.max_total_cn = int(max(sum(p) for p in self.pair_states))

        # Training history
        self.loss_history = {"epoch": [], "elbo": []}

        # Define which sites to expose to the guide (continuous latent variables)
        self.latent_sites = ["bin_bias", "sample_var", "bin_var", "pair_state_probs"]

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
    def model(self, depth: torch.Tensor, interval_sizes: torch.Tensor, n_bins: int = None, n_samples: int = None):
        """
        Probabilistic model for CNV detection.

        Args:
            depth: Observed normalized read depth (n_bins x n_samples)
            interval_sizes: Bin sizes in bp (n_bins x 1) for variance scaling
        """

        if self.debug:
            print("\n=== MODEL DEBUG ===")
            print(f"depth.shape: {depth.shape}")

        zero_t = torch.zeros(1, device=self.device, dtype=self.dtype)
        one_t = torch.ones(1, device=self.device, dtype=self.dtype)

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

            # Pair-state prior (Dirichlet-Categorical)
            # Heavily weight the diploid reference state (1,1)
            alpha_pair = self.alpha_non_ref * one_t.expand(self.n_states)
            alpha_pair = alpha_pair.clone()
            alpha_pair[self.ref_state_idx] = torch.tensor(
                self.alpha_ref, device=self.device, dtype=self.dtype,
            )
            pair_state_probs = pyro.sample("pair_state_probs", dist.Dirichlet(alpha_pair))
        if self.debug:
            print(f"pair_state_probs.shape: {pair_state_probs.shape}")

        # Per-bin, per-sample pair state and observation
        with plate_bins, plate_samples:
            # Sample pair state (discrete latent variable)
            # Shape: (n_bins, n_samples)
            pair_state = pyro.sample("pair_state", dist.Categorical(pair_state_probs))
            if self.debug:
                print(f"pair_state.shape: {pair_state.shape}")

            # Expected depth depends on total CN implied by the pair state.
            if self.debug:
                print(f"bin_bias.shape: {bin_bias.shape}")
            expected_total_cn = Vindex(self.total_cn_by_state)[pair_state]
            expected_depth = expected_total_cn * bin_bias
            if self.debug:
                print(f"expected_depth.shape: {expected_depth.shape}")

            # Variance: combination of sample and bin variance,
            # scaled inversely by bin size so that smaller bins have
            # proportionally higher variance.
            if self.debug:
                print(f"bin_var).shape: {bin_var.shape}")
            base_variance = sample_var + bin_var
            size_modifier = self.bin_size_factor / interval_sizes
            variance = base_variance * size_modifier
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

            # Optional BAF observation on the minor-allele fraction.
            # The observed variance is estimated upstream from the number of
            # SNP sites per bin, so bins with more sites contribute a tighter
            # likelihood. Missing / unsupported bins are masked out.
            if hasattr(self, "current_data") and getattr(self.current_data, "has_baf", False):
                baf_obs = self.current_data.minor_baf_median
                baf_var = self.current_data.baf_variance
                baf_sites = self.current_data.baf_n_sites

                valid_mask = (torch.isfinite(baf_obs) & torch.isfinite(baf_var) & (baf_sites > 0) & (baf_var > 0))
                if torch.any(valid_mask):
                    expected_minor_baf = Vindex(self.minor_baf_by_state)[pair_state]
                    safe_baf_var = torch.where(valid_mask, baf_var, torch.ones_like(baf_var))
                    baf_std = torch.sqrt(torch.clamp(safe_baf_var, min=1e-6))
                    safe_baf_obs = torch.where(valid_mask, baf_obs, expected_minor_baf)
                    with poutine.mask(mask=valid_mask):
                        pyro.sample(
                            "baf_obs",
                            dist.Normal(expected_minor_baf, baf_std),
                            obs=safe_baf_obs,
                        )
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
            self.current_data = data
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
                        depth=data.depth,
                        interval_sizes=data.interval_sizes,
                        n_bins=data.n_bins,
                        n_samples=data.n_samples,
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
        finally:
            self.current_data = None

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
            depth=data.depth,
            interval_sizes=data.interval_sizes,
            n_bins=data.n_bins,
            n_samples=data.n_samples,
        )

        # Get model trace conditioned on guide
        trained_model = poutine.replay(self.model, trace=guide_trace)

        # Get discrete MAP using infer_discrete
        inferred_model = infer_discrete(
            trained_model, temperature=0, first_available_dim=-3
        )
        trace = poutine.trace(inferred_model).get_trace(
            depth=data.depth,
            interval_sizes=data.interval_sizes,
            n_bins=data.n_bins,
            n_samples=data.n_samples,
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
        map_estimates["pair_state_probs"] = (
            guide_trace.nodes["pair_state_probs"]["value"].detach().cpu().numpy()
        )

        # Discrete variable (pair state), plus compatibility total CN.
        pair_state_idx = trace.nodes["pair_state"]["value"].detach().cpu().numpy().squeeze()
        map_estimates["pair_state"] = pair_state_idx
        total_cn_lookup = pair_state_total_cn(self.pair_states)
        map_estimates["cn"] = total_cn_lookup[pair_state_idx]
        map_estimates["pair_state_labels"] = self.pair_states

        pair_state_probs = np.asarray(map_estimates["pair_state_probs"]).squeeze()
        if pair_state_probs.ndim == 1:
            pair_state_probs = pair_state_probs.reshape(1, -1)
        cn_probs = np.zeros((pair_state_probs.shape[0], self.max_total_cn + 1), dtype=np.float32)
        for pair_idx, total_cn in enumerate(total_cn_lookup):
            cn_probs[:, total_cn] += pair_state_probs[:, pair_idx]
        map_estimates["cn_probs"] = cn_probs

        return map_estimates

    def run_discrete_inference(self, data, **kwargs):
        """
        Compute exact posterior probabilities analytically using Bayes' rule.

        Because the model uses a discrete Categorical hidden state with a
        small state space (K=6) and a Gaussian observation likelihood, the
        posterior can be computed exactly once the MAP estimates for the
        continuous latent variables are available.  For each bin *b* and
        sample *s*:

        .. math::

            P(CN_{b,s} = k \\mid x_{b,s})
            = \\frac{P(x_{b,s} \\mid CN_{b,s}=k) \\cdot P(CN_b=k)}
                    {\\sum_{j} P(x_{b,s} \\mid CN_{b,s}=j) \\cdot P(CN_b=j)}

        This replaces the previous Monte Carlo sampling approach
        (``infer_discrete`` with thousands of samples) and is both faster
        and perfectly accurate.

        Args:
            data: DepthData object.
            **kwargs: Accepted for backward compatibility (``n_samples``,
                ``log_freq``) but ignored.

        Returns:
            Dictionary containing pair-state posteriors and total-CN
            marginals.
        """
        print("Computing exact discrete marginal posteriors analytically...")

        # 1. Get MAP estimates of continuous variables
        maps = self.get_map_estimates(data)
        bin_bias = maps["bin_bias"]       # (n_bins,) or (n_bins, 1)
        sample_var = maps["sample_var"]   # (n_samples,) or (1, n_samples)
        bin_var = maps["bin_var"]         # (n_bins,) or (n_bins, 1)
        pair_state_probs = maps["pair_state_probs"]

        # Flatten to 1-D / 2-D where needed (guide may add singleton
        # plate dimensions).
        bin_bias = bin_bias.squeeze()
        sample_var = sample_var.squeeze()
        bin_var = bin_var.squeeze()
        pair_state_probs = pair_state_probs.squeeze()   # (n_bins, n_pair_states)

        # 2. Prepare data matrices
        obs = data.depth.detach().cpu().numpy()  # (n_bins, n_samples)
        interval_sizes = data.interval_sizes.detach().cpu().numpy().squeeze()  # (n_bins,)
        pair_total_cn = pair_state_total_cn(self.pair_states)
        pair_minor_baf = pair_state_minor_baf(self.pair_states)

        # Combined variance: (sample_var + bin_var) * bin_size_factor / interval_size
        # → shape (n_bins, n_samples)
        base_variance = sample_var[np.newaxis, :] + bin_var[:, np.newaxis]
        size_modifier = self.bin_size_factor / interval_sizes[:, np.newaxis]
        variance = base_variance * size_modifier
        std = np.sqrt(variance)

        # 3. Compute log-likelihoods for all states simultaneously
        # states: (n_states, 1, 1);  bin_bias: (1, n_bins, 1)
        states_total_cn = pair_total_cn.reshape(-1, 1, 1)
        expected_depth = states_total_cn * bin_bias[np.newaxis, :, np.newaxis]

        # Broadcast obs and std to (1, n_bins, n_samples)
        obs_b = obs[np.newaxis, :, :]
        std_b = std[np.newaxis, :, :]

        # Gaussian log-PDF
        log_lik = -0.5 * np.log(2 * np.pi * std_b ** 2) - (
            (obs_b - expected_depth) ** 2
        ) / (2 * std_b ** 2)

        # 4. Add log-prior over pair states.
        log_prior = np.log(np.maximum(pair_state_probs.T[:, :, np.newaxis], 1e-10))
        if self.state_prior_weight != 1.0:
            log_prior = self.state_prior_weight * log_prior
        log_unnormalized = log_lik + log_prior

        # 4b. Optional BAF log-likelihood contribution.
        if getattr(data, "has_baf", False):
            minor_baf = data.minor_baf_median.detach().cpu().numpy()
            baf_var = data.baf_variance.detach().cpu().numpy()
            baf_sites = data.baf_n_sites.detach().cpu().numpy()

            valid = (np.isfinite(minor_baf) & np.isfinite(baf_var) & (baf_sites > 0) &
                     (baf_var > 0))
            if np.any(valid):
                exp_minor_baf = pair_minor_baf.reshape(-1, 1, 1)
                baf_obs = minor_baf[np.newaxis, :, :]
                baf_std = np.sqrt(np.maximum(baf_var[np.newaxis, :, :], 1e-6))
                baf_log_lik = -0.5 * np.log(2 * np.pi * baf_std ** 2) - (
                    (baf_obs - exp_minor_baf) ** 2
                ) / (2 * baf_std ** 2)
                baf_log_lik = np.where(valid[np.newaxis, :, :], baf_log_lik, 0.0)
                log_unnormalized += baf_log_lik

        # 5. Log-sum-exp softmax across the state dimension (axis 0)
        max_log = np.max(log_unnormalized, axis=0, keepdims=True)
        exp_vals = np.exp(log_unnormalized - max_log)
        posterior = exp_vals / np.sum(exp_vals, axis=0, keepdims=True)

        # Transpose from (n_states, n_bins, n_samples) → (n_bins, n_samples, n_states)
        pair_posterior = np.transpose(posterior, (1, 2, 0))

        cn_posterior = np.zeros(
            (pair_posterior.shape[0], pair_posterior.shape[1], self.max_total_cn + 1),
            dtype=np.float32,
        )
        for pair_idx, total_cn in enumerate(pair_total_cn):
            cn_posterior[:, :, total_cn] += pair_posterior[:, :, pair_idx]

        print("Exact analytical inference complete.")
        return {
            "cn_posterior": cn_posterior,
            "pair_state_posterior": pair_posterior,
            "pair_state_labels": self.pair_states,
        }


