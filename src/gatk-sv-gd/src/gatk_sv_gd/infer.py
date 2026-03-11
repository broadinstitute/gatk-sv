"""
Genomic Disorder (GD) CNV Detection from Binned Read Counts

This script detects copy number variants at known genomic disorder loci using
a hierarchical Bayesian model. Each locus (cluster) is defined by one or more
breakpoints, and the script determines which breakpoint set best explains the
observed depth signal for each sample.

Input:
    - Binned read count file (TSV)
    - GD table with locus definitions (TSV)
    - Optional exclusion interval BED files for masking (e.g. segdups, centromeres)

Output:
    - Per-locus CNV calls with breakpoint assignments
    - Copy number posteriors for each sample at each locus
"""

import argparse
import os
from typing import Optional

import numpy as np
import pandas as pd
import pyro
import torch

from gatk_sv_gd import _util
from gatk_sv_gd._util import get_sample_columns, setup_logging
from gatk_sv_gd.bins import filter_low_quality_bins, read_data
from gatk_sv_gd.depth import CNVModel, DepthData, ExclusionMask
from gatk_sv_gd.models import GDTable
from gatk_sv_gd.output import estimate_ploidy, write_locus_metadata, write_posterior_tables
from gatk_sv_gd.preprocess import collect_all_locus_bins, load_preprocessed_data


def run_gd_analysis(
    df: pd.DataFrame,
    gd_table: GDTable,
    exclusion_mask: Optional[ExclusionMask],
    args: argparse.Namespace,
    device: str = "cpu",
    column_medians: Optional[np.ndarray] = None,
    lowres_median_bin_size: Optional[float] = None,
    preprocessed_bins: Optional[pd.DataFrame] = None,
    preprocessed_mappings=None,
):
    """
    Run GD CNV analysis on all loci using a single unified model.

    This function performs model training and inference only.
    CNV calling is handled by downstream scripts (plot_gd_cnv_output.py).

    When *preprocessed_bins* and *preprocessed_mappings* are supplied (from
    the ``preprocess`` subcommand), bin collection is skipped entirely and
    the provided data is used directly.

    Args:
        df: DataFrame with normalized read depth
        gd_table: GDTable with locus definitions
        exclusion_mask: Optional ExclusionMask for filtering
        args: Command line arguments
        device: Torch device
        column_medians: Per-sample autosomal median raw counts (before
            normalisation).  Needed when ``args.high_res_counts`` is set.
        lowres_median_bin_size: Median bin size (bp) of the low-res file.
            Needed when ``args.high_res_counts`` is set.
        preprocessed_bins: Optional combined DataFrame from the
            ``preprocess`` subcommand.  When set, *preprocessed_mappings*
            must also be provided.
        preprocessed_mappings: Optional list of LocusBinMapping objects
            from the ``preprocess`` subcommand.
    """
    if preprocessed_bins is not None and preprocessed_mappings is not None:
        # Use preprocessed data directly — skip bin collection
        combined_df = preprocessed_bins
        mappings = preprocessed_mappings
        included_loci = None  # already written by preprocess
    else:
        # Build quality-filter params dict.
        filter_params: dict = {
            "median_min": args.median_min,
            "median_max": args.median_max,
            "mad_max": args.mad_max,
        }
        highres_path: Optional[str] = getattr(args, "high_res_counts", None)

        # Collect all bins across all loci
        combined_df, mappings, included_loci = collect_all_locus_bins(
            df, gd_table, exclusion_mask,
            exclusion_threshold=args.exclusion_threshold,
            locus_padding=args.locus_padding,
            min_bins_per_interval=args.min_bins_per_interval,
            max_bins_per_interval=args.max_bins_per_interval,
            highres_counts_path=highres_path,
            column_medians=column_medians,
            lowres_median_bin_size=lowres_median_bin_size,
            filter_params=filter_params,
            exclusion_bypass_threshold=args.exclusion_bypass_threshold,
            min_rebin_coverage=args.min_rebin_coverage,
            min_flank_bases=args.min_flank_bases,
            min_flank_bins=args.min_flank_bins,
            min_flank_coverage=args.min_flank_coverage,
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
        bin_size_factor=args.bin_size_factor,
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
        jit=not args.disable_jit,
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

    # Write locus metadata for downstream calling/plotting (skip if
    # preprocessed data was loaded — preprocess already wrote these files).
    if included_loci is not None:
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
        required=False,
        help="Input TSV file with normalized read depth (bins x samples). "
             "Not required when --preprocessed-dir is set.",
    )
    parser.add_argument(
        "-g", "--gd-table",
        required=False,
        help="GD locus definition table (TSV). "
             "Not required when --preprocessed-dir is set.",
    )
    parser.add_argument(
        "-e", "--exclusion-intervals",
        nargs="*",
        default=[],
        help="One or more BED files (plain or .bed.gz) of genomic regions "
             "to mask (e.g. segmental duplications, centromeres, satellites)."
             "  Intervals from all files are merged into a single exclusion "
             "mask.  Only the first three columns (chr, start, end) are "
             "required; additional columns are ignored.  May be specified "
             "multiple times or as a space-separated list.",
    )
    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Output directory for results",
    )
    parser.add_argument(
        "--high-res-counts",
        required=False,
        help="Optional bgzipped, tabix-indexed high-resolution read count "
             "file (.tsv.gz + .tbi).  When provided, loci with any body "
             "interval below the target bin count (--min-bins-per-interval) "
             "are re-queried at this finer resolution before the hard check "
             "is enforced.  The file must have the same sample columns as "
             "the low-res input and contain raw (un-normalised) counts.",
    )
    parser.add_argument(
        "--preprocessed-dir",
        required=False,
        help="Directory produced by 'gatk-sv-gd preprocess'.  When set, "
             "bins and mappings are loaded from preprocessed_bins.tsv.gz "
             "and bin_mappings.tsv.gz instead of re-running preprocessing. "
             "The -i/--input, -g/--gd-table and -e/--exclusion-intervals "
             "flags are not required in this mode.",
    )

    # Locus processing
    parser.add_argument(
        "--locus-padding",
        type=int,
        default=10000,
        help="Padding around locus boundaries (bp)",
    )
    parser.add_argument(
        "--exclusion-threshold",
        type=float,
        default=0.5,
        help="Minimum overlap fraction with exclusion regions to mask a bin",
    )
    parser.add_argument(
        "--exclusion-bypass-threshold",
        type=float,
        default=0.6,
        help="If a body interval (region between adjacent breakpoints) is "
             "at least this fraction overlapped by exclusion regions, "
             "masking is skipped for bins in that interval.  This "
             "prevents intervals that are entirely within excluded "
             "regions from losing all their bins.  Set to 1.0 to "
             "disable bypass.",
    )
    parser.add_argument(
        "--min-bins-per-interval",
        type=int,
        default=10,
        help="Hard-failure minimum bins per body interval.  Intervals "
             "below this count after all processing cause a hard failure.",
    )
    parser.add_argument(
        "--max-bins-per-interval",
        type=int,
        default=20,
        help="Maximum bins per body interval after rebinning "
             "(0 = no rebinning)",
    )
    parser.add_argument(
        "--min-rebin-coverage",
        type=float,
        default=0.5,
        help="Minimum fraction of each new rebinned bin's width that must "
             "be covered by original bins (0.0–1.0, default 0.5).  Putative "
             "bins with less coverage are discarded to avoid biased estimates.",
    )
    parser.add_argument(
        "--min-flank-bases",
        type=int,
        default=50000,
        help="Minimum cumulative base pairs each flank must cover, "
             "regardless of locus size.  For small loci (e.g. < 1 kb) this "
             "ensures flanks are wide enough to establish a reliable "
             "baseline depth.  Flanks keep growing outward until this AND "
             "the locus-size target AND --min-flank-bins are all satisfied.",
    )
    parser.add_argument(
        "--min-flank-bins",
        type=int,
        default=10,
        help="Minimum number of bins each flank must contain, regardless "
             "of the base-pair thresholds.  Flanks keep growing outward "
             "until this AND the base-pair targets are all satisfied.",
    )
    parser.add_argument(
        "--min-flank-coverage",
        type=float,
        default=0.5,
        help="Minimum fraction of the effective bp target that a flank's "
             "accumulated bin coverage should reach.  The flank generator "
             "keeps extending outward until this threshold (and the other "
             "criteria) are met; it only stops when it runs out of bins "
             "at the chromosome boundary.  When the threshold cannot be "
             "met, a warning is logged but the flank is kept "
             "(default 0.5 = 50%%).",
    )

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
        "--bin-size-factor",
        type=float,
        default=10000.0,
        help="Reference bin size (bp) for variance scaling.  The total "
             "variance is multiplied by bin_size_factor / interval_size "
             "so that smaller bins have proportionally higher variance.  "
             "Note that this is redundant with the other scale factors "
             "and is only exposed for debugging. Set to 0 to disable "
             "bin-size variance scaling.",
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
        default=5000,
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
        "--disable-jit",
        action="store_true",
        default=False,
        help="Disable JIT compilation",
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
        default=1000.0,
        help="Minimum improvement for early stopping",
    )

    # Inference parameters
    parser.add_argument(
        "--n-discrete-samples",
        type=int,
        default=1000,
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

    # Verbosity
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        default=False,
        help="Enable detailed per-bin diagnostic logging for normalisation, "
             "quality filtering, exclusion masking, high-res replacement, and "
             "rebinning.  Useful for investigating discrepancies between "
             "raw and normalised depth signals.",
    )

    return parser.parse_args()


def main():
    """Main function to run GD CNV detection pipeline."""
    args = parse_args()

    _util.VERBOSE = args.verbose

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Mirror all stdout/stderr to a log file in the output directory
    log_fh = setup_logging(args.output_dir)

    print(f"Output directory: {args.output_dir}")

    # ------------------------------------------------------------------
    # Branch: load preprocessed data or run full preprocessing
    # ------------------------------------------------------------------
    if args.preprocessed_dir:
        print(f"\nLoading preprocessed data from: {args.preprocessed_dir}")
        preprocessed_bins, preprocessed_mappings = load_preprocessed_data(
            args.preprocessed_dir
        )

        # Set up Pyro
        pyro.enable_validation(True)
        pyro.distributions.enable_validation(True)
        pyro.set_rng_seed(42)
        torch.manual_seed(42)
        np.random.seed(42)

        # Run GD analysis with preprocessed data (no bin collection)
        run_gd_analysis(
            pd.DataFrame(), GDTable.__new__(GDTable), None, args,
            device=args.device,
            preprocessed_bins=preprocessed_bins,
            preprocessed_mappings=preprocessed_mappings,
        )
    else:
        # Validate that required args are present
        if not args.input:
            print("Error: --input is required unless --preprocessed-dir is set.")
            raise SystemExit(1)
        if not args.gd_table:
            print("Error: --gd-table is required unless --preprocessed-dir is set.")
            raise SystemExit(1)

        # Load GD table
        print(f"\nLoading GD table: {args.gd_table}")
        gd_table = GDTable(args.gd_table)
        print(f"Loaded {len(gd_table.loci)} loci")
        for cluster, locus in gd_table.loci.items():
            if locus.breakpoints:
                overall_start = min(bp[0] for bp in locus.breakpoints)
                overall_end = max(bp[1] for bp in locus.breakpoints)
                print(f"  {cluster}: {locus.chrom}:{overall_start}-{overall_end} "
                      f"({len(locus.gd_entries)} entries, {locus.n_breakpoints} breakpoints)")
            else:
                print(f"  {cluster}: {locus.chrom} - NO BREAKPOINTS DEFINED")

        # Load exclusion mask
        exclusion_mask = None
        for bed_path in args.exclusion_intervals:
            bed_label = os.path.basename(bed_path)
            if exclusion_mask is None:
                print(f"\nLoading exclusion mask: {bed_path}")
                exclusion_mask = ExclusionMask(bed_path, label=bed_label)
            else:
                print(f"\nMerging exclusion intervals: {bed_path}")
                exclusion_mask.merge(bed_path, label=bed_label)

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

        if _util.VERBOSE:
            print("\n  [verbose] Per-sample autosomal median raw counts:")
            for i, s in enumerate(sample_cols):
                print(f"    {s}: {column_medians[i]:.3f}")
            raw_depths = df[sample_cols].values
            print(f"\n  [verbose] Pre-normalisation depth summary "
                  f"({len(df)} bins x {len(sample_cols)} samples):")
            print(f"    global mean = {np.nanmean(raw_depths):.4f}")
            print(f"    global median = {np.nanmedian(raw_depths):.4f}")
            print(f"    per-sample means: min={np.nanmean(raw_depths, axis=0).min():.4f}, "
                  f"max={np.nanmean(raw_depths, axis=0).max():.4f}")

        lowres_bin_sizes = (df["End"] - df["Start"]).values
        lowres_median_bin_size = float(np.median(lowres_bin_sizes))
        print(f"Low-res median bin size: {lowres_median_bin_size:,.0f} bp")

        # Normalize such that CN=2 corresponds to depth of 2.0
        df[sample_cols] = 2.0 * df[sample_cols] / column_medians[np.newaxis, :]

        if _util.VERBOSE:
            norm_depths = df[sample_cols].values
            print(f"\n  [verbose] Post-normalisation depth summary:")
            print(f"    global mean = {np.nanmean(norm_depths):.4f}")
            print(f"    global median = {np.nanmedian(norm_depths):.4f}")
            print(f"    per-sample means: min={np.nanmean(norm_depths, axis=0).min():.4f}, "
                  f"max={np.nanmean(norm_depths, axis=0).max():.4f}")

        # Filter low quality bins
        if not args.skip_bin_filter:
            df = filter_low_quality_bins(
                df,
                median_min=args.median_min,
                median_max=args.median_max,
                mad_max=args.mad_max,
            )

        # Estimate ploidy
        ploidy_df = estimate_ploidy(df, args.output_dir)

        # Set up Pyro
        pyro.enable_validation(True)
        pyro.distributions.enable_validation(True)
        pyro.set_rng_seed(42)
        torch.manual_seed(42)
        np.random.seed(42)

        # Log high-res counts file if provided
        if args.high_res_counts:
            print(f"\nHigh-resolution counts file: {args.high_res_counts}")
        else:
            print("\nNo high-resolution counts file provided (--high-res-counts)")

        # Run GD analysis (training and inference only)
        run_gd_analysis(
            df, gd_table, exclusion_mask, args, device=args.device,
            column_medians=column_medians,
            lowres_median_bin_size=lowres_median_bin_size,
        )

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print("\nOutput files:")
    print("  - ploidy_estimates.tsv")
    print("  - cn_posteriors.tsv.gz")
    print("  - sample_posteriors.tsv.gz")
    print("  - bin_posteriors.tsv.gz")
    print("  - bin_mappings.tsv.gz")
    print("  - locus_intervals.tsv.gz")
    print("\nNext step: Run plot_gd_cnv_output.py to call CNVs and generate plots.")
    print("=" * 80)

    print("\n" + "=" * 80)
    print("Analysis complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
