"""
Output writers for GD CNV inference results.

Functions for writing posterior tables, locus metadata, and ploidy estimates
to disk after model inference.
"""

import os
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from gatk_sv_gd._util import get_sample_columns
from gatk_sv_gd.bins import LocusBinMapping
from gatk_sv_gd.models import GDLocus


def write_posterior_tables(
    combined_data,
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

    # 3. Write GD-entry → interval mapping
    #
    # This is the crucial bridge between the GD table (identified by GD_ID)
    # and the posteriors / bin_mappings (identified by cluster + interval).
    #
    # A single GD entry may span *multiple* sub-intervals in the posteriors.
    # For example, a BP1→BP3 entry in a cluster that also defines BP2 will
    # cover both the "1-2" and "2-3" intervals, neither of which is named
    # "1-3".  Without this file users have no way to map a GD_ID to the
    # posterior rows that carry its copy-number signal.
    #
    # Columns:
    #   GD_ID        – identifier from the input GD table
    #   cluster      – cluster key (matches cn_posteriors.cluster)
    #   svtype       – DEL / DUP
    #   BP1, BP2     – breakpoint names for this GD entry
    #   interval     – sub-interval name (matches cn_posteriors.interval)
    #   chr          – chromosome
    #   interval_start, interval_end – genomic coordinates of the sub-interval
    gd_entry_rows = []
    for cluster, locus in included_loci.items():
        for entry in locus.gd_entries:
            covered = locus.get_intervals_between(entry["BP1"], entry["BP2"])
            if not covered:
                # BP1/BP2 names not found in this locus (should not happen for
                # well-formed input, but guard against it gracefully).
                print(f"  WARNING: GD entry {entry['GD_ID']} has BP1={entry['BP1']}, "
                      f"BP2={entry['BP2']} which do not resolve to any interval in "
                      f"cluster {cluster} — omitting from gd_entry_intervals.tsv.gz")
                continue
            for iv_start, iv_end, iv_name in covered:
                gd_entry_rows.append({
                    "GD_ID": entry["GD_ID"],
                    "cluster": cluster,
                    "svtype": entry["svtype"],
                    "BP1": entry["BP1"],
                    "BP2": entry["BP2"],
                    "interval": iv_name,
                    "chr": locus.chrom,
                    "interval_start": iv_start,
                    "interval_end": iv_end,
                })

    gd_entry_df = pd.DataFrame(gd_entry_rows)
    gd_entry_output = os.path.join(output_dir, "gd_entry_intervals.tsv.gz")
    gd_entry_df.to_csv(gd_entry_output, sep="\t", index=False, compression="gzip")
    print(f"  Saved: {gd_entry_output}")
    print(f"  Rows: {len(gd_entry_df):,} (GD entry × interval)")


def estimate_ploidy(
    df: pd.DataFrame,
    output_dir: str,
) -> pd.DataFrame:
    """
    Estimate ploidy for each sample/contig pair from the filtered bin set.

    Uses the median normalized depth across all filtered bins on each
    chromosome for each sample. Since the data is normalized so that CN=2
    corresponds to a depth of 2.0, the rounded median gives the ploidy.

    Args:
        df: Filtered DataFrame with bins as rows and samples as columns.
            Expected to already be normalized so diploid depth ≈ 2.0.
        output_dir: Directory to write the ploidy table.

    Returns:
        DataFrame with columns: sample, contig, median_depth, ploidy
    """
    sample_cols = get_sample_columns(df)

    print(f"\n{'=' * 80}")
    print("ESTIMATING PLOIDY PER SAMPLE / CONTIG")
    print(f"{'=' * 80}")

    rows = []
    for contig, contig_df in df.groupby("Chr"):
        depths = contig_df[sample_cols].values  # bins × samples
        medians = np.median(depths, axis=0)     # per-sample median
        for sample_id, med in zip(sample_cols, medians):
            ploidy = int(np.round(med))
            rows.append({
                "sample": sample_id,
                "contig": contig,
                "median_depth": float(med),
                "ploidy": ploidy,
            })

    ploidy_df = pd.DataFrame(rows)

    # Summary
    n_samples = len(sample_cols)
    n_contigs = ploidy_df["contig"].nunique()
    print(f"  Samples: {n_samples}")
    print(f"  Contigs: {n_contigs}")
    for contig in sorted(ploidy_df["contig"].unique(),
                         key=lambda c: (0, int(c.replace('chr', ''))) if c.replace('chr', '').isdigit() else (1, c)):
        sub = ploidy_df[ploidy_df["contig"] == contig]
        counts = sub["ploidy"].value_counts().sort_index()
        dist_str = ", ".join(f"ploidy {p}: {n}" for p, n in counts.items())
        print(f"    {contig}: {dist_str}")

    # Write to disk
    output_path = os.path.join(output_dir, "ploidy_estimates.tsv")
    ploidy_df.to_csv(output_path, sep="\t", index=False)
    print(f"  Saved: {output_path}")
    print(f"  Rows: {len(ploidy_df):,}")
    print(f"{'=' * 80}\n")

    return ploidy_df


def build_ploidy_map(ploidy_df: pd.DataFrame) -> Dict[Tuple[str, str], int]:
    """Build a ``{(sample, chrom): ploidy}`` lookup from the ploidy table.

    This is the format consumed by the quality-filtering helpers in
    :mod:`gatk_sv_gd.bins` for ploidy-adjusted median/MAD computation.
    """
    return {
        (str(row["sample"]), str(row["contig"])): int(row["ploidy"])
        for _, row in ploidy_df.iterrows()
    }
