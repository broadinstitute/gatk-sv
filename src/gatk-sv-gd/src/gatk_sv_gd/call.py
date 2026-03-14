"""
GD CNV calling from model posteriors.

Supports two calling modes:
- ``viterbi``: smooth segmentation using transition matrices
- ``posterior-marginal``: direct scoring from pair-state posterior marginals
"""

import argparse
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from gatk_sv_gd import _util
from gatk_sv_gd.models import GDTable
from gatk_sv_gd.viterbi import (
    load_transition_matrix,
    viterbi_call_gd_cnv,
)


def setup_logging(output_dir: str, filename: str = "call_log.txt"):
    """Redirect stdout/stderr to console and log file."""
    log_path = os.path.join(output_dir, filename)
    log_fh = open(log_path, "w")
    _util._log_fh = log_fh
    sys.stdout = _util.TeeStream(sys.__stdout__, log_fh)
    sys.stderr = _util.TeeStream(sys.__stderr__, log_fh)
    return log_fh


def get_locus_interval_bins(
    bin_mappings_df: pd.DataFrame,
    cluster: str,
) -> Dict[str, List[int]]:
    """Return array_idx values grouped by interval for one locus."""
    locus_bins = bin_mappings_df[bin_mappings_df["cluster"] == cluster]
    interval_bins: Dict[str, List[int]] = {}
    for interval_name, group in locus_bins.groupby("interval"):
        interval_bins[interval_name] = group["array_idx"].tolist()
    return interval_bins


def get_call_confidence(call: dict) -> float:
    """Return the preferred confidence score for a call."""
    confidence = call.get("confidence_score", np.nan)
    if not pd.isna(confidence):
        return float(confidence)
    log_prob_score = call.get("log_prob_score", np.nan)
    if not pd.isna(log_prob_score):
        return float(log_prob_score)
    return float("nan")


def determine_best_breakpoints(
    calls: List[dict],
    calling_mode: str = "viterbi",
    carrier_only: bool = True,
) -> Dict[str, Optional[str]]:
    """Pick the best GD_ID per svtype."""
    best_by_svtype: Dict[str, Optional[str]] = {}
    for svtype in ["DEL", "DUP"]:
        sv_calls = [c for c in calls if c.get("svtype") == svtype]
        if carrier_only:
            sv_calls = [c for c in sv_calls if bool(c.get("is_carrier", False))]

        if not sv_calls:
            best_by_svtype[svtype] = None
            continue

        if calling_mode == "posterior-marginal":
            best = max(
                sv_calls,
                key=lambda c: (
                    get_call_confidence(c),
                    c.get("matched_interval_bp", 0),
                    c.get("end", 0) - c.get("start", 0),
                    str(c.get("GD_ID", "")),
                ),
            )
        else:
            best = max(
                sv_calls,
                key=lambda c: (
                    c.get("matched_interval_bp", 0),
                    c.get("interval_coverage", c.get("reciprocal_overlap", 0.0)),
                    get_call_confidence(c),
                    c.get("end", 0) - c.get("start", 0),
                    str(c.get("GD_ID", "")),
                ),
            )

        best_by_svtype[svtype] = best.get("GD_ID")
    return best_by_svtype


def determine_posterior_carrier_breakpoints(
    calls: List[dict],
    min_interval_confidence: float,
    min_flank_non_event_confidence: float,
) -> Dict[str, Optional[str]]:
    """Pick at most one qualifying posterior-marginal GD_ID per SV type.

    A call qualifies only when each covered interval's mean event probability
    meets the minimum threshold. Among qualifying calls, choose the largest one.
    """
    selected_by_svtype: Dict[str, Optional[str]] = {}
    for svtype in ["DEL", "DUP"]:
        qualifying_calls: List[dict] = []
        for call in calls:
            if call.get("svtype") != svtype:
                continue
            interval_confidences = call.get("interval_confidences", [])
            flank_non_event_medians = [
                call.get("left_flank_non_event_median", np.nan),
                call.get("right_flank_non_event_median", np.nan),
            ]
            flank_pass = all(
                ((not pd.notna(flank_median)) or (
                    float(flank_median) >= min_flank_non_event_confidence
                ))
                for flank_median in flank_non_event_medians
            )
            if interval_confidences and all(
                float(confidence) >= min_interval_confidence
                for confidence in interval_confidences
            ) and flank_pass:
                qualifying_calls.append(call)

        if not qualifying_calls:
            selected_by_svtype[svtype] = None
            continue

        best_call = max(
            qualifying_calls,
            key=lambda call: (
                call.get("matched_interval_bp", 0),
                call.get("end", 0) - call.get("start", 0),
                call.get("n_bins", 0),
                get_call_confidence(call),
                str(call.get("GD_ID", "")),
            ),
        )
        selected_by_svtype[svtype] = str(best_call.get("GD_ID", ""))
    return selected_by_svtype


def get_pair_state_columns(
    cn_posteriors_df: pd.DataFrame,
) -> Tuple[List[str], List[Tuple[int, int]]]:
    """Return pair-state columns and canonicalized labels."""
    pair_cols = [
        column for column in cn_posteriors_df.columns
        if column.startswith("prob_pair_")
    ]
    if not pair_cols:
        raise ValueError(
            "cn_posteriors.tsv.gz is missing pair-state posterior columns "
            "(expected columns like prob_pair_0_1). Re-run infer first."
        )

    canonical_seen: Dict[Tuple[int, int], str] = {}
    canonical_labels: List[Tuple[int, int]] = []
    for column in pair_cols:
        match = re.fullmatch(r"prob_pair_(\d+)_(\d+)", column)
        if match is None:
            raise ValueError(f"Unrecognized pair-state column name: {column}")
        pair = tuple(sorted((int(match.group(1)), int(match.group(2)))))
        if pair in canonical_seen:
            raise ValueError(
                f"Duplicate canonical pair-state labels: {canonical_seen[pair]} and "
                f"{column} both map to {pair}"
            )
        canonical_seen[pair] = column
        canonical_labels.append(pair)

    return pair_cols, canonical_labels


def build_event_pair_mask(
    pair_states: List[Tuple[int, int]],
    svtype: str,
    sample_ploidy: int,
) -> np.ndarray:
    """Return a boolean mask of pair states supporting the event class."""
    if svtype == "DEL":
        return np.array(
            [(h1 + h2) < sample_ploidy for h1, h2 in pair_states],
            dtype=bool,
        )
    if svtype == "DUP":
        return np.array(
            [(h1 + h2) > sample_ploidy for h1, h2 in pair_states],
            dtype=bool,
        )
    raise ValueError(f"Unsupported svtype: {svtype}")


def build_flank_non_event_pair_mask(
    pair_states: List[Tuple[int, int]],
    svtype: str,
    sample_ploidy: int,
) -> np.ndarray:
    """Return a boolean mask of pair states consistent with no event in a flank."""
    if svtype == "DEL":
        return np.array(
            [(h1 + h2) >= sample_ploidy for h1, h2 in pair_states],
            dtype=bool,
        )
    if svtype == "DUP":
        return np.array(
            [(h1 + h2) <= sample_ploidy for h1, h2 in pair_states],
            dtype=bool,
        )
    raise ValueError(f"Unsupported svtype: {svtype}")


def compute_event_marginal_probabilities(
    pair_prob_matrix: np.ndarray,
    pair_states: List[Tuple[int, int]],
    sample_ploidy: int,
) -> Dict[str, np.ndarray]:
    """Return per-bin event marginal probabilities for DEL and DUP."""
    event_probs: Dict[str, np.ndarray] = {}
    for svtype in ("DEL", "DUP"):
        event_mask = build_event_pair_mask(pair_states, svtype, sample_ploidy)
        if pair_prob_matrix.size == 0 or not np.any(event_mask):
            event_probs[svtype] = np.zeros(pair_prob_matrix.shape[0], dtype=float)
            continue
        event_probs[svtype] = np.clip(
            pair_prob_matrix[:, event_mask].sum(axis=1),
            0.0,
            1.0,
        )
    return event_probs


def score_call_from_posterior_marginals(
    locus,
    entry: dict,
    sample_pair_probs: np.ndarray,
    pair_states: List[Tuple[int, int]],
    interval_bin_arrays: Dict[str, np.ndarray],
    sample_ploidy: int,
) -> dict:
    """Score one GD entry directly from pair-state posterior marginals."""
    svtype = str(entry["svtype"])
    bp1 = str(entry["BP1"])
    bp2 = str(entry["BP2"])
    covered_tuples = locus.get_intervals_between(bp1, bp2)
    covered_intervals = [name for _, _, name in covered_tuples]

    covered_bin_indices: List[int] = []
    for _, _, interval_name in covered_tuples:
        covered_bin_indices.extend(interval_bin_arrays.get(interval_name, np.array([], dtype=int)).tolist())

    event_mask = build_event_pair_mask(pair_states, svtype, sample_ploidy)
    flank_non_event_mask = build_flank_non_event_pair_mask(
        pair_states,
        svtype,
        sample_ploidy,
    )
    interval_confidences: List[float] = []
    if np.any(event_mask):
        for interval_name in covered_intervals:
            interval_bin_indices = interval_bin_arrays.get(
                interval_name,
                np.array([], dtype=int),
            )
            if len(interval_bin_indices) == 0:
                interval_confidences.append(0.0)
                continue
            interval_probs = np.clip(
                sample_pair_probs[interval_bin_indices][:, event_mask].sum(axis=1),
                0.0,
                1.0,
            )
            interval_confidences.append(float(np.mean(interval_probs)))

    if covered_bin_indices and np.any(event_mask):
        per_bin_event_probs = np.clip(
            sample_pair_probs[covered_bin_indices][:, event_mask].sum(axis=1),
            0.0,
            1.0,
        )
        confidence = float(np.mean(per_bin_event_probs))
    else:
        confidence = 0.0

    covered_bp_total = int(
        sum(max(0, int(end) - int(start)) for start, end, _ in covered_tuples)
    )

    flank_non_event_medians: Dict[str, float] = {}
    for flank_name in ("left_flank", "right_flank"):
        flank_bin_indices = interval_bin_arrays.get(flank_name, np.array([], dtype=int))
        if len(flank_bin_indices) == 0 or not np.any(flank_non_event_mask):
            flank_non_event_medians[flank_name] = np.nan
            continue
        flank_non_event_probs = np.clip(
            sample_pair_probs[flank_bin_indices][:, flank_non_event_mask].sum(axis=1),
            0.0,
            1.0,
        )
        flank_non_event_medians[flank_name] = float(np.median(flank_non_event_probs))

    valid_flank_medians = [
        value for value in flank_non_event_medians.values()
        if pd.notna(value)
    ]

    return {
        "GD_ID": entry["GD_ID"],
        "chrom": locus.chrom,
        "start": int(entry["start_GRCh38"]),
        "end": int(entry["end_GRCh38"]),
        "svtype": svtype,
        "BP1": bp1,
        "BP2": bp2,
        "is_terminal": locus.is_terminal,
        "n_bins": len(covered_bin_indices),
        "sample_ploidy": sample_ploidy,
        "haplotype": np.nan,
        "hap_cn_state": np.nan,
        "matched_seg_start": np.nan,
        "matched_seg_end": np.nan,
        "matched_seg_n_bins": 0,
        "matched_interval_bp": covered_bp_total,
        "interval_coverage": confidence,
        "reciprocal_overlap": confidence,
        "intervals": covered_intervals,
        "interval_confidences": interval_confidences,
        "min_interval_confidence": (
            float(min(interval_confidences)) if interval_confidences else 0.0
        ),
        "left_flank_non_event_median": flank_non_event_medians.get("left_flank", np.nan),
        "right_flank_non_event_median": flank_non_event_medians.get("right_flank", np.nan),
        "min_flank_non_event_confidence": (
            float(min(valid_flank_medians)) if valid_flank_medians else np.nan
        ),
        "log_prob_score": confidence,
        "confidence_score": confidence,
        "is_carrier": False,
    }


def call_cnvs_from_posteriors(
    cn_posteriors_df: pd.DataFrame,
    bin_mappings_df: pd.DataFrame,
    gd_table: GDTable,
    transition_matrix: Optional[np.ndarray] = None,
    ploidy_df: Optional[pd.DataFrame] = None,
    verbose: bool = False,
    min_mean_coverage: float = 0.90,
    breakpoint_transition_matrix: Optional[np.ndarray] = None,
    calling_mode: str = "viterbi",
    min_posterior_interval_confidence: float = 0.80,
    min_flank_non_event_confidence: float = 0.90,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Call GD CNVs from posterior probabilities."""
    if calling_mode not in {"viterbi", "posterior-marginal"}:
        raise ValueError(
            f"Unsupported calling_mode: {calling_mode}. "
            "Expected 'viterbi' or 'posterior-marginal'."
        )
    if calling_mode == "viterbi" and transition_matrix is None:
        raise ValueError("transition_matrix is required for calling_mode='viterbi'")

    print("\n" + "=" * 80)
    print("CALLING CNVs FROM POSTERIORS")
    if calling_mode == "viterbi":
        bp_str = " + breakpoint matrix" if breakpoint_transition_matrix is not None else ""
        print(
            f"  Calling mode: Viterbi segmentation{bp_str}  "
            f"(minimum per-interval coverage={min_mean_coverage:.0%})"
        )
    else:
        print(
            "  Calling mode: posterior-marginal scoring  "
            "(minimum per-interval confidence="
            f"{min_posterior_interval_confidence:.0%}, "
            "minimum flank non-event confidence="
            f"{min_flank_non_event_confidence:.0%})"
        )
    print("=" * 80)

    all_results: List[dict] = []
    all_path_records: List[dict] = []
    all_event_records: List[dict] = []
    sample_ids = cn_posteriors_df["sample"].unique()

    ploidy_lookup: Dict[Tuple[str, str], int] = {}
    if ploidy_df is not None:
        for _, row in ploidy_df.iterrows():
            ploidy_lookup[(str(row["sample"]), str(row["contig"]))] = int(row["ploidy"])
        print(f"  Loaded ploidy for {len(ploidy_lookup)} sample/contig pairs")
    else:
        print("  No ploidy table provided; assuming diploid (ploidy=2) everywhere")

    n_bins = len(bin_mappings_df)
    n_samples = len(sample_ids)
    expected_rows = n_bins * n_samples
    if len(cn_posteriors_df) != expected_rows:
        print(
            f"  WARNING: cn_posteriors has {len(cn_posteriors_df)} rows, "
            f"expected {expected_rows} ({n_bins} bins × {n_samples} samples)"
        )

    first_sample = sample_ids[0]
    first_sample_rows = cn_posteriors_df[cn_posteriors_df["sample"] == first_sample]
    if len(first_sample_rows) == n_bins:
        post_coords = list(
            zip(
                first_sample_rows["chr"].values,
                first_sample_rows["start"].values,
                first_sample_rows["end"].values,
            )
        )
        map_coords = list(
            zip(
                bin_mappings_df["chr"].values,
                bin_mappings_df["start"].values,
                bin_mappings_df["end"].values,
            )
        )
        if post_coords != map_coords:
            raise ValueError(
                "Bin coordinates in cn_posteriors do not match bin_mappings. "
                "Please re-run infer to regenerate both files."
            )
        print(f"  Validated: bin coordinates match between posteriors and mappings ({n_bins} bins)")
    else:
        print(
            f"  WARNING: sample {first_sample} has {len(first_sample_rows)} bins, "
            f"expected {n_bins}; skipping alignment check"
        )

    print("  Organizing data for fast access...")
    pair_prob_cols, pair_state_labels = get_pair_state_columns(cn_posteriors_df)
    n_pair_states = len(pair_prob_cols)
    pair_prob_3d = np.empty((n_samples, n_bins, n_pair_states))
    depth_2d = np.empty((n_samples, n_bins))
    for s_idx, sample_id in enumerate(sample_ids):
        mask = cn_posteriors_df["sample"] == sample_id
        sample_rows = cn_posteriors_df.loc[mask]
        pair_prob_3d[s_idx] = sample_rows[pair_prob_cols].values
        depth_2d[s_idx] = sample_rows["depth"].values
    print(
        f"    Extracted {n_samples} x {n_bins} x {n_pair_states} "
        "pair-state probability array"
    )

    bin_coords_by_idx: Dict[int, Tuple[int, int]] = dict(
        zip(
            bin_mappings_df["array_idx"].astype(int),
            zip(
                bin_mappings_df["start"].astype(int),
                bin_mappings_df["end"].astype(int),
            ),
        )
    )

    for cluster, locus in gd_table.loci.items():
        print(f"\nCalling CNVs for locus: {cluster}")

        cluster_bin_rows = (
            bin_mappings_df[bin_mappings_df["cluster"] == cluster]
            .drop_duplicates(subset=["array_idx"])
            .sort_values("array_idx")
        )
        cluster_bin_indices = cluster_bin_rows["array_idx"].astype(int).to_numpy()

        interval_bins = get_locus_interval_bins(bin_mappings_df, cluster)
        bp_masked = interval_bins.pop("breakpoint_ranges", [])
        if bp_masked:
            print(f"  Masking {len(bp_masked)} breakpoint-range bin(s)")

        if not interval_bins:
            print("  Skipping — no bins in mappings for this locus")
            continue

        for interval_name, bins in interval_bins.items():
            print(f"  {interval_name}: {len(bins)} bins")

        for flank_name in ("left_flank", "right_flank"):
            if flank_name not in interval_bins or len(interval_bins[flank_name]) == 0:
                warning = f"  WARNING: no {flank_name} bins for {cluster}"
                if calling_mode == "viterbi":
                    warning += " — Viterbi trace will not cover that flank"
                print(warning)

        interval_bin_arrays = {
            name: np.array(bins, dtype=int)
            for name, bins in interval_bins.items()
        }

        all_cluster_idxs = set(
            bin_mappings_df[
                bin_mappings_df["cluster"] == cluster
            ]["array_idx"].astype(int)
        )
        bp_set = set(int(b) for b in bp_masked)
        all_cluster_bins = sorted(all_cluster_idxs - bp_set)

        for s_idx, sample_id in enumerate(sample_ids):
            sample_ploidy = ploidy_lookup.get((str(sample_id), locus.chrom), 2)

            cluster_pair_probs = pair_prob_3d[s_idx, cluster_bin_indices, :]
            cluster_event_probs = compute_event_marginal_probabilities(
                cluster_pair_probs,
                pair_state_labels,
                sample_ploidy,
            )
            for (_, bin_row), del_prob, dup_prob in zip(
                cluster_bin_rows.iterrows(),
                cluster_event_probs["DEL"],
                cluster_event_probs["DUP"],
            ):
                all_event_records.append(
                    {
                        "sample": sample_id,
                        "cluster": cluster,
                        "chrom": locus.chrom,
                        "start": int(bin_row["start"]),
                        "end": int(bin_row["end"]),
                        "prob_del_event": float(del_prob),
                        "prob_dup_event": float(dup_prob),
                    }
                )

            if verbose:
                _util.vlog(
                    f"\n  [{calling_mode}] Sample: {sample_id}  "
                    f"({locus.chrom}, ploidy={sample_ploidy})"
                )

            if calling_mode == "viterbi":
                calls, path_records = viterbi_call_gd_cnv(
                    locus,
                    pair_prob_3d[s_idx],
                    pair_state_labels,
                    transition_matrix,
                    interval_bin_arrays,
                    ploidy=sample_ploidy,
                    min_mean_coverage=min_mean_coverage,
                    verbose=verbose,
                    sample_id=str(sample_id),
                    breakpoint_transition_matrix=breakpoint_transition_matrix,
                    bin_coords=bin_coords_by_idx,
                    all_cluster_bins=all_cluster_bins,
                )
                for start, end, cn_state, category, haplotype in path_records:
                    all_path_records.append(
                        {
                            "sample": sample_id,
                            "cluster": cluster,
                            "start": start,
                            "end": end,
                            "cn_state": cn_state,
                            "category": category,
                            "haplotype": haplotype,
                        }
                    )
                best_by_svtype = determine_best_breakpoints(
                    calls,
                    calling_mode="viterbi",
                    carrier_only=True,
                )
            else:
                calls = [
                    score_call_from_posterior_marginals(
                        locus=locus,
                        entry=entry,
                        sample_pair_probs=pair_prob_3d[s_idx],
                        pair_states=pair_state_labels,
                        interval_bin_arrays=interval_bin_arrays,
                        sample_ploidy=sample_ploidy,
                    )
                    for entry in locus.gd_entries
                ]
                selected_by_svtype = determine_posterior_carrier_breakpoints(
                    calls,
                    min_interval_confidence=min_posterior_interval_confidence,
                    min_flank_non_event_confidence=min_flank_non_event_confidence,
                )
                best_by_svtype = {}
                for call in calls:
                    is_selected = str(call["GD_ID"]) == selected_by_svtype.get(
                        call["svtype"]
                    )
                    call["is_carrier"] = bool(is_selected)
                    call["is_best_match"] = bool(is_selected)

            if verbose:
                for call in calls:
                    tag = " [CARRIER]" if call.get("is_carrier", False) else ""
                    if calling_mode == "viterbi":
                        coverage = call.get(
                            "interval_coverage",
                            call.get("reciprocal_overlap", 0.0),
                        )
                        _util.vlog(
                            f"    [VIT] {sample_id:30s}  "
                            f"{call['GD_ID']:25s}  "
                            f"{call['svtype']:4s}  ploidy={sample_ploidy}  "
                            f"score={get_call_confidence(call):+.4f}  "
                            f"coverage={coverage:.2%}  "
                            f"intervals={','.join(call['intervals'])}"
                            f"{tag}"
                        )
                    else:
                        _util.vlog(
                            f"    [POST] {sample_id:30s}  "
                            f"{call['GD_ID']:25s}  "
                            f"{call['svtype']:4s}  ploidy={sample_ploidy}  "
                            f"confidence={get_call_confidence(call):.4f}  "
                            f"min_interval={call.get('min_interval_confidence', np.nan):.4f}  "
                            f"min_flank_non_event={call.get('min_flank_non_event_confidence', np.nan):.4f}  "
                            f"intervals={','.join(call['intervals'])}"
                            f"{tag}"
                        )

            for call in calls:
                svtype = call["svtype"]
                best_gd_for_svtype = best_by_svtype.get(svtype)

                covered_bin_indices: List[int] = []
                for interval_name in call["intervals"]:
                    if interval_name in interval_bin_arrays:
                        covered_bin_indices.extend(interval_bin_arrays[interval_name].tolist())

                if covered_bin_indices:
                    mean_depth = float(depth_2d[s_idx, covered_bin_indices].mean())
                else:
                    mean_depth = np.nan

                confidence_score = get_call_confidence(call)
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
                    "is_terminal": call["is_terminal"],
                    "n_bins": call["n_bins"],
                    "mean_depth": mean_depth,
                    "sample_ploidy": call.get("sample_ploidy", sample_ploidy),
                    "matched_haplotype": call.get("haplotype", np.nan),
                    "hap_cn_state": call.get("hap_cn_state", np.nan),
                    "matched_seg_start": call.get("matched_seg_start", np.nan),
                    "matched_seg_end": call.get("matched_seg_end", np.nan),
                    "matched_seg_n_bins": call.get("matched_seg_n_bins", 0),
                    "matched_interval_bp": call.get("matched_interval_bp", 0),
                    "interval_coverage": call.get("interval_coverage", np.nan),
                    "reciprocal_overlap": call.get("reciprocal_overlap", np.nan),
                    "min_interval_confidence": call.get("min_interval_confidence", np.nan),
                    "left_flank_non_event_median": call.get(
                        "left_flank_non_event_median",
                        np.nan,
                    ),
                    "right_flank_non_event_median": call.get(
                        "right_flank_non_event_median",
                        np.nan,
                    ),
                    "min_flank_non_event_confidence": call.get(
                        "min_flank_non_event_confidence",
                        np.nan,
                    ),
                    "is_carrier": bool(call.get("is_carrier", False)),
                    "is_best_match": (
                        bool(call.get("is_best_match", False))
                        if calling_mode == "posterior-marginal"
                        else (
                            call["GD_ID"] == best_gd_for_svtype
                            if best_gd_for_svtype is not None else False
                        )
                    ),
                    "log_prob_score": call.get("log_prob_score", confidence_score),
                    "confidence_score": confidence_score,
                    "calling_method": calling_mode,
                }
                all_results.append(result)

    calls_df = pd.DataFrame(
        all_results,
        columns=[
            "sample",
            "cluster",
            "GD_ID",
            "chrom",
            "start",
            "end",
            "svtype",
            "BP1",
            "BP2",
            "is_terminal",
            "n_bins",
            "mean_depth",
            "sample_ploidy",
            "matched_haplotype",
            "hap_cn_state",
            "matched_seg_start",
            "matched_seg_end",
            "matched_seg_n_bins",
            "matched_interval_bp",
            "interval_coverage",
            "reciprocal_overlap",
            "min_interval_confidence",
            "left_flank_non_event_median",
            "right_flank_non_event_median",
            "min_flank_non_event_confidence",
            "is_carrier",
            "is_best_match",
            "log_prob_score",
            "confidence_score",
            "calling_method",
        ],
    )
    paths_df = pd.DataFrame(
        all_path_records,
        columns=[
            "sample",
            "cluster",
            "start",
            "end",
            "cn_state",
            "category",
            "haplotype",
        ],
    )
    event_marginals_df = pd.DataFrame(
        all_event_records,
        columns=[
            "sample",
            "cluster",
            "chrom",
            "start",
            "end",
            "prob_del_event",
            "prob_dup_event",
        ],
    )
    return calls_df, paths_df, event_marginals_df


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Call GD CNVs from model posterior probabilities",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--cn-posteriors",
        required=True,
        help="CN posteriors file (cn_posteriors.tsv.gz) with depth values",
    )
    parser.add_argument(
        "--bin-mappings",
        required=True,
        help="Bin mappings file (bin_mappings.tsv.gz) from gd_cnv_pyro.py",
    )
    parser.add_argument(
        "--gd-table", "-g",
        required=True,
        help="GD locus definition table (TSV)",
    )
    parser.add_argument(
        "--ploidy-table",
        required=False,
        help="Ploidy estimates TSV (ploidy_estimates.tsv) from gd_cnv_pyro.py. "
             "Columns: sample, contig, median_depth, ploidy. "
             "If not provided, ploidy=2 is assumed for all sample/contig pairs.",
    )
    parser.add_argument(
        "--output-dir", "-o",
        required=True,
        help="Output directory for calls",
    )
    parser.add_argument(
        "--transition-matrix",
        required=False,
        help="CN-state transition probability matrix (TSV) for Viterbi calling only. "
             "Required when --calling-mode=viterbi.",
    )
    parser.add_argument(
        "--breakpoint-transition-matrix",
        required=False,
        help="CN-state transition probability matrix (TSV) applied at known "
             "breakpoint boundaries during Viterbi calling only.",
    )
    parser.add_argument(
        "--min-mean-coverage",
        type=float,
        default=0.50,
        help="Minimum per-interval coverage required for a breakpoint pair "
             "to participate in the selected contiguous run in Viterbi mode.",
    )
    parser.add_argument(
        "--calling-mode",
        choices=["viterbi", "posterior-marginal"],
        default="posterior-marginal",
        help="Calling strategy to use.",
    )
    parser.add_argument(
        "--min-posterior-interval-confidence",
        type=float,
        default=0.80,
        help="In posterior-marginal mode, mark every breakpoint combination as a "
             "carrier when each covered interval's mean event probability meets "
             "or exceeds this threshold.",
    )
    parser.add_argument(
        "--min-flank-non-event-confidence",
        type=float,
        default=0.90,
        help="In posterior-marginal mode, reject a breakpoint combination when any "
             "available flank has median non-event probability below this threshold.",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print detailed per-sample scores for all GD entries at every locus.",
    )
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    setup_logging(args.output_dir)

    print(f"Output directory: {args.output_dir}")
    print(f"Calling mode: {args.calling_mode}")

    print("\nLoading data...")

    print(f"  Loading CN posteriors: {args.cn_posteriors}")
    cn_posteriors_df = pd.read_csv(args.cn_posteriors, sep="\t", compression="infer")
    print(f"    {len(cn_posteriors_df)} bin-sample records")

    print(f"  Loading bin mappings: {args.bin_mappings}")
    bin_mappings_df = pd.read_csv(args.bin_mappings, sep="\t", compression="infer")
    print(f"    {len(bin_mappings_df)} bin mappings")

    print(f"  Loading GD table: {args.gd_table}")
    gd_table = GDTable(args.gd_table)
    print(f"    {len(gd_table.loci)} loci")

    ploidy_df = None
    if args.ploidy_table:
        print(f"  Loading ploidy table: {args.ploidy_table}")
        ploidy_df = pd.read_csv(args.ploidy_table, sep="\t")
        print(f"    {len(ploidy_df)} sample/contig ploidy records")

    transition_matrix = None
    breakpoint_transition_matrix = None
    if args.calling_mode == "viterbi":
        if not args.transition_matrix:
            raise ValueError("--transition-matrix is required when --calling-mode=viterbi")
        print(f"\n  Loading transition matrix: {args.transition_matrix}")
        transition_matrix = load_transition_matrix(args.transition_matrix)

        if args.breakpoint_transition_matrix:
            print(
                "  Loading breakpoint transition matrix: "
                f"{args.breakpoint_transition_matrix}"
            )
            breakpoint_transition_matrix = load_transition_matrix(
                args.breakpoint_transition_matrix
            )

    calls_df, paths_df, event_marginals_df = call_cnvs_from_posteriors(
        cn_posteriors_df,
        bin_mappings_df,
        gd_table,
        transition_matrix=transition_matrix,
        ploidy_df=ploidy_df,
        verbose=args.verbose,
        min_mean_coverage=args.min_mean_coverage,
        breakpoint_transition_matrix=breakpoint_transition_matrix,
        calling_mode=args.calling_mode,
        min_posterior_interval_confidence=args.min_posterior_interval_confidence,
        min_flank_non_event_confidence=args.min_flank_non_event_confidence,
    )

    output_file = os.path.join(args.output_dir, "gd_cnv_calls.tsv.gz")
    calls_df.to_csv(output_file, sep="\t", index=False, compression="gzip")
    print(f"\n  Saved calls to: {output_file}")
    print(f"    {len(calls_df)} call records")

    paths_file = os.path.join(args.output_dir, "viterbi_paths.tsv.gz")
    paths_df.to_csv(paths_file, sep="\t", index=False, compression="gzip")
    print(f"  Saved Viterbi paths to: {paths_file}")
    print(f"    {len(paths_df)} path records")

    event_marginals_file = os.path.join(args.output_dir, "event_marginals.tsv.gz")
    event_marginals_df.to_csv(
        event_marginals_file,
        sep="\t",
        index=False,
        compression="gzip",
    )
    print(f"  Saved event marginals to: {event_marginals_file}")
    print(f"    {len(event_marginals_df)} bin-sample records")

    if len(calls_df) > 0:
        carriers = calls_df[calls_df["is_carrier"]]
        n_carriers = carriers["sample"].nunique()
        n_sites = carriers["GD_ID"].nunique()
        print(f"\n  {n_carriers} carrier samples across {n_sites} GD sites")
    else:
        print(
            "\n  No calls produced — check that the GD table, bin mappings, "
            "and CN posteriors refer to the same loci."
        )

    print("\n" + "=" * 80)
    print("Calling complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
