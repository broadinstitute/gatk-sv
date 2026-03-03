#!/bin/python

"""
GD CNV Calling from Model Posteriors

Call genomic disorder (GD) copy-number variants from posterior probabilities
produced by gd_cnv_pyro.py.

Two calling strategies are available:

1. **Mean log-probability** (original) – for each GD entry, averages the
   posterior probability across bins in the affected interval and thresholds
   on a log-probability score.

2. **Viterbi segmentation** – runs the Viterbi algorithm on per-bin CN
   posteriors using a user-supplied state transition matrix, then checks
   whether the resulting copy-number segmentation exactly matches the
   expected breakpoint pattern for each GD entry (including flanks).
   A ``--viterbi-confidence-threshold`` filters on the Viterbi path
   confidence.

When ``--transition-matrix`` is provided the Viterbi strategy is used;
otherwise the mean log-probability strategy is used.

Usage:
    python gd_cnv_call.py \
        --cn-posteriors cn_posteriors.tsv.gz \
        --bin-mappings bin_mappings.tsv.gz \
        --gd-table gd_table.tsv \
        --transition-matrix transition_probs.tsv \
        --output-dir results/
"""

import argparse
import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

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

    def __getattr__(self, name):
        return getattr(self.original_stream, name)


def setup_logging(output_dir: str, filename: str = "call_log.txt"):
    """Redirect stdout and stderr to both the console and a log file.
    Also sets the module-level _log_fh used by _vlog() for log-only output.
    """
    global _log_fh
    log_path = os.path.join(output_dir, filename)
    log_fh = open(log_path, "w")
    _log_fh = log_fh
    sys.stdout = TeeStream(sys.__stdout__, log_fh)
    sys.stderr = TeeStream(sys.__stderr__, log_fh)
    return log_fh


# Module-level handle for log-only (verbose) output.
# Set by setup_logging(); written to directly by _vlog() so that
# verbose diagnostics appear in the log file but NOT on stdout.
_log_fh = None


def _vlog(msg: str) -> None:
    """Write *msg* to the log file only (never to stdout).

    Used for verbose diagnostic output so that it doesn't flood the
    terminal while still being fully captured in the run log.
    """
    if _log_fh is not None:
        _log_fh.write(msg + "\n")
        _log_fh.flush()


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

    sample_posteriors = (
        cn_posteriors_df[cn_posteriors_df["sample"] == sample_id]
        .reset_index(drop=True)
    )

    prob_cols = [c for c in sample_posteriors.columns if c.startswith("prob_cn_")]
    n_states = len(prob_cols)

    for interval_name, bin_indices in interval_bins.items():
        interval_data = sample_posteriors.iloc[bin_indices]

        if len(interval_data) == 0:
            result[interval_name] = {
                "n_bins": 0,
                "cn_probs": np.zeros(n_states),
            }
            continue

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
    verbose: bool = False,
    sample_id: str = "",
) -> List[dict]:
    """
    Call GD CNVs based on interval copy number statistics.

    For each GD entry in the locus (each DEL/DUP definition), check if the
    copy number in the corresponding region supports a CNV call.

    Args:
        locus: GDLocus object
        interval_stats: Dict mapping interval name to CN statistics
        log_prob_threshold: Minimum log probability score to call a CNV
        flanking_log_prob_threshold: Minimum log probability score in flanking
            regions to classify as spanning
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

        covered_tuples = locus.get_intervals_between(bp1, bp2)
        covered_intervals = [name for _, _, name in covered_tuples]

        if len(covered_intervals) == 0:
            print(f"  WARNING [MLP]: {gd_id} ({svtype}) BP1={bp1}, BP2={bp2} "
                  f"does not resolve to any interval in cluster '{locus.cluster}' "
                  f"(breakpoints: {locus.breakpoint_names}) — entry skipped.")
            continue

        if verbose:
            _vlog(f"  [MLP] {sample_id}  {locus.cluster}  {gd_id} "
                  f"({svtype}  BP1={bp1}, BP2={bp2}  ploidy={ploidy}):")
            _vlog(f"    Covered intervals: {sorted(covered_intervals)}")

        # Compute weighted average log probability score
        log_prob_score = 0.0
        total_weight = 0.0

        if svtype == "DEL":
            for interval in covered_intervals:
                if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                    probs = interval_stats[interval]["cn_probs"]
                    weight = interval_stats[interval]["n_bins"]
                    raw_p = probs[:ploidy].sum()
                    p_del = max(raw_p, 1e-5)
                    if p_del > 0:
                        log_prob_score += weight * np.log(p_del)
                        total_weight += weight
                    if verbose:
                        _vlog(f"      {interval}: {weight} bins  "
                              f"p_del={raw_p:.4f}  "
                              f"log(p_del)={np.log(p_del):+.4f}  "
                              f"weighted={weight * np.log(p_del):+.4f}")
                elif verbose:
                    n = interval_stats[interval]["n_bins"] if interval in interval_stats else 0
                    _vlog(f"      {interval}: {n} bins -> skipped (no data)")

        elif svtype == "DUP":
            for interval in covered_intervals:
                if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                    probs = interval_stats[interval]["cn_probs"]
                    weight = interval_stats[interval]["n_bins"]
                    raw_p = probs[ploidy + 1:].sum()
                    p_dup = max(raw_p, 1e-5)
                    if p_dup > 0:
                        log_prob_score += weight * np.log(p_dup)
                        total_weight += weight
                    if verbose:
                        _vlog(f"      {interval}: {weight} bins  "
                              f"p_dup={raw_p:.4f}  "
                              f"log(p_dup)={np.log(p_dup):+.4f}  "
                              f"weighted={weight * np.log(p_dup):+.4f}")
                elif verbose:
                    n = interval_stats[interval]["n_bins"] if interval in interval_stats else 0
                    _vlog(f"      {interval}: {n} bins -> skipped (no data)")

        if total_weight > 0:
            log_prob_score = log_prob_score / total_weight
        else:
            log_prob_score = np.nan

        if verbose:
            score_str = f"{log_prob_score:+.4f}" if not np.isnan(log_prob_score) else "NaN"
            if np.isnan(log_prob_score):
                verdict = "FAIL (no data)"
            elif log_prob_score > log_prob_threshold:
                verdict = f"PASS (> {log_prob_threshold:+.4f})"
            else:
                verdict = f"FAIL (<= {log_prob_threshold:+.4f})"
            _vlog(f"    Score: {score_str} -> {verdict}")

        is_cnv = log_prob_score > log_prob_threshold

        # Check flanking regions for evidence of a large spanning CNV
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
                if flanking_log_prob_score > flanking_log_prob_threshold:
                    is_spanning = True
                    is_cnv = False
                    if verbose:
                        _vlog(f"    Spanning: flank_score={flanking_log_prob_score:+.4f} > "
                              f"threshold={flanking_log_prob_threshold:+.4f} "
                              f"-> reclassified as SPANNING (not CARRIER)")

        if verbose:
            result = "CARRIER" if is_cnv else ("SPANNING" if is_spanning else "NOT_CARRIER")
            _vlog(f"    -> {result}")

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
    Determine the most likely breakpoint pair for a GD CNV, separately for
    DEL and DUP.

    Args:
        locus: GDLocus object
        interval_stats: Dict mapping interval name to CN statistics
        calls: List of CNV call dictionaries
        ploidy: Expected copy number for this sample on this chromosome

    Returns:
        Dict mapping svtype ("DEL", "DUP") to best matching GD_ID or None
    """
    best_by_svtype = {}

    all_intervals = set(name for _, _, name in locus.get_intervals())

    for svtype in ["DEL", "DUP"]:
        carrier_calls = [c for c in calls if c["is_carrier"] and c["svtype"] == svtype]

        if len(carrier_calls) == 0:
            best_by_svtype[svtype] = None
            continue

        if len(carrier_calls) == 1:
            best_by_svtype[svtype] = carrier_calls[0]["GD_ID"]
            continue

        best_score = -np.inf
        best_size = -1
        best_gd_id = None

        for call in carrier_calls:
            covered_intervals = set(call["intervals"])
            uncovered_intervals = all_intervals - covered_intervals

            score = 0.0
            total_weight = 0.0

            if svtype == "DEL":
                for interval in covered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = 1.
                        p_del = probs[:ploidy].sum()
                        if p_del > 0:
                            score += weight * np.log(p_del)
                            total_weight += weight

                for interval in uncovered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = 1.
                        p_normal = probs[ploidy:].sum()
                        if p_normal > 0:
                            score += weight * np.log(p_normal)
                            total_weight += weight

            else:  # DUP
                for interval in covered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = 1.
                        p_dup = probs[ploidy + 1:].sum()
                        if p_dup > 0:
                            score += weight * np.log(p_dup)
                            total_weight += weight

                for interval in uncovered_intervals:
                    if interval in interval_stats and interval_stats[interval]["n_bins"] > 0:
                        probs = interval_stats[interval]["cn_probs"]
                        weight = 1.
                        p_normal = probs[:ploidy + 1].sum()
                        if p_normal > 0:
                            score += weight * np.log(p_normal)
                            total_weight += weight

            if total_weight > 0:
                score = score / total_weight

            call_size = call["end"] - call["start"]
            if score > best_score or (score == best_score and call_size > best_size):
                best_score = score
                best_size = call_size
                best_gd_id = call["GD_ID"]

        best_by_svtype[svtype] = best_gd_id

    return best_by_svtype


# =============================================================================
# Viterbi HMM Segmentation
# =============================================================================


def load_transition_matrix(filepath: str) -> np.ndarray:
    """Load a CN-state transition probability matrix from a TSV file.

    The file should be a square matrix of transition probabilities where
    rows are the *from* state and columns are the *to* state.  The matrix
    must be square and each row must sum to 1 (within tolerance).

    Accepted formats:
      - Plain NxN numeric matrix (no header, no index).
      - NxN matrix with a header row and/or an index column whose values
        match ``CN0, CN1, ...`` or ``0, 1, ...``.

    Returns:
        (n_states, n_states) numpy array of transition probabilities.
    """
    # Try to read with a header first; if the first row is all numeric,
    # fall back to header=None so we don't lose a data row.
    df_with_header = pd.read_csv(filepath, sep="\t")
    first_row_numeric = pd.to_numeric(
        df_with_header.iloc[0], errors="coerce"
    ).notna().all()

    if first_row_numeric:
        # No header row — re-read without one
        df = pd.read_csv(filepath, sep="\t", header=None)
    else:
        df = df_with_header

    # Drop a leading index column if present:
    #   - string column whose name/values look like a state label
    #   - first column is object dtype (state names) while rest are numeric
    first_col = df.columns[0]
    is_string_index = (
        isinstance(first_col, str)
        and (first_col == "" or first_col.lower() in ("state", "cn", "index"))
    ) or df.iloc[:, 0].dtype == object
    if is_string_index:
        df = df.iloc[:, 1:]

    mat = df.values.astype(float)
    if mat.shape[0] != mat.shape[1]:
        raise ValueError(
            f"Transition matrix must be square, got shape {mat.shape}"
        )

    # Ensure no zeros (add tiny floor for numerical stability), then normalize
    mat = np.maximum(mat, 1e-30)
    row_sums = mat.sum(axis=1, keepdims=True)
    mat = mat / row_sums

    print(f"  Loaded {mat.shape[0]}x{mat.shape[1]} transition matrix from {filepath}")
    return mat


def run_viterbi(
    state_log_priors: np.ndarray,
    transition_matrix: np.ndarray,
    initial_log_probs: Optional[np.ndarray] = None,
    breakpoint_mask: Optional[np.ndarray] = None,
    breakpoint_transition_matrix: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, float]:
    """Run the Viterbi algorithm using per-bin CN posterior probabilities as
    state priors.

    The Pyro model posterior p(CN=k | data) at each bin is used directly as
    the prior probability of the hidden CN state at that bin.  The transition
    matrix controls the stickiness of the segmentation: high diagonal values
    keep the state constant across adjacent bins; off-diagonal values allow
    state changes.

    The joint probability being maximised is:

        P(path) ∝ ∏_t prior_t(s_t)  ×  ∏_t A(s_{t-1}, s_t)

    where prior_t(s_t) = p(CN=s_t | data, bin t)  (the Pyro posterior).

    When *breakpoint_mask* and *breakpoint_transition_matrix* are both
    supplied, the breakpoint matrix is used in place of *transition_matrix*
    for every step ``t`` where ``breakpoint_mask[t-1]`` is True.  This
    allows a less sticky (more permissive) transition at known recurrent
    breakpoint boundaries so that CN state changes are actively encouraged
    at those positions.

    Args:
        state_log_priors: (T, K) array of log posterior probabilities from
            the Pyro model — i.e. log p(CN=k | data) for each of T bins and
            K copy-number states.  These act as per-bin priors on the hidden
            state, not as emission likelihoods.
        transition_matrix: (K, K) transition probability matrix where entry
            (i, j) is the probability of moving from state i to state j
            between adjacent bins.  High diagonal entries = sticky/smooth
            segmentation.
        initial_log_probs: (K,) log probabilities for the initial hidden
            state.  If None, a uniform prior is used.
        breakpoint_mask: Boolean array of shape (T-1,).  Element ``t-1`` is
            True when the transition from bin ``t-1`` to bin ``t`` crosses a
            known recurrent breakpoint boundary.  Used together with
            *breakpoint_transition_matrix*.
        breakpoint_transition_matrix: (K, K) transition matrix used at steps
            flagged by *breakpoint_mask*.  Should have lower diagonal values
            than *transition_matrix* to reflect a higher prior probability of
            a CN-state change at known breakpoint positions.

    Returns:
        Tuple of:
          - (T,) integer array of the most-likely CN state at each bin.
          - float: mean per-bin log-prior along the Viterbi path, used as
            a confidence score.
    """
    T, K = state_log_priors.shape
    log_trans = np.log(transition_matrix)  # (K, K)

    # Precompute breakpoint-specific log-transition matrix if requested.
    log_bp_trans: Optional[np.ndarray] = None
    if breakpoint_mask is not None and breakpoint_transition_matrix is not None:
        log_bp_trans = np.log(breakpoint_transition_matrix)

    if initial_log_probs is None:
        initial_log_probs = np.full(K, -np.log(K))

    # Viterbi tables
    V = np.empty((T, K))               # best log-prob ending in state k at time t
    ptr = np.empty((T, K), dtype=int)  # back-pointer

    # t=0: joint = initial_prior * state_prior_at_bin_0
    V[0] = initial_log_probs + state_log_priors[0]

    for t in range(1, T):
        # For each current state k, find the best previous state j:
        #   V[t-1, j] + log A(j→k)  for all j
        # Shape: (K_prev, 1) + (K_prev, K_cur) -> (K_prev, K_cur)
        # Use the breakpoint matrix when crossing a known BP boundary.
        lt = (
            log_bp_trans
            if (log_bp_trans is not None and breakpoint_mask[t - 1])
            else log_trans
        )
        scores = V[t - 1, :, np.newaxis] + lt
        ptr[t] = np.argmax(scores, axis=0)
        V[t] = scores[ptr[t], np.arange(K)] + state_log_priors[t]

    # Back-trace to recover the most-likely path
    path = np.empty(T, dtype=int)
    path[-1] = np.argmax(V[-1])
    for t in range(T - 2, -1, -1):
        path[t] = ptr[t + 1, path[t + 1]]

    # Confidence: mean per-bin log-prior along the chosen path.
    # Higher (less negative) = the Pyro posteriors strongly supported
    # the called states at every bin.
    path_log_priors = state_log_priors[np.arange(T), path]
    confidence = float(np.mean(path_log_priors))

    return path, confidence


def _extract_segments(
    path: np.ndarray,
) -> List[Tuple[int, int, int]]:
    """Extract contiguous segments from a Viterbi state path.

    Returns:
        List of (start_bin_idx, end_bin_idx_exclusive, state) tuples.
    """
    segments = []
    if len(path) == 0:
        return segments
    seg_start = 0
    seg_state = int(path[0])
    for i in range(1, len(path)):
        if int(path[i]) != seg_state:
            segments.append((seg_start, i, seg_state))
            seg_start = i
            seg_state = int(path[i])
    segments.append((seg_start, len(path), seg_state))
    return segments


def viterbi_call_gd_cnv(
    locus: GDLocus,
    sample_probs: np.ndarray,
    transition_matrix: np.ndarray,
    interval_bin_arrays: Dict[str, np.ndarray],
    ploidy: int = 2,
    confidence_threshold: float = -1.0,
    flank_coverage_threshold: float = 0.70,
    verbose: bool = False,
    sample_id: str = "",
    breakpoint_transition_matrix: Optional[np.ndarray] = None,
    bin_coords: Optional[Dict[int, Tuple[int, int]]] = None,
) -> List[dict]:
    """Call GD CNVs using Viterbi segmentation of CN posteriors.

    The Pyro model posteriors p(CN=k | data) are used as per-bin priors on
    the hidden CN state.  The transition matrix controls segmentation
    stickiness — high diagonal values prevent single-bin state flips and
    produce smooth copy-number segments.

    The Viterbi path represents the most probable smooth copy-number
    segmentation consistent with both the per-bin posteriors (from Pyro) and
    the expected smoothness encoded in the transition matrix.

    For each GD entry the segmentation is checked against the expected
    breakpoint pattern:

    - **Covered intervals** (between BP1 and BP2):
      - DEL: every bin's Viterbi state must be < ploidy.
      - DUP: every bin's Viterbi state must be > ploidy.
      Multi-step duplications (e.g. CN 3→4) are allowed.

    - **Flanking / uncovered intervals** are checked with a fractional
      threshold: at least ``flank_coverage_threshold`` fraction of bins
      must satisfy the reference-CN requirement.
      - DEL: fraction of bins with state ≥ ploidy must be ≥ threshold.
      - DUP: fraction of bins with state ≤ ploidy must be ≥ threshold.

    Args:
        locus: GDLocus object.
        sample_probs: (n_bins, n_states) posterior probabilities for one
            sample across all bins at this locus (in array_idx order).
            These are p(CN=k | data) from the Pyro model.
        transition_matrix: (n_states, n_states) transition matrix.  High
            diagonal = sticky segmentation.
        interval_bin_arrays: Dict mapping interval name to array of bin
            array_idx values.
        ploidy: Expected reference copy number for this sample/chromosome.
        confidence_threshold: Minimum mean per-bin log-prior along the
            Viterbi path to make a carrier call.  More negative = more
            permissive.
        flank_coverage_threshold: Fraction of bins in each flanking /
            uncovered interval that must be at the reference CN for the
            flank check to pass.  1.0 = all bins must pass (strict);
            0.0 = no bins need to pass (disabled).  Default 0.70.
        breakpoint_transition_matrix: Optional (n_states, n_states) transition
            matrix applied at known recurrent breakpoint boundaries (both the
            start and end coordinate of each breakpoint's SD-block range).
            Should have lower diagonal values than *transition_matrix* so that
            CN-state changes are more favoured at these positions.  Requires
            *bin_coords* to locate boundaries.
        bin_coords: Mapping from array_idx (int) to (genomic_start, genomic_end)
            for each bin.  Used to determine which bin-to-bin transitions cross
            a known breakpoint boundary.  Required when
            *breakpoint_transition_matrix* is set.

    Returns:
        List of call dicts (one per GD entry).
    """
    # Collect all bin indices for this locus, sorted by array_idx
    all_bins = []
    for bin_idx in interval_bin_arrays.values():
        all_bins.extend(int(b) for b in bin_idx)
    all_bins = sorted(set(all_bins))

    if len(all_bins) == 0:
        if verbose:
            _vlog(f"  [VIT] {sample_id}  {locus.cluster}: no bins -> skip")
        return []

    n_states = sample_probs.shape[1]

    # Build per-bin state log-priors from Pyro posteriors
    priors = sample_probs[all_bins]          # (n_locus_bins, n_states)
    priors = np.maximum(priors, 1e-30)       # numerical floor
    state_log_priors = np.log(priors)        # log p(CN=k | data)

    # Subset transition matrix to the number of CN states in this run
    tm = transition_matrix[:n_states, :n_states]

    # Build a per-step breakpoint mask: True where a transition crosses the
    # START or END coordinate of a known breakpoint SD-block range.  At those
    # positions the breakpoint_transition_matrix (less sticky) will be used
    # instead of the regular transition_matrix, encoding a higher prior on
    # CN-state changes at known recurrent breakpoints.
    breakpoint_mask: Optional[np.ndarray] = None
    if (
        breakpoint_transition_matrix is not None
        and bin_coords is not None
        and len(all_bins) > 1
    ):
        # Collect both the start and end coordinate of every breakpoint range.
        bp_boundaries: List[int] = []
        for bp_start, bp_end in locus.breakpoints:
            bp_boundaries.append(bp_start)
            bp_boundaries.append(bp_end)
        bp_boundaries_sorted = sorted(set(bp_boundaries))

        breakpoint_mask = np.zeros(len(all_bins) - 1, dtype=bool)
        for t in range(1, len(all_bins)):
            coords_prev = bin_coords.get(all_bins[t - 1])
            coords_curr = bin_coords.get(all_bins[t])
            if coords_prev is None or coords_curr is None:
                continue
            prev_mid = (coords_prev[0] + coords_prev[1]) / 2.0
            curr_mid = (coords_curr[0] + coords_curr[1]) / 2.0
            lo, hi = min(prev_mid, curr_mid), max(prev_mid, curr_mid)
            for boundary in bp_boundaries_sorted:
                if lo < boundary <= hi:
                    breakpoint_mask[t - 1] = True
                    break

        if verbose:
            n_bp = int(breakpoint_mask.sum())
            if n_bp > 0:
                bp_steps = [
                    (all_bins[t], all_bins[t + 1])
                    for t in range(len(all_bins) - 1)
                    if breakpoint_mask[t]
                ]
                _vlog(
                    f"  [VIT] {sample_id}  {locus.cluster}: "
                    f"{n_bp} transition(s) will use breakpoint matrix "
                    f"(array_idx pairs: "
                    f"{bp_steps[:5]}{'...' if len(bp_steps) > 5 else ''})"
                )

    # Initial hidden-state prior: strong prior on reference CN (ploidy)
    initial_log_probs = np.full(n_states, np.log(1e-6))
    if ploidy < n_states:
        initial_log_probs[ploidy] = np.log(1.0 - 1e-6 * (n_states - 1))

    # Subset breakpoint matrix to n_states if provided
    bp_tm = (
        breakpoint_transition_matrix[:n_states, :n_states]
        if breakpoint_transition_matrix is not None else None
    )

    path, confidence = run_viterbi(
        state_log_priors, tm, initial_log_probs,
        breakpoint_mask=breakpoint_mask,
        breakpoint_transition_matrix=bp_tm,
    )

    if verbose:
        segments = _extract_segments(path)
        seg_str = " -> ".join(
            f"CN{state}x{end - start}" for start, end, state in segments
        )
        conf_pass = "PASS" if confidence > confidence_threshold else "FAIL"
        _vlog(f"  [VIT] {sample_id}  {locus.cluster}  "
              f"{len(all_bins)} bins  ploidy={ploidy}  "
              f"confidence={confidence:+.4f} "
              f"(threshold={confidence_threshold:+.4f}) [{conf_pass}]")
        _vlog(f"    Segments: {seg_str}")

    # Map Viterbi path back to per-interval bin states
    local_idx_for_bin = {b: i for i, b in enumerate(all_bins)}
    interval_bin_states: Dict[str, np.ndarray] = {}
    for interval_name, bin_idx in interval_bin_arrays.items():
        local_indices = [local_idx_for_bin[int(b)] for b in bin_idx
                         if int(b) in local_idx_for_bin]
        if local_indices:
            interval_bin_states[interval_name] = path[local_indices]
        else:
            interval_bin_states[interval_name] = np.array([], dtype=int)

    all_interval_names = set(interval_bin_arrays.keys())
    flank_names = {name for _, _, name in locus.get_flanking_regions()}
    body_interval_names = all_interval_names - flank_names

    calls = []
    for entry in locus.gd_entries:
        gd_id    = entry["GD_ID"]
        svtype   = entry["svtype"]
        gd_start = entry["start_GRCh38"]
        gd_end   = entry["end_GRCh38"]
        bp1      = entry["BP1"]
        bp2      = entry["BP2"]

        covered_tuples    = locus.get_intervals_between(bp1, bp2)
        covered_intervals = set(name for _, _, name in covered_tuples)
        if not covered_intervals:
            print(f"  WARNING [Viterbi]: {gd_id} ({svtype}) BP1={bp1}, BP2={bp2} "
                  f"does not resolve to any interval in cluster '{locus.cluster}' "
                  f"(breakpoints: {locus.breakpoint_names}) — entry skipped.")
            continue

        uncovered_body = body_interval_names - covered_intervals

        if verbose:
            _vlog(f"    Entry {gd_id} ({svtype}, BP1={bp1}, BP2={bp2}  ploidy={ploidy}):")
            _vlog(f"      Covered ({len(covered_intervals)} interval(s)): "
                  f"{sorted(covered_intervals)}")

        # --- Check covered intervals ---
        covered_ok     = True
        n_covered_bins = 0
        covered_cn_sum = 0.0
        for iv in covered_intervals:
            states = interval_bin_states.get(iv, np.array([], dtype=int))
            if len(states) == 0:
                if verbose:
                    _vlog(f"        {iv}: 0 bins -> FAIL (no bins found for this interval)")
                covered_ok = False
                break
            n_covered_bins += len(states)
            covered_cn_sum += float(states.sum())
            if svtype == "DEL":
                iv_pass = bool(np.all(states < ploidy))
                if verbose:
                    bad = states[states >= ploidy]
                    status = ("PASS" if iv_pass else
                              f"FAIL  ({len(bad)} bin(s) >= {ploidy}; "
                              f"offending CN states: {np.unique(bad).tolist()})")
                    _vlog(f"        {iv}: {len(states)} bins  "
                          f"CN=[{int(states.min())},{int(states.max())}]  "
                          f"all<{ploidy} -> {status}")
                if not iv_pass:
                    covered_ok = False
                    break
            if svtype == "DUP":
                iv_pass = bool(np.all(states > ploidy))
                if verbose:
                    bad = states[states <= ploidy]
                    status = ("PASS" if iv_pass else
                              f"FAIL  ({len(bad)} bin(s) <= {ploidy}; "
                              f"offending CN states: {np.unique(bad).tolist()})")
                    _vlog(f"        {iv}: {len(states)} bins  "
                          f"CN=[{int(states.min())},{int(states.max())}]  "
                          f"all>{ploidy} -> {status}")
                if not iv_pass:
                    covered_ok = False
                    break

        # --- Check flanks and uncovered body intervals ---
        # At least flank_coverage_threshold fraction of bins must be at the
        # reference CN; remaining bins may be off (e.g. noise at boundaries).
        flanks_ok = True
        if verbose:
            flanks_to_check = sorted(uncovered_body | flank_names)
            _vlog(f"      Flanks / uncovered ({len(flanks_to_check)} region(s)): "
                  f"{flanks_to_check}  (coverage threshold={flank_coverage_threshold:.0%})")
        for iv in (uncovered_body | flank_names):
            states = interval_bin_states.get(iv, np.array([], dtype=int))
            if len(states) == 0:
                if verbose:
                    _vlog(f"        {iv}: 0 bins -> skip (missing flanks OK)")
                continue  # missing flanks are OK
            if svtype == "DEL":
                ok_frac = float(np.sum(states >= ploidy)) / len(states)
                iv_pass = ok_frac >= flank_coverage_threshold
                if verbose:
                    bad = states[states < ploidy]
                    status = (f"PASS  ({ok_frac:.0%} >= {flank_coverage_threshold:.0%})" if iv_pass else
                              f"FAIL  ({ok_frac:.0%} < {flank_coverage_threshold:.0%}; "
                              f"{len(bad)} offending bin(s): CN {np.unique(bad).tolist()})")
                    _vlog(f"        {iv}: {len(states)} bins  "
                          f"CN=[{int(states.min())},{int(states.max())}]  "
                          f">={ploidy} fraction={ok_frac:.0%} -> {status}")
            elif svtype == "DUP":
                ok_frac = float(np.sum(states <= ploidy)) / len(states)
                iv_pass = ok_frac >= flank_coverage_threshold
                if verbose:
                    bad = states[states > ploidy]
                    status = (f"PASS  ({ok_frac:.0%} >= {flank_coverage_threshold:.0%})" if iv_pass else
                              f"FAIL  ({ok_frac:.0%} < {flank_coverage_threshold:.0%}; "
                              f"{len(bad)} offending bin(s): CN {np.unique(bad).tolist()})")
                    _vlog(f"        {iv}: {len(states)} bins  "
                          f"CN=[{int(states.min())},{int(states.max())}]  "
                          f"<={ploidy} fraction={ok_frac:.0%} -> {status}")
            else:
                continue
            if not iv_pass:
                flanks_ok = False
                break

        is_cnv  = covered_ok and flanks_ok and (confidence > confidence_threshold)
        if verbose:
            cov_str  = "PASS" if covered_ok else "FAIL"
            flk_str  = "PASS" if flanks_ok else "FAIL"
            conf_str = "PASS" if confidence > confidence_threshold else "FAIL"
            result_str = "CARRIER" if is_cnv else "NOT_CARRIER"
            _vlog(f"      -> covered={cov_str}  flanks={flk_str}  "
                  f"confidence={conf_str}  [{result_str}]")
        mean_cn = covered_cn_sum / n_covered_bins if n_covered_bins > 0 else float(ploidy)
        n_bins  = sum(len(interval_bin_arrays.get(iv, [])) for iv in covered_intervals)

        calls.append({
            "GD_ID":                   gd_id,
            "cluster":                 locus.cluster,
            "chrom":                   locus.chrom,
            "start":                   gd_start,
            "end":                     gd_end,
            "svtype":                  svtype,
            "BP1":                     bp1,
            "BP2":                     bp2,
            "is_nahr":                 locus.is_nahr,
            "is_terminal":             locus.is_terminal,
            "log_prob_score":          confidence,
            "flanking_log_prob_score": np.nan,
            "is_carrier":              is_cnv,
            "is_spanning":             False,
            "intervals":               list(covered_intervals),
            "n_bins":                  n_bins,
            "mean_cn":                 mean_cn,
            "viterbi_covered_ok":         covered_ok,
            "viterbi_flanks_ok":           flanks_ok,
            "flank_coverage_threshold":    flank_coverage_threshold,
        })

    return calls


# =============================================================================
# Orchestrator
# =============================================================================


def call_cnvs_from_posteriors(
    cn_posteriors_df: pd.DataFrame,
    bin_mappings_df: pd.DataFrame,
    gd_table: GDTable,
    log_prob_threshold: float = -0.5,
    flanking_log_prob_threshold: float = -1.0,
    ploidy_df: Optional[pd.DataFrame] = None,
    verbose: bool = False,
    transition_matrix: Optional[np.ndarray] = None,
    viterbi_confidence_threshold: float = -1.0,
    viterbi_flank_coverage_threshold: float = 0.70,
    breakpoint_transition_matrix: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """
    Call CNVs from posterior probabilities.

    When *transition_matrix* is provided, the Viterbi segmentation strategy
    is used.  Otherwise the original mean log-probability strategy is used.

    Args:
        cn_posteriors_df: DataFrame with CN posterior probabilities per bin/sample
        bin_mappings_df: DataFrame with bin-to-interval mappings
        gd_table: GDTable with locus definitions
        log_prob_threshold: Minimum log probability score to call a CNV
            (mean-log-prob strategy only)
        flanking_log_prob_threshold: Minimum log probability score in flanking
            regions to classify as spanning (mean-log-prob strategy only)
        ploidy_df: Optional DataFrame with columns (sample, contig, ploidy).
            If None, ploidy=2 is assumed for all sample/contig pairs.
        verbose: If True, print per-sample log probability scores for every
            GD entry at every locus, including flanking region scores.
        transition_matrix: Optional (n_states, n_states) CN-state transition
            probability matrix.  When provided, the Viterbi strategy is used.
        viterbi_confidence_threshold: Minimum mean per-bin log-posterior
            along the Viterbi path to make a carrier call (Viterbi strategy only).
        viterbi_flank_coverage_threshold: Fraction of bins in each flanking /
            uncovered interval that must be at the reference CN for the flank
            check to pass (Viterbi strategy only).  Default 0.70.
        breakpoint_transition_matrix: Optional (n_states, n_states) transition
            matrix applied at known recurrent breakpoint boundaries (start and
            end of each breakpoint's SD-block range) during Viterbi.  Should
            have lower diagonal values than *transition_matrix* to encourage
            CN state changes at those positions.  Viterbi strategy only.

    Returns:
        DataFrame with CNV calls for all samples and loci
    """
    use_viterbi = transition_matrix is not None

    print("\n" + "=" * 80)
    print("CALLING CNVs FROM POSTERIORS")
    if use_viterbi:
        bp_str = " + breakpoint matrix" if breakpoint_transition_matrix is not None else ""
        print(f"  Strategy: Viterbi segmentation{bp_str}  "
              f"(confidence threshold={viterbi_confidence_threshold}, "
              f"flank coverage threshold={viterbi_flank_coverage_threshold:.0%})")
    else:
        print(f"  Strategy: Mean log-probability  "
              f"(threshold={log_prob_threshold})")
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

    # Pre-extract numpy arrays for fast access.
    print("  Organizing data for fast access...")
    prob_cols = sorted(c for c in cn_posteriors_df.columns if c.startswith("prob_cn_"))
    n_states = len(prob_cols)
    prob_3d = np.empty((n_samples, n_bins, n_states))
    depth_2d = np.empty((n_samples, n_bins))
    for s_idx, sample_id in enumerate(sample_ids):
        mask = cn_posteriors_df["sample"] == sample_id
        sample_rows = cn_posteriors_df.loc[mask]
        prob_3d[s_idx] = sample_rows[prob_cols].values
        depth_2d[s_idx] = sample_rows["depth"].values
    print(f"    Extracted {n_samples} x {n_bins} x {n_states} probability array")

    # Build a per-bin coordinate lookup used by viterbi_call_gd_cnv to identify
    # which transitions cross known breakpoint boundaries.
    bin_coords_by_idx: Dict[int, Tuple[int, int]] = dict(zip(
        bin_mappings_df["array_idx"].astype(int),
        zip(
            bin_mappings_df["start"].astype(int),
            bin_mappings_df["end"].astype(int),
        ),
    ))

    for cluster, locus in gd_table.loci.items():
        print(f"\nCalling CNVs for locus: {cluster}")

        interval_bins = get_locus_interval_bins(bin_mappings_df, cluster)

        # Breakpoint-range bins (inside SD-block ranges, not in any body
        # interval or flank) carry unreliable depth signal and must be
        # excluded from CNV calling.  They should already be absent from
        # the mappings file but this is a safety guard.
        bp_masked = interval_bins.pop("breakpoint_ranges", [])
        if bp_masked:
            print(f"  Masking {len(bp_masked)} breakpoint-range bin(s)")

        for interval_name, bins in interval_bins.items():
            print(f"  {interval_name}: {len(bins)} bins")

        interval_bin_arrays = {
            name: np.array(bins, dtype=int)
            for name, bins in interval_bins.items()
        }

        for s_idx, sample_id in enumerate(sample_ids):
            sample_ploidy = ploidy_lookup.get((str(sample_id), locus.chrom), 2)

            if verbose:
                strategy = "VIT" if use_viterbi else "MLP"
                _vlog(f"\n  [{strategy}] Sample: {sample_id}  "
                      f"({locus.chrom}, ploidy={sample_ploidy})")

            if use_viterbi:
                # ---- Viterbi strategy ----
                calls = viterbi_call_gd_cnv(
                    locus,
                    prob_3d[s_idx],  # (n_bins, n_states)
                    transition_matrix,
                    interval_bin_arrays,
                    ploidy=sample_ploidy,
                    confidence_threshold=viterbi_confidence_threshold,
                    flank_coverage_threshold=viterbi_flank_coverage_threshold,
                    verbose=verbose,
                    sample_id=str(sample_id),
                    breakpoint_transition_matrix=breakpoint_transition_matrix,
                    bin_coords=bin_coords_by_idx,
                )
                best_by_svtype = determine_best_breakpoints(
                    locus,
                    # Build interval_stats for best-breakpoint scoring
                    {name: {"n_bins": len(idx),
                            "cn_probs": prob_3d[s_idx, idx].mean(axis=0)
                            if len(idx) > 0 else np.zeros(n_states)}
                     for name, idx in interval_bin_arrays.items()},
                    calls, ploidy=sample_ploidy,
                )
            else:
                # ---- Mean log-probability strategy ----
                interval_stats = {}
                for interval_name, bin_idx in interval_bin_arrays.items():
                    if len(bin_idx) == 0:
                        interval_stats[interval_name] = {
                            "n_bins": 0,
                            "cn_probs": np.zeros(n_states),
                        }
                    else:
                        mean_probs = prob_3d[s_idx, bin_idx].mean(axis=0)
                        interval_stats[interval_name] = {
                            "n_bins": len(bin_idx),
                            "cn_probs": mean_probs,
                        }

                calls = call_gd_cnv(
                    locus, interval_stats, log_prob_threshold,
                    flanking_log_prob_threshold,
                    ploidy=sample_ploidy,
                    verbose=verbose,
                    sample_id=str(sample_id),
                )
                best_by_svtype = determine_best_breakpoints(
                    locus, interval_stats, calls, ploidy=sample_ploidy)

            if verbose:
                strategy = "VIT" if use_viterbi else "MLP"
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
                    extra = ""
                    if use_viterbi:
                        cov = "Y" if call.get("viterbi_covered_ok") else "N"
                        flk = "Y" if call.get("viterbi_flanks_ok") else "N"
                        extra = f"  cov={cov} flk={flk}"
                    _vlog(
                        f"    [{strategy}] {sample_id:30s}  "
                        f"{call['GD_ID']:25s}  "
                        f"{call['svtype']:4s}  ploidy={sample_ploidy}  "
                        f"score={call['log_prob_score']:+.4f}  "
                        f"flank={flank_str}  "
                        f"intervals={','.join(call['intervals'])}"
                        f"{extra}{tag}"
                    )

            for call in calls:
                svtype = call["svtype"]
                best_gd_for_svtype = best_by_svtype.get(svtype)

                covered_bin_indices = []
                for interval_name in call["intervals"]:
                    if interval_name in interval_bin_arrays:
                        covered_bin_indices.extend(interval_bin_arrays[interval_name])

                if len(covered_bin_indices) > 0:
                    mean_depth = float(depth_2d[s_idx, covered_bin_indices].mean())
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
                    "mean_cn": call.get("mean_cn", np.nan),
                    "is_carrier": call["is_carrier"],
                    "is_best_match": (
                        call["GD_ID"] == best_gd_for_svtype
                        if best_gd_for_svtype else False
                    ),
                    "log_prob_score": call["log_prob_score"],
                }
                all_results.append(result)

    return pd.DataFrame(all_results)


# =============================================================================
# Main
# =============================================================================


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
        "--log-prob-threshold",
        type=float,
        default=-0.3,
        help="Minimum log probability score to call a CNV",
    )
    parser.add_argument(
        "--flanking-log-prob-threshold",
        type=float,
        default=-1.0,
        help="Minimum log probability score in flanking regions to classify "
             "a call as spanning",
    )
    parser.add_argument(
        "--transition-matrix",
        required=False,
        help="CN-state transition probability matrix (TSV).  When provided, "
             "the Viterbi segmentation strategy is used instead of the "
             "mean log-probability strategy.  The file should be a square "
             "K×K matrix where entry (i,j) is the probability of "
             "transitioning from CN state i to CN state j between adjacent "
             "bins.  Rows will be automatically normalized to sum to 1.",
    )
    parser.add_argument(
        "--breakpoint-transition-matrix",
        required=False,
        help="CN-state transition probability matrix (TSV) applied "
             "specifically at known recurrent breakpoint boundaries — both "
             "the start and end coordinate of each breakpoint's SD-block "
             "range.  Should have lower diagonal values than "
             "--transition-matrix to encourage CN state changes at these "
             "positions.  Same TSV format as --transition-matrix.  "
             "Only used when --transition-matrix is also set (Viterbi only).",
    )
    parser.add_argument(
        "--viterbi-confidence-threshold",
        type=float,
        default=-1.0,
        help="Minimum mean per-bin log-posterior along the Viterbi path to "
             "call a carrier (Viterbi strategy only).  More negative values "
             "are more permissive.",
    )
    parser.add_argument(
        "--viterbi-flank-coverage-threshold",
        type=float,
        default=0.70,
        help="Fraction of bins in each flanking / uncovered interval that must "
             "be at the reference copy number for the flank check to pass "
             "(Viterbi strategy only).  1.0 = all bins must pass (strict); "
             "0.0 = disabled.  Default: 0.70.",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print detailed per-sample log probability scores for all GD "
             "entries at every locus, including flanking region scores.",
    )
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    log_fh = setup_logging(args.output_dir)

    print(f"Output directory: {args.output_dir}")

    # Load data
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

    # Load ploidy table if provided
    ploidy_df = None
    if args.ploidy_table:
        print(f"  Loading ploidy table: {args.ploidy_table}")
        ploidy_df = pd.read_csv(args.ploidy_table, sep="\t")
        print(f"    {len(ploidy_df)} sample/contig ploidy records")

    # Load transition matrix if provided
    transition_matrix = None
    if args.transition_matrix:
        print(f"\n  Loading transition matrix: {args.transition_matrix}")
        transition_matrix = load_transition_matrix(args.transition_matrix)

    # Load breakpoint-specific transition matrix if provided
    breakpoint_transition_matrix = None
    if args.breakpoint_transition_matrix:
        if transition_matrix is None:
            print("  WARNING: --breakpoint-transition-matrix requires "
                  "--transition-matrix to be set; ignoring")
        else:
            print(f"  Loading breakpoint transition matrix: "
                  f"{args.breakpoint_transition_matrix}")
            breakpoint_transition_matrix = load_transition_matrix(
                args.breakpoint_transition_matrix)

    # Call CNVs from posteriors
    calls_df = call_cnvs_from_posteriors(
        cn_posteriors_df,
        bin_mappings_df,
        gd_table,
        log_prob_threshold=args.log_prob_threshold,
        flanking_log_prob_threshold=args.flanking_log_prob_threshold,
        ploidy_df=ploidy_df,
        verbose=args.verbose,
        transition_matrix=transition_matrix,
        viterbi_confidence_threshold=args.viterbi_confidence_threshold,
        viterbi_flank_coverage_threshold=args.viterbi_flank_coverage_threshold,
        breakpoint_transition_matrix=breakpoint_transition_matrix,
    )

    # Save calls
    output_file = os.path.join(args.output_dir, "gd_cnv_calls.tsv.gz")
    calls_df.to_csv(output_file, sep="\t", index=False, compression="gzip")
    print(f"\n  Saved calls to: {output_file}")
    print(f"    {len(calls_df)} call records")

    # Print carrier summary
    carriers = calls_df[calls_df["is_carrier"]]
    n_carriers = carriers["sample"].nunique()
    n_sites = carriers["GD_ID"].nunique()
    print(f"\n  {n_carriers} carrier samples across {n_sites} GD sites")

    print("\n" + "=" * 80)
    print("Calling complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
