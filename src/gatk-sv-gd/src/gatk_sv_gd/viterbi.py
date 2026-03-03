"""
Viterbi HMM Segmentation for GD CNV Calling

Implements the Viterbi algorithm on per-bin CN posterior probabilities to
produce smooth copy-number segmentations, then checks each GD entry's
breakpoint pattern against the resulting path to call carriers.

This module is used by :mod:`gatk_sv_gd.call` when the ``--transition-matrix``
option is provided.
"""

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from gatk_sv_gd import _util
from gatk_sv_gd.models import GDLocus


# =============================================================================
# Transition matrix I/O
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


# =============================================================================
# Viterbi algorithm
# =============================================================================


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


# =============================================================================
# Viterbi-based GD CNV calling
# =============================================================================


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
            _util.vlog(f"  [VIT] {sample_id}  {locus.cluster}: no bins -> skip")
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
                _util.vlog(
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
        _util.vlog(f"  [VIT] {sample_id}  {locus.cluster}  "
                   f"{len(all_bins)} bins  ploidy={ploidy}  "
                   f"confidence={confidence:+.4f} "
                   f"(threshold={confidence_threshold:+.4f}) [{conf_pass}]")
        _util.vlog(f"    Segments: {seg_str}")

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
            _util.vlog(f"    Entry {gd_id} ({svtype}, BP1={bp1}, BP2={bp2}  ploidy={ploidy}):")
            _util.vlog(f"      Covered ({len(covered_intervals)} interval(s)): "
                       f"{sorted(covered_intervals)}")

        # --- Check covered intervals ---
        covered_ok     = True
        n_covered_bins = 0
        covered_cn_sum = 0.0
        for iv in covered_intervals:
            states = interval_bin_states.get(iv, np.array([], dtype=int))
            if len(states) == 0:
                if verbose:
                    _util.vlog(f"        {iv}: 0 bins -> FAIL (no bins found for this interval)")
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
                    _util.vlog(f"        {iv}: {len(states)} bins  "
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
                    _util.vlog(f"        {iv}: {len(states)} bins  "
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
            _util.vlog(f"      Flanks / uncovered ({len(flanks_to_check)} region(s)): "
                       f"{flanks_to_check}  (coverage threshold={flank_coverage_threshold:.0%})")
        for iv in (uncovered_body | flank_names):
            states = interval_bin_states.get(iv, np.array([], dtype=int))
            if len(states) == 0:
                if verbose:
                    _util.vlog(f"        {iv}: 0 bins -> skip (missing flanks OK)")
                continue  # missing flanks are OK
            if svtype == "DEL":
                ok_frac = float(np.sum(states >= ploidy)) / len(states)
                iv_pass = ok_frac >= flank_coverage_threshold
                if verbose:
                    bad = states[states < ploidy]
                    status = (f"PASS  ({ok_frac:.0%} >= {flank_coverage_threshold:.0%})" if iv_pass else
                              f"FAIL  ({ok_frac:.0%} < {flank_coverage_threshold:.0%}; "
                              f"{len(bad)} offending bin(s): CN {np.unique(bad).tolist()})")
                    _util.vlog(f"        {iv}: {len(states)} bins  "
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
                    _util.vlog(f"        {iv}: {len(states)} bins  "
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
            _util.vlog(f"      -> covered={cov_str}  flanks={flk_str}  "
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
