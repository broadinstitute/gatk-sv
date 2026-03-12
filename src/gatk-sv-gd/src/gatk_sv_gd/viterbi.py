"""
Viterbi HMM Segmentation for GD CNV Calling

Implements the Viterbi algorithm on per-bin CN posterior probabilities to
produce smooth copy-number segmentations, then checks each GD entry's
breakpoint pattern against the resulting path to call carriers.

For **diploid** samples the hidden state is a *pair* of per-haplotype copy
numbers ``(h1, h2)`` with ``h1 <= h2`` (canonical ordering).  The observed
total CN posterior ``p(CN=k | data)`` acts as the prior on each pair whose
``h1 + h2 == k``.  Transitions factorise as the product of independent
per-haplotype transitions.  Any bias toward or against specific
haplotype states should be encoded in the transition matrix itself.

For **haploid** samples the algorithm is unchanged.

Samples with **ploidy > 2** are skipped.

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
# Diploid pair-state helpers
# =============================================================================

# Maximum per-haplotype copy number.  For a diploid sample the hidden state
# is (h1, h2) with 0 <= h1 <= h2 <= _MAX_HAP_CN.  Total CN = h1 + h2.
_MAX_HAP_CN = 3


def _build_pair_states(max_hap_cn: int = _MAX_HAP_CN) -> List[Tuple[int, int]]:
    """Return the canonical list of (h1, h2) pairs with h1 <= h2."""
    return [
        (h1, h2)
        for h1 in range(max_hap_cn + 1)
        for h2 in range(h1, max_hap_cn + 1)
    ]


def _build_diploid_obs_log_priors(
    total_cn_log_priors: np.ndarray,
    pair_states: List[Tuple[int, int]],
) -> np.ndarray:
    """Derive per-bin log priors for each (h1, h2) pair from total CN posteriors.

    For each pair (h1, h2) with total = h1 + h2, the observation prior is
    simply ``log p(total_CN = h1+h2 | data)``.  All pairs that sum to the
    same total share the same observation prior — the transition matrix
    alone determines which decomposition is preferred.

    No degeneracy bonus is applied for heterozygous pairs (h1 != h2).
    Although there are two physical orderings of such pairs, adding a
    per-bin log(2) bonus would overwhelm the transition penalty after
    only a few bins and systematically favour heterozygous states like
    (0, 2) over the normal (1, 1).

    Args:
        total_cn_log_priors: (T, K_total) array of log p(CN=k | data) for
            each of T bins.
        pair_states: list of (h1, h2) pairs from ``_build_pair_states()``.

    Returns:
        (T, n_pairs) array of log observation priors for each pair state.
    """
    T, K_total = total_cn_log_priors.shape
    n_pairs = len(pair_states)
    obs = np.full((T, n_pairs), -1e30)

    for p_idx, (h1, h2) in enumerate(pair_states):
        total = h1 + h2
        if total >= K_total:
            continue  # total CN outside the model's state range
        obs[:, p_idx] = total_cn_log_priors[:, total]

    return obs


def _equalize_hap_diagonals(hap_trans: np.ndarray) -> np.ndarray:
    """Equalize diagonal (self-transition) entries of a per-haplotype matrix.

    When independent per-haplotype transitions are combined into a diploid
    pair-state matrix, rows with *lower* diagonal values produce less sticky
    pair states.  This creates an artifact: state 1 (per-haplotype
    reference) typically has the most off-diagonal weight — easy transitions
    to CN 0, 2, 3, etc. — so its diagonal is lower after normalization.
    The reference pair ``(1, 1)`` therefore becomes less sticky than
    non-reference pairs like ``(0, 2)`` whose constituent states 0 and 2
    each have fewer exit paths.

    This function rescales each row so that **all diagonal entries equal the
    maximum diagonal** in the matrix.  Relative off-diagonal ratios are
    preserved::

        new_off_diag[i, j] = old_off_diag[i, j] × (1 − d_max) / (1 − d_i)

    where ``d_i`` is the original diagonal and ``d_max`` is the target.

    Returns:
        A copy of *hap_trans* with equalized diagonals.
    """
    n = hap_trans.shape[0]
    diags = np.diag(hap_trans).copy()
    d_max = diags.max()

    result = hap_trans.copy()
    for i in range(n):
        if diags[i] >= d_max - 1e-12:
            continue
        off_diag_sum = 1.0 - diags[i]
        new_off_diag_sum = 1.0 - d_max
        if off_diag_sum > 1e-30:
            scale = new_off_diag_sum / off_diag_sum
            result[i] = hap_trans[i] * scale
            result[i, i] = d_max
    return result


def _build_diploid_transition_matrix(
    hap_trans: np.ndarray,
    pair_states: List[Tuple[int, int]],
) -> np.ndarray:
    """Build a pair-state transition matrix from a per-haplotype matrix.

    Before combining, the per-haplotype diagonals are equalized via
    :func:`_equalize_hap_diagonals` so that no pair state is artificially
    stickier than the reference pair ``(1, 1)``.

    The two haplotypes then transition independently::

        P((h1', h2') | (h1, h2)) = P_hap(h1' | h1) * P_hap(h2' | h2)

    When (h1, h2) → (h1', h2') requires remapping to canonical order
    (h1' <= h2'), we sum over both orderings.

    Args:
        hap_trans: (K_hap, K_hap) per-haplotype transition matrix.
        pair_states: list of (h1, h2) canonical pairs.

    Returns:
        (n_pairs, n_pairs) transition matrix over pair states.
    """
    # Equalize per-haplotype diagonals so that the reference pair (1,1)
    # is not disadvantaged relative to pairs like (0,2).
    hap_trans = _equalize_hap_diagonals(hap_trans)

    n_pairs = len(pair_states)
    pair_idx = {pair: i for i, pair in enumerate(pair_states)}
    trans = np.zeros((n_pairs, n_pairs))

    max_cn = hap_trans.shape[0]
    for i, (a1, a2) in enumerate(pair_states):
        if a1 >= max_cn or a2 >= max_cn:
            trans[i, i] = 1.0
            continue
        for b1 in range(max_cn):
            for b2 in range(max_cn):
                # Canonical ordering
                lo, hi = (b1, b2) if b1 <= b2 else (b2, b1)
                j = pair_idx.get((lo, hi))
                if j is None:
                    continue
                prob = hap_trans[a1, b1] * hap_trans[a2, b2]
                trans[i, j] += prob

    # Normalize rows
    row_sums = trans.sum(axis=1, keepdims=True)
    row_sums = np.maximum(row_sums, 1e-30)
    trans = trans / row_sums
    return trans


def _build_diploid_initial_log_probs(
    pair_states: List[Tuple[int, int]],
) -> np.ndarray:
    """Initial log probabilities for pair states.

    Strong prior on (1, 1) — the normal diploid state.  Small uniform
    probability on all other pairs.
    """
    n_pairs = len(pair_states)
    log_probs = np.full(n_pairs, np.log(1e-6))
    for i, (h1, h2) in enumerate(pair_states):
        if h1 == 1 and h2 == 1:
            log_probs[i] = np.log(1.0 - 1e-6 * (n_pairs - 1))
    return log_probs


def _run_diploid_viterbi(
    total_cn_log_priors: np.ndarray,
    hap_transition_matrix: np.ndarray,
    breakpoint_mask: Optional[np.ndarray],
    hap_breakpoint_transition_matrix: Optional[np.ndarray],
) -> Tuple[np.ndarray, np.ndarray, float, List[Tuple[int, int]]]:
    """Run the Viterbi algorithm over diploid pair states.

    Args:
        total_cn_log_priors: (T, K_total) log posteriors from the Pyro model.
        hap_transition_matrix: (K_hap, K_hap) per-haplotype transition matrix.
        breakpoint_mask: optional (T-1,) boolean mask.
        hap_breakpoint_transition_matrix: optional per-haplotype breakpoint
            transition matrix.

    Returns:
        Tuple of:
          - hap1_path: (T,) int array of haplotype-1 CN states.
          - hap2_path: (T,) int array of haplotype-2 CN states.
          - confidence: float mean per-bin log-prior along the path.
          - pair_states: the list of (h1, h2) pairs used.
    """
    pair_states = _build_pair_states()

    # Observation priors
    obs = _build_diploid_obs_log_priors(
        total_cn_log_priors, pair_states,
    )

    # Transition matrices
    pair_trans = _build_diploid_transition_matrix(
        hap_transition_matrix, pair_states,
    )
    pair_bp_trans = None
    if hap_breakpoint_transition_matrix is not None:
        pair_bp_trans = _build_diploid_transition_matrix(
            hap_breakpoint_transition_matrix, pair_states,
        )

    # Initial probs
    initial = _build_diploid_initial_log_probs(
        pair_states,
    )

    # Run standard Viterbi on pair states
    pair_path, confidence = run_viterbi(
        obs, pair_trans, initial,
        breakpoint_mask=breakpoint_mask,
        breakpoint_transition_matrix=pair_bp_trans,
    )

    # Decompose into per-haplotype paths
    T = len(pair_path)
    hap1 = np.empty(T, dtype=int)
    hap2 = np.empty(T, dtype=int)
    for t in range(T):
        h1, h2 = pair_states[pair_path[t]]
        hap1[t] = h1
        hap2[t] = h2

    return hap1, hap2, confidence, pair_states


# =============================================================================
# Category-segment builder
# =============================================================================


def _build_category_segments(
    path: np.ndarray,
    all_bins: List[int],
    bin_coords: Dict[int, Tuple[int, int]],
    ploidy: int,
) -> List[dict]:
    """Partition the Viterbi path into contiguous genomic segments by exact
    CN state.

    Adjacent bins with the **same integer CN state** are merged into a
    single segment.  Each segment is labelled with a category relative to
    *ploidy*:

    - ``DEL``: CN < ploidy
    - ``REF``: CN == ploidy
    - ``DUP``: CN > ploidy

    Returns:
        List of dicts with keys ``category``, ``start``, ``end``,
        ``cn_state`` (int), ``n_bins``.
    """
    if len(path) == 0:
        return []

    def _cat(cn: int) -> str:
        if cn < ploidy:
            return "DEL"
        elif cn > ploidy:
            return "DUP"
        return "REF"

    segments: List[dict] = []
    first_coords = bin_coords.get(all_bins[0])
    if first_coords is None:
        return []

    seg_start = first_coords[0]
    seg_end = first_coords[1]
    seg_cn = int(path[0])
    seg_n = 1

    for i in range(1, len(path)):
        coords = bin_coords.get(all_bins[i])
        if coords is None:
            continue
        cn = int(path[i])
        if cn == seg_cn:
            seg_end = coords[1]
            seg_n += 1
        else:
            segments.append({
                "category": _cat(seg_cn),
                "start": seg_start,
                "end": seg_end,
                "cn_state": seg_cn,
                "n_bins": seg_n,
            })
            seg_start = coords[0]
            seg_end = coords[1]
            seg_cn = cn
            seg_n = 1

    segments.append({
        "category": _cat(seg_cn),
        "start": seg_start,
        "end": seg_end,
        "cn_state": seg_cn,
        "n_bins": seg_n,
    })
    return segments


def _combine_hap_match_segments(
    hap1_segs: List[dict],
    hap2_segs: List[dict],
) -> List[dict]:
    """Merge per-haplotype non-REF segments for GD entry matching.

    For diploid calling we want to match GD entries against per-haplotype
    segments rather than total-CN segments, so that overlapping events on
    different haplotypes don't cancel each other out.

    The logic:
      - Collect all DEL segments from either haplotype.
      - Collect all DUP segments from either haplotype.
      - Merge overlapping/adjacent segments of the same svtype.
      - Return the merged list.  ``mean_cn`` is set to the per-haplotype
        value (relative to ploidy=1 per haplotype).

    This lets e.g. a DUP spanning A→E on haplotype 1 be matched even when
    haplotype 2 has a DEL in the D-E region that lowers total CN.
    """
    combined: List[dict] = []
    for svtype in ("DEL", "DUP"):
        segs = sorted(
            [s for s in hap1_segs + hap2_segs if s["category"] == svtype],
            key=lambda s: s["start"],
        )
        if not segs:
            continue
        # Merge overlapping/adjacent
        merged_start = segs[0]["start"]
        merged_end = segs[0]["end"]
        merged_cn = segs[0]["cn_state"]
        merged_n = segs[0]["n_bins"]
        for s in segs[1:]:
            if s["start"] <= merged_end:
                merged_end = max(merged_end, s["end"])
                # Keep the cn_state from whichever segment has more bins
                if s["n_bins"] > merged_n:
                    merged_cn = s["cn_state"]
                merged_n += s["n_bins"]
            else:
                combined.append({
                    "category": svtype,
                    "start": merged_start,
                    "end": merged_end,
                    "cn_state": merged_cn,
                    "n_bins": merged_n,
                })
                merged_start = s["start"]
                merged_end = s["end"]
                merged_cn = s["cn_state"]
                merged_n = s["n_bins"]
        combined.append({
            "category": svtype,
            "start": merged_start,
            "end": merged_end,
            "cn_state": merged_cn,
            "n_bins": merged_n,
        })
    return combined


# =============================================================================
# Viterbi-based GD CNV calling
# =============================================================================


def viterbi_call_gd_cnv(
    locus: GDLocus,
    sample_probs: np.ndarray,
    transition_matrix: np.ndarray,
    interval_bin_arrays: Dict[str, np.ndarray],
    ploidy: int = 2,
    reciprocal_overlap_threshold: float = 0.90,
    verbose: bool = False,
    sample_id: str = "",
    breakpoint_transition_matrix: Optional[np.ndarray] = None,
    bin_coords: Optional[Dict[int, Tuple[int, int]]] = None,
    all_cluster_bins: Optional[List[int]] = None,
) -> List[dict]:
    """Call GD CNVs using Viterbi segmentation of CN posteriors.

    The Viterbi algorithm is run **once** per locus/sample to produce a
    smooth copy-number segmentation.  The resulting path is then
    partitioned into contiguous regions of altered copy number relative to
    *ploidy*:

    - **DEL** segments: CN < ploidy
    - **DUP** segments: CN > ploidy

    Each altered-CN segment is tested for **reciprocal overlap** against
    every GD entry of matching svtype.  A carrier call is made when

        min(overlap / segment_length, overlap / entry_length)
            >= reciprocal_overlap_threshold

    Multiple calls per locus are possible (though uncommon).

    Args:
        locus: GDLocus object.
        sample_probs: (n_bins, n_states) posterior probabilities for one
            sample across all bins at this locus (in array_idx order).
            These are p(CN=k | data) from the Pyro model.
        transition_matrix: (n_states, n_states) transition matrix.  High
            diagonal = sticky segmentation.
        interval_bin_arrays: Dict mapping interval name to array of bin
            array_idx values.  Used for matching GD entries.
        ploidy: Expected reference copy number for this sample/chromosome.
        reciprocal_overlap_threshold: Minimum reciprocal overlap between a
            Viterbi segment and a GD entry to call a carrier.
            Default 0.90 (90 %).
        verbose: If True, emit diagnostic messages via ``_util.vlog``.
        sample_id: Sample identifier for diagnostic logging.
        breakpoint_transition_matrix: Optional (n_states, n_states) transition
            matrix applied at known recurrent breakpoint boundaries (both the
            start and end coordinate of each breakpoint's SD-block range).
            Should have lower diagonal values than *transition_matrix* so that
            CN-state changes are more favoured at these positions.  Requires
            *bin_coords* to locate boundaries.
        bin_coords: Mapping from array_idx (int) to (genomic_start, genomic_end)
            for each bin.  Required for computing reciprocal overlap with GD
            entries and for the breakpoint mask.
        all_cluster_bins: Optional pre-computed list of ALL bin array_idx
            values for this cluster (excluding breakpoint-range bins).  When
            provided, the Viterbi runs over this full set rather than only
            the bins found in *interval_bin_arrays*, ensuring flanks are
            always included even if they are not assigned to a named interval.
            If *None*, falls back to the union of *interval_bin_arrays* values.

    Returns:
        Tuple of:
          - List of call dicts (one per GD entry).
          - List of (genomic_start, genomic_end, cn_state, category,
            haplotype) tuples, one per contiguous CN-state segment.  Empty
            when *bin_coords* is None or no bins are available.
    """
    # Collect all bin indices for this locus.
    # Prefer the caller-supplied complete set if available, since it
    # guarantees flanks are included even when they are not assigned to a
    # named interval in the bin-mappings file.
    if all_cluster_bins is not None:
        all_bins_set = set(int(b) for b in all_cluster_bins)
    else:
        all_bins_set = {int(b) for bins in interval_bin_arrays.values()
                        for b in bins}

    if len(all_bins_set) == 0:
        if verbose:
            _util.vlog(f"  [VIT] {sample_id}  {locus.cluster}: no bins -> skip")
        return [], []

    if bin_coords is None:
        if verbose:
            _util.vlog(f"  [VIT] {sample_id}  {locus.cluster}: "
                       f"no bin_coords -> skip")
        return [], []

    # Sort bins by **genomic start coordinate** (not array_idx).
    # Flank bins may have array_idx values that don't match their genomic
    # position (e.g. left_flank is genomically first but has high array_idx).
    # Sorting by coordinate ensures the Viterbi processes spatially adjacent
    # bins sequentially, and that _build_category_segments produces segments
    # in proper genomic order.
    all_bins = sorted(
        all_bins_set,
        key=lambda idx: bin_coords.get(idx, (float("inf"), float("inf")))[0],
    )

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

    # Subset breakpoint matrix to n_states if provided
    bp_tm = (
        breakpoint_transition_matrix[:n_states, :n_states]
        if breakpoint_transition_matrix is not None else None
    )

    # ==========================================================================
    # Dispatch: diploid pair-state model  vs  haploid (unchanged)
    # ==========================================================================
    if ploidy == 2:
        # ----- Diploid: run Viterbi on (h1, h2) pair states -----
        hap1_path, hap2_path, confidence, pair_states = _run_diploid_viterbi(
            state_log_priors,
            tm,
            breakpoint_mask,
            bp_tm,
        )
        total_path = hap1_path + hap2_path  # total CN per bin

        if verbose:
            # Log pair-state segments
            pair_segs = _extract_segments(total_path)
            seg_str = " -> ".join(
                f"CN{state}x{end - start}" for start, end, state in pair_segs
            )
            _util.vlog(f"  [VIT-DIP] {sample_id}  {locus.cluster}  "
                       f"{len(all_bins)} bins  ploidy={ploidy}  "
                       f"confidence={confidence:+.4f}")
            _util.vlog(f"    Total CN segments: {seg_str}")
            h1_segs = _extract_segments(hap1_path)
            h2_segs = _extract_segments(hap2_path)
            _util.vlog(f"    Hap1: {' -> '.join(f'CN{s}x{e-b}' for b,e,s in h1_segs)}")
            _util.vlog(f"    Hap2: {' -> '.join(f'CN{s}x{e-b}' for b,e,s in h2_segs)}")

        # Build per-haplotype category segments (ploidy_per_hap = 1)
        hap1_cat_segs = _build_category_segments(
            hap1_path, all_bins, bin_coords, ploidy=1,
        )
        hap2_cat_segs = _build_category_segments(
            hap2_path, all_bins, bin_coords, ploidy=1,
        )
        # Also build total-CN category segments for GD entry matching
        total_cat_segs = _build_category_segments(
            total_path, all_bins, bin_coords, ploidy=ploidy,
        )

    elif ploidy == 1:
        # ----- Haploid: unchanged single-chain Viterbi -----
        initial_log_probs = np.full(n_states, np.log(1e-6))
        if ploidy < n_states:
            initial_log_probs[ploidy] = np.log(1.0 - 1e-6 * (n_states - 1))

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
            _util.vlog(f"  [VIT] {sample_id}  {locus.cluster}  "
                       f"{len(all_bins)} bins  ploidy={ploidy}  "
                       f"confidence={confidence:+.4f}")
            _util.vlog(f"    Segments: {seg_str}")

        total_path = path
        hap1_path = path
        hap2_path = None
        hap1_cat_segs = _build_category_segments(path, all_bins, bin_coords, ploidy)
        hap2_cat_segs = []
        total_cat_segs = hap1_cat_segs

    else:
        # Ploidy > 2 is not supported — return empty results.
        if verbose:
            _util.vlog(f"  [VIT] {sample_id}  {locus.cluster}: "
                       f"ploidy={ploidy} > 2, skipping")
        return [], []

    if verbose:
        for seg in total_cat_segs:
            if seg["category"] != "REF":
                _util.vlog(
                    f"    {seg['category']} segment: "
                    f"{seg['start']:,}-{seg['end']:,} "
                    f"({seg['n_bins']} bins, CN={seg['cn_state']})"
                )

    # ----------------------------------------------------------------
    # Match each GD entry against category segments by reciprocal overlap.
    # For diploid we match per-haplotype segments: a DEL is any haplotype
    # dropping below reference (CN=1) and a DUP is any haplotype rising
    # above reference.
    # ----------------------------------------------------------------

    # Combine per-haplotype non-REF segments for matching.  For haploid
    # this is identical to total_cat_segs.
    match_segments = total_cat_segs if ploidy != 2 else _combine_hap_match_segments(
        hap1_cat_segs, hap2_cat_segs,
    )

    calls: List[dict] = []
    for entry in locus.gd_entries:
        gd_id    = entry["GD_ID"]
        svtype   = entry["svtype"]
        gd_start = entry["start_GRCh38"]
        gd_end   = entry["end_GRCh38"]
        bp1      = entry["BP1"]
        bp2      = entry["BP2"]

        covered_tuples = locus.get_intervals_between(bp1, bp2)
        covered_intervals = [name for _, _, name in covered_tuples]

        entry_len = gd_end - gd_start
        best_ro = 0.0
        best_seg = None
        if entry_len > 0:
            for seg in match_segments:
                if seg["category"] != svtype:
                    continue
                seg_len = seg["end"] - seg["start"]
                if seg_len == 0:
                    continue
                overlap = max(
                    0, min(seg["end"], gd_end) - max(seg["start"], gd_start),
                )
                ro = min(overlap / seg_len, overlap / entry_len)
                if ro > best_ro:
                    best_ro = ro
                    best_seg = seg

        is_carrier = best_ro >= reciprocal_overlap_threshold
        cn_state = best_seg["cn_state"] if best_seg else ploidy

        if verbose:
            tag = " [CARRIER]" if is_carrier else ""
            _util.vlog(
                f"    Entry {gd_id} ({svtype}, BP1={bp1}, BP2={bp2}  "
                f"ploidy={ploidy}):  "
                f"best_RO={best_ro:.2%}  "
                f"(threshold={reciprocal_overlap_threshold:.0%}){tag}"
            )

        n_bins = sum(
            len(interval_bin_arrays.get(iv, []))
            for iv in covered_intervals
        )

        calls.append({
            "GD_ID":                   gd_id,
            "cluster":                 locus.cluster,
            "chrom":                   locus.chrom,
            "start":                   gd_start,
            "end":                     gd_end,
            "svtype":                  svtype,
            "BP1":                     bp1,
            "BP2":                     bp2,
            "is_terminal":             locus.is_terminal,
            "log_prob_score":          confidence,
            "is_carrier":              is_carrier,
            "reciprocal_overlap":      best_ro,
            "intervals":               covered_intervals,
            "n_bins":                  n_bins,
            "cn_state":                cn_state,
        })

    # ----------------------------------------------------------------
    # Build segment records for downstream use (plotting).
    # For diploid, emit per-haplotype records tagged with haplotype=1/2
    # as well as a total-CN record (haplotype=0).
    # For haploid, emit haplotype=0 only.
    # Each entry is (genomic_start, genomic_end, cn_state, category, haplotype).
    # ----------------------------------------------------------------
    segment_records: List[Tuple[int, int, int, str, int]] = []
    for seg in total_cat_segs:
        segment_records.append(
            (seg["start"], seg["end"], seg["cn_state"], seg["category"], 0)
        )
    if ploidy == 2:
        for seg in hap1_cat_segs:
            segment_records.append(
                (seg["start"], seg["end"], seg["cn_state"], seg["category"], 1)
            )
        for seg in hap2_cat_segs:
            segment_records.append(
                (seg["start"], seg["end"], seg["cn_state"], seg["category"], 2)
            )

    return calls, segment_records
