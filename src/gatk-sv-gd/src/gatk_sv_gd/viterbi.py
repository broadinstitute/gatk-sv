"""
Viterbi HMM Segmentation for GD CNV Calling

Implements the Viterbi algorithm on per-bin CN posterior probabilities to
produce smooth copy-number segmentations, then checks each GD entry's
breakpoint pattern against the resulting path to call carriers.

The hidden state is a *pair* of per-haplotype copy numbers ``(h1, h2)`` with
``h1 <= h2`` (canonical ordering).  The observed per-bin posterior over these
pair states acts directly as the Viterbi observation prior.  Transitions
factorise as the product of independent per-haplotype transitions.  Any bias
toward or against specific haplotype states should be encoded in the
transition matrix itself.

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
        isinstance(first_col, str) and (first_col == "" or first_col.lower() in ("state", "cn", "index"))
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


def run_list_viterbi(
    state_log_priors: np.ndarray,
    transition_matrix: np.ndarray,
    initial_log_probs: Optional[np.ndarray] = None,
    breakpoint_mask: Optional[np.ndarray] = None,
    breakpoint_transition_matrix: Optional[np.ndarray] = None,
    n_best: int = 10,
) -> List[Tuple[np.ndarray, float, float]]:
    """Parallel List Viterbi Algorithm — return the *n_best* most-likely paths.

    This is the PLVA of Seshadri & Sundberg (1994).  At every time step and
    for every state, the algorithm maintains the top *n_best* partial paths
    (ranked by cumulative log-probability).  Back-tracing then yields up to
    *n_best* globally distinct complete paths.

    Interface mirrors :func:`run_viterbi` except that it returns a list of
    ``(path, confidence, mean_path_score)`` tuples sorted from best to worst.

    Args:
        state_log_priors: (T, K) log posteriors per bin.
        transition_matrix: (K, K) transition probability matrix.
        initial_log_probs: (K,) log probs for the initial state.
        breakpoint_mask: Optional (T-1,) boolean mask for breakpoint steps.
        breakpoint_transition_matrix: Optional (K, K) breakpoint matrix.
        n_best: Number of best paths to return.

    Returns:
        List of ``(path, confidence, mean_path_score)`` tuples, length
        <= *n_best*, sorted by descending full Viterbi path score (best
        first). ``confidence`` is the mean observation log-prior along
        the path, while ``mean_path_score`` is the mean total Viterbi
        log-score per bin (observations + transitions + initial state).
    """
    T, K = state_log_priors.shape
    log_trans = np.log(transition_matrix)

    log_bp_trans: Optional[np.ndarray] = None
    if breakpoint_mask is not None and breakpoint_transition_matrix is not None:
        log_bp_trans = np.log(breakpoint_transition_matrix)

    if initial_log_probs is None:
        initial_log_probs = np.full(K, -np.log(K))

    # V[t, k, r] = log-probability of the r-th best partial path ending in
    # state k at time t.
    NEG_INF = -1e30
    V = np.full((T, K, n_best), NEG_INF)
    # Back-pointers: (previous state, rank in that state's list)
    ptr_state = np.zeros((T, K, n_best), dtype=np.int32)
    ptr_rank = np.zeros((T, K, n_best), dtype=np.int32)

    # t = 0: only rank 0 is populated
    V[0, :, 0] = initial_log_probs + state_log_priors[0]

    for t in range(1, T):
        lt = (
            log_bp_trans
            if (log_bp_trans is not None and breakpoint_mask[t - 1])
            else log_trans
        )
        for k in range(K):
            # Collect all (score, prev_state, prev_rank) candidates
            # that could lead to state k at time t.
            # There are up to K * n_best candidates.
            candidates = []
            for j in range(K):
                trans_score = lt[j, k]
                for r in range(n_best):
                    score = V[t - 1, j, r]
                    if score <= NEG_INF + 1:
                        break  # remaining ranks for state j are -inf
                    candidates.append((score + trans_score, j, r))

            # Keep top n_best by score (partial sort)
            if len(candidates) <= n_best:
                candidates.sort(key=lambda x: x[0], reverse=True)
                top = candidates
            else:
                # Use argpartition for efficiency on large lists
                top = sorted(candidates, key=lambda x: x[0], reverse=True)[:n_best]

            for rank, (sc, prev_j, prev_r) in enumerate(top):
                V[t, k, rank] = sc + state_log_priors[t, k]
                ptr_state[t, k, rank] = prev_j
                ptr_rank[t, k, rank] = prev_r

    # ---- Extract globally best n_best paths ----
    # Collect all (score, state, rank) at the final time step, pick top n_best.
    final_candidates = []
    for k in range(K):
        for r in range(n_best):
            if V[T - 1, k, r] > NEG_INF + 1:
                final_candidates.append((V[T - 1, k, r], k, r))
    final_candidates.sort(key=lambda x: x[0], reverse=True)

    results: List[Tuple[np.ndarray, float, float]] = []
    seen_paths: set = set()
    for final_score, end_k, end_r in final_candidates[:n_best * 2]:
        # Back-trace
        path = np.empty(T, dtype=int)
        path[T - 1] = end_k
        cur_k, cur_r = end_k, end_r
        for t in range(T - 2, -1, -1):
            prev_k = int(ptr_state[t + 1, cur_k, cur_r])
            prev_r = int(ptr_rank[t + 1, cur_k, cur_r])
            path[t] = prev_k
            cur_k, cur_r = prev_k, prev_r

        # De-duplicate identical paths
        path_key = path.tobytes()
        if path_key in seen_paths:
            continue
        seen_paths.add(path_key)

        path_log_priors = state_log_priors[np.arange(T), path]
        confidence = float(np.mean(path_log_priors))
        mean_path_score = float(final_score / T)
        results.append((path, confidence, mean_path_score))
        if len(results) >= n_best:
            break

    return results


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

    **Combinatorial correction:** when the source pair is *homozygous*
    (a1 == a2) and the destination is *heterozygous* (lo != hi), the inner
    loop visits the same canonical destination twice (once where haplotype 1
    mutates and once where haplotype 2 mutates).  This doubles the raw
    probability and gives homozygous→heterozygous transitions an unfair 2×
    advantage — e.g. (1,1)→(0,1) appears twice as likely as
    (0,1)→(0,2).  We neutralise this by halving the contribution when
    leaving a homozygous state for a heterozygous one.

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

                # Neutralise the combinatorial bonus: when we start in a
                # homozygous state (a1 == a2) and land in a heterozygous
                # state (lo != hi), the inner loop visits this outcome
                # twice (once for Hap1 mutating, once for Hap2).  Halving
                # removes the 2× advantage that would otherwise make the
                # model prefer bouncing back to (1,1) to start new
                # variants rather than chaining them (e.g. 0,1 → 0,2).
                if a1 == a2 and lo != hi:
                    prob *= 0.5

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
    n_best: int = 20,
    score_tolerance: float = 1e-6,
    verbose: bool = False,
    sample_id: str = "",
    cluster: str = "",
) -> Tuple[List[Tuple[np.ndarray, np.ndarray, float]], List[Tuple[int, int]]]:
    """Run the Viterbi algorithm over diploid pair states.

    Uses the Parallel List Viterbi Algorithm (PLVA) to find the top
    *n_best* highest-scoring paths.  Returns ALL paths whose per-bin
    score is within *score_tolerance* of the best, decomposed into
    per-haplotype CN sequences.  The caller can then try each candidate
    against GD entries and take the best match.

    Args:
        total_cn_log_priors: (T, K_total) log posteriors from the Pyro model.
        hap_transition_matrix: (K_hap, K_hap) per-haplotype transition matrix.
        breakpoint_mask: optional (T-1,) boolean mask.
        hap_breakpoint_transition_matrix: optional per-haplotype breakpoint
            transition matrix.
        n_best: Number of candidate paths to consider (default 20).
        score_tolerance: Maximum per-bin log-probability difference from the
            best path for a candidate to be included (default 0.001).
        verbose: Emit diagnostic messages.
        sample_id: For logging.
        cluster: For logging.

    Returns:
        Tuple of:
          - candidates: list of (hap1_path, hap2_path, confidence) tuples,
            ordered from best to worst score.
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

    # Run list Viterbi to get candidate paths
    candidates = run_list_viterbi(
        obs, pair_trans, initial,
        breakpoint_mask=breakpoint_mask,
        breakpoint_transition_matrix=pair_bp_trans,
        n_best=n_best,
    )

    if not candidates:
        # Fallback: return all-reference
        T = total_cn_log_priors.shape[0]
        fallback = [(np.ones(T, dtype=int), np.ones(T, dtype=int), -1e30)]
        return fallback, pair_states

    # Collect all near-tied candidate decompositions.
    # Each candidate within score_tolerance of the best is decomposed
    # into per-haplotype paths and returned.  The caller tries all of
    # them against GD entries and takes the best match.
    #
    # GUARD: only return multiple candidates when the best path already
    # contains variant pair states (at least one bin is NOT (1,1)).
    # If the best path is entirely reference, alternative decompositions
    # (e.g. (0,2) instead of (1,1)) are not meaningful variants — they
    # just split reference CN=2 differently between haplotypes.  Including
    # them would create spurious DUP/DEL segments that match GD entries
    # and cause false positive calls.
    best_score = candidates[0][2]

    # Determine the (1,1) pair-state index for the guard check.
    ref_pair_idx = None
    for pi, (h1, h2) in enumerate(pair_states):
        if h1 == 1 and h2 == 1:
            ref_pair_idx = pi
            break

    best_has_variants = (
        ref_pair_idx is not None and not np.all(candidates[0][0] == ref_pair_idx)
    )

    result_candidates: List[Tuple[np.ndarray, np.ndarray, float]] = []

    limit = len(candidates) if best_has_variants else 1
    for cand_path, cand_conf, cand_score in candidates[:limit]:
        score_gap = best_score - cand_score  # always >= 0
        if score_gap > score_tolerance:
            break  # remaining candidates are even worse
        # Decompose pair-state path into per-haplotype CN paths
        T = len(cand_path)
        hap1 = np.empty(T, dtype=int)
        hap2 = np.empty(T, dtype=int)
        for t in range(T):
            h1, h2 = pair_states[cand_path[t]]
            hap1[t] = h1
            hap2[t] = h2
        result_candidates.append((hap1, hap2, cand_conf))

    if verbose:
        _util.vlog(f"    [LIST-VIT] {sample_id}  {cluster}: "
                   f"{len(result_candidates)} near-tied candidate(s) "
                   f"from {len(candidates)} total")

    return result_candidates, pair_states


def _run_pair_state_viterbi(
    pair_log_priors: np.ndarray,
    pair_states: List[Tuple[int, int]],
    hap_transition_matrix: np.ndarray,
    breakpoint_mask: Optional[np.ndarray],
    hap_breakpoint_transition_matrix: Optional[np.ndarray],
    n_best: int = 20,
    score_tolerance: float = 1e-6,
    verbose: bool = False,
    sample_id: str = "",
    cluster: str = "",
) -> List[Tuple[np.ndarray, np.ndarray, float]]:
    """Run Viterbi directly on emitted pair-state posteriors."""
    pair_trans = _build_diploid_transition_matrix(
        hap_transition_matrix, pair_states,
    )
    pair_bp_trans = None
    if hap_breakpoint_transition_matrix is not None:
        pair_bp_trans = _build_diploid_transition_matrix(
            hap_breakpoint_transition_matrix, pair_states,
        )

    initial = _build_diploid_initial_log_probs(pair_states)
    candidates = run_list_viterbi(
        pair_log_priors,
        pair_trans,
        initial,
        breakpoint_mask=breakpoint_mask,
        breakpoint_transition_matrix=pair_bp_trans,
        n_best=n_best,
    )

    if not candidates:
        n_bins = pair_log_priors.shape[0]
        return [(np.ones(n_bins, dtype=int), np.ones(n_bins, dtype=int), -1e30)]

    ref_pair_idx = None
    for pair_idx, (h1, h2) in enumerate(pair_states):
        if h1 == 1 and h2 == 1:
            ref_pair_idx = pair_idx
            break

    best_score = candidates[0][2]
    best_has_variants = (
        ref_pair_idx is not None and not np.all(candidates[0][0] == ref_pair_idx)
    )

    result_candidates: List[Tuple[np.ndarray, np.ndarray, float]] = []
    limit = len(candidates) if best_has_variants else 1
    for cand_path, cand_conf, cand_score in candidates[:limit]:
        score_gap = best_score - cand_score
        if score_gap > score_tolerance:
            break

        n_bins = len(cand_path)
        hap1 = np.empty(n_bins, dtype=int)
        hap2 = np.empty(n_bins, dtype=int)
        for bin_idx, pair_idx in enumerate(cand_path):
            h1, h2 = pair_states[int(pair_idx)]
            hap1[bin_idx] = h1
            hap2[bin_idx] = h2
        result_candidates.append((hap1, hap2, cand_conf))

    if verbose:
        _util.vlog(f"    [PAIR-VIT] {sample_id}  {cluster}: "
                   f"{len(result_candidates)} near-tied candidate(s) "
                   f"from {len(candidates)} total")

    return result_candidates


def _haplotype_balance_score(
    pair_path: np.ndarray,
    pair_states: List[Tuple[int, int]],
) -> int:
    """Score a pair-state path by haplotype imbalance.

    Decomposes the path into per-haplotype CN sequences, counts the number
    of contiguous non-reference (CN != 1) segments on each haplotype, and
    returns ``max(n_variant_segs_h1, n_variant_segs_h2)``.

    Lower is more balanced.  A path with 1 variant on each haplotype
    scores 1, while 2 variants on one haplotype and 0 on the other
    scores 2.

    Returns:
        Maximum number of variant segments on either haplotype.
    """
    if len(pair_path) == 0:
        return 0

    def _count_variant_segments(hap_path: np.ndarray, ref_cn: int = 1) -> int:
        n_segs = 0
        in_variant = False
        for cn in hap_path:
            if cn != ref_cn:
                if not in_variant:
                    n_segs += 1
                    in_variant = True
            else:
                in_variant = False
        return n_segs

    T = len(pair_path)
    hap1 = np.empty(T, dtype=int)
    hap2 = np.empty(T, dtype=int)
    for t in range(T):
        h1, h2 = pair_states[pair_path[t]]
        hap1[t] = h1
        hap2[t] = h2

    return max(_count_variant_segments(hap1), _count_variant_segments(hap2))


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


# Maximum gap (bp) to bridge when merging nearby non-REF segments of the
# same svtype for GD-entry matching.  The diploid pair-state Viterbi can
# briefly bounce back to reference mid-event (e.g. (1,2) → (1,1) → (1,2))
# when the data in a small region looks reference-like, splitting a single
# large event into two smaller segments.  Bridging short gaps during the
# matching step allows these fragments to be reconnected without altering
# the Viterbi path itself.
_MAX_MATCH_GAP_BP: int = 500_000


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
      - Merge overlapping, adjacent, or nearby (within
        :data:`_MAX_MATCH_GAP_BP`) segments of the same svtype.
      - Return the merged list.  ``cn_state`` is set to the per-haplotype
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
        # Merge overlapping/adjacent/nearby
        merged_start = segs[0]["start"]
        merged_end = segs[0]["end"]
        merged_cn = segs[0]["cn_state"]
        merged_n = segs[0]["n_bins"]
        for s in segs[1:]:
            if s["start"] <= merged_end + _MAX_MATCH_GAP_BP:
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


def _build_haplotype_match_segments(
    hap_segments: List[dict],
    haplotype: int,
) -> List[dict]:
    """Merge nearby same-svtype segments on a single haplotype."""
    merged_segments: List[dict] = []
    for svtype in ("DEL", "DUP"):
        segs = sorted(
            [segment for segment in hap_segments if segment["category"] == svtype],
            key=lambda segment: segment["start"],
        )
        if not segs:
            continue

        merged = dict(segs[0])
        for segment in segs[1:]:
            if segment["start"] <= merged["end"] + _MAX_MATCH_GAP_BP:
                merged["end"] = max(merged["end"], segment["end"])
                if segment["n_bins"] > merged["n_bins"]:
                    merged["cn_state"] = segment["cn_state"]
                merged["n_bins"] += segment["n_bins"]
            else:
                merged["haplotype"] = haplotype
                merged_segments.append(merged)
                merged = dict(segment)
        merged["haplotype"] = haplotype
        merged_segments.append(merged)

    return merged_segments


# =============================================================================
# Viterbi-based GD CNV calling
# =============================================================================


def viterbi_call_gd_cnv(
    locus: GDLocus,
    sample_pair_probs: np.ndarray,
    pair_states: List[Tuple[int, int]],
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
        sample_pair_probs: (n_bins, n_pair_states) pair-state posterior
            probabilities for one sample across all bins at this locus.
        pair_states: Canonical pair-state labels aligned to
            ``sample_pair_probs`` columns.
        transition_matrix: Per-haplotype transition matrix. High diagonal =
            sticky segmentation.
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

    n_pair_states = sample_pair_probs.shape[1]
    pair_priors = sample_pair_probs[all_bins]
    pair_priors = np.maximum(pair_priors, 1e-30)
    pair_log_priors = np.log(pair_priors)

    max_hap_state = max(max(h1, h2) for h1, h2 in pair_states) + 1
    tm = transition_matrix[:max_hap_state, :max_hap_state]

    # Build a per-step breakpoint mask: True where a transition crosses the
    # START or END coordinate of a known breakpoint SD-block range.  At those
    # positions the breakpoint_transition_matrix (less sticky) will be used
    # instead of the regular transition_matrix, encoding a higher prior on
    # CN-state changes at known recurrent breakpoints.
    breakpoint_mask: Optional[np.ndarray] = None
    if breakpoint_transition_matrix is not None and len(all_bins) > 1:
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
        breakpoint_transition_matrix[:max_hap_state, :max_hap_state]
        if breakpoint_transition_matrix is not None else None
    )

    pair_candidates = _run_pair_state_viterbi(
        pair_log_priors,
        pair_states[:n_pair_states],
        tm,
        breakpoint_mask,
        bp_tm,
        verbose=verbose,
        sample_id=sample_id,
        cluster=locus.cluster,
    )

    hap1_path, hap2_path, confidence = pair_candidates[0]
    total_path = hap1_path + hap2_path

    if verbose:
        pair_state_path = np.argmax(pair_priors, axis=1)
        pair_segs = _extract_segments(pair_state_path)
        seg_str = " -> ".join(
            f"PAIR{pair_states[state]}x{end - start}" for start, end, state in pair_segs
        )
        _util.vlog(f"  [PAIR-VIT] {sample_id}  {locus.cluster}  "
                   f"{len(all_bins)} bins  sample_ploidy={ploidy}  "
                   f"confidence={confidence:+.4f}")
        _util.vlog(f"    Pair posterior peaks: {seg_str}")
        h1_segs = _extract_segments(hap1_path)
        h2_segs = _extract_segments(hap2_path)
        _util.vlog(f"    Hap1: {' -> '.join(f'CN{s}x{e-b}' for b,e,s in h1_segs)}")
        _util.vlog(f"    Hap2: {' -> '.join(f'CN{s}x{e-b}' for b,e,s in h2_segs)}")

    hap1_cat_segs = _build_category_segments(hap1_path, all_bins, bin_coords, ploidy=1)
    hap2_cat_segs = _build_category_segments(hap2_path, all_bins, bin_coords, ploidy=1)
    total_cat_segs = _build_category_segments(total_path, all_bins, bin_coords, ploidy=2)

    all_candidate_match_segs: List[List[dict]] = []
    for cand_h1, cand_h2, _ in pair_candidates:
        cand_h1_segs = _build_category_segments(cand_h1, all_bins, bin_coords, ploidy=1)
        cand_h2_segs = _build_category_segments(cand_h2, all_bins, bin_coords, ploidy=1)
        cand_match = _build_haplotype_match_segments(
            cand_h1_segs, haplotype=1,
        ) + _build_haplotype_match_segments(cand_h2_segs, haplotype=2)
        all_candidate_match_segs.append(cand_match)

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
    # For diploid, try ALL near-tied candidate decompositions and take the
    # best reciprocal overlap across all of them.  This handles cases where
    # the best-scoring path has a reference bounce that splits an event,
    # but an equally-likely alternative path has a contiguous event that
    # matches the GD entry better.
    # For haploid, use a single set of match segments.
    # ----------------------------------------------------------------

    match_segments_list = all_candidate_match_segs

    entry_results: List[dict] = []
    for entry in locus.gd_entries:
        gd_id = entry["GD_ID"]
        svtype = entry["svtype"]
        gd_start = entry["start_GRCh38"]
        gd_end = entry["end_GRCh38"]
        bp1 = entry["BP1"]
        bp2 = entry["BP2"]

        covered_tuples = locus.get_intervals_between(bp1, bp2)
        covered_intervals = [name for _, _, name in covered_tuples]

        entry_len = gd_end - gd_start
        best_ro = 0.0
        best_seg = None
        if entry_len > 0:
            # First pass: check the best path's match segments.
            best_path_segs = match_segments_list[0] if match_segments_list else []
            for seg in best_path_segs:
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

            # Second pass: for DEL entries, also inspect alternative
            # near-tied pair-state decompositions even when the top path
            # has no overlap signal. In practice, the best-scoring path can
            # stay reference-like through a shallow deletion while an
            # equally likely alternative carries a contiguous haploid DEL.
            # Keep the stricter gating for DUP entries to avoid inflating
            # spurious duplication matches from phantom decompositions.
            should_try_alternatives = (
                len(match_segments_list) > 1 and (best_ro > 0 or svtype == "DEL")
            )
            if should_try_alternatives:
                for match_segments in match_segments_list[1:]:
                    for seg in match_segments:
                        if seg["category"] != svtype:
                            continue
                        seg_len = seg["end"] - seg["start"]
                        if seg_len == 0:
                            continue
                        overlap = max(
                            0,
                            min(seg["end"], gd_end) - max(seg["start"], gd_start),
                        )
                        ro = min(overlap / seg_len, overlap / entry_len)
                        if ro > best_ro:
                            best_ro = ro
                            best_seg = seg

        n_bins = sum(
            len(interval_bin_arrays.get(iv, []))
            for iv in covered_intervals
        )

        overlap_bp = 0
        if best_seg is not None:
            overlap_bp = max(
                0,
                min(best_seg["end"], gd_end) - max(best_seg["start"], gd_start),
            )

        entry_results.append({
            "GD_ID": gd_id,
            "cluster": locus.cluster,
            "chrom": locus.chrom,
            "start": gd_start,
            "end": gd_end,
            "svtype": svtype,
            "BP1": bp1,
            "BP2": bp2,
            "is_terminal": locus.is_terminal,
            "log_prob_score": confidence,
            "is_carrier": False,
            "reciprocal_overlap": best_ro,
            "intervals": covered_intervals,
            "n_bins": n_bins,
            "sample_ploidy": ploidy,
            "haplotype": best_seg["haplotype"] if best_seg else np.nan,
            "hap_cn_state": best_seg["cn_state"] if best_seg else np.nan,
            "matched_seg_start": best_seg["start"] if best_seg else np.nan,
            "matched_seg_end": best_seg["end"] if best_seg else np.nan,
            "matched_seg_n_bins": best_seg["n_bins"] if best_seg else 0,
            "entry_len": gd_end - gd_start,
            "overlap_bp": overlap_bp,
        })

    # Keep only the best carrier per svtype for this locus.  This avoids
    # emitting multiple nested GD calls for the same underlying event when
    # alternative near-tied diploid decompositions can explain several
    # overlapping entries.  Rank by reciprocal overlap first, then overlap
    # length, then GD span length so that larger events win close ties.
    best_carrier_idx_by_svtype: Dict[str, int] = {}
    for idx, call in enumerate(entry_results):
        if call["reciprocal_overlap"] < reciprocal_overlap_threshold:
            continue
        svtype = call["svtype"]
        prev_idx = best_carrier_idx_by_svtype.get(svtype)
        if prev_idx is None:
            best_carrier_idx_by_svtype[svtype] = idx
            continue
        prev = entry_results[prev_idx]
        key = (call["reciprocal_overlap"], call["overlap_bp"], call["entry_len"])
        prev_key = (
            prev["reciprocal_overlap"],
            prev["overlap_bp"],
            prev["entry_len"],
        )
        if key > prev_key:
            best_carrier_idx_by_svtype[svtype] = idx

    calls: List[dict] = []
    for idx, call in enumerate(entry_results):
        call["is_carrier"] = best_carrier_idx_by_svtype.get(call["svtype"]) == idx

        if verbose:
            tag = " [CARRIER]" if call["is_carrier"] else ""
            _util.vlog(
                f"    Entry {call['GD_ID']} ({call['svtype']}, BP1={call['BP1']}, "
                f"BP2={call['BP2']}  ploidy={ploidy}):  "
                f"best_RO={call['reciprocal_overlap']:.2%}  "
                f"(threshold={reciprocal_overlap_threshold:.0%}){tag}"
            )

        del call["entry_len"]
        del call["overlap_bp"]
        calls.append(call)

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
