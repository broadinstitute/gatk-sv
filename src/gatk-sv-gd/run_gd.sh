#!/usr/bin/env bash
set -euo pipefail

# ── Constants ──────────────────────────────────────────────────────────────
# Required inputs — set these before running.
INPUT_DEPTH="/path/to/normalized_depth.tsv"
HIGH_RESOLUTION_DEPTH="/path/to/normalized_depth.tsv"
GD_TABLE="/path/to/gd_table.tsv"

SEG_DUP_BED=""
CENTROMERE_BED=""
ACROCENTRIC_ARM_BED=""
GAPS_BED=""
GTF="/path/to/genes.gtf"

TRANSITION_MATRIX="/path/to/transition_matrix.tsv"
BREAKPOINT_TRANSITION_MATRIX="/path/to/breakpoint_transition_matrix.tsv"
NON_NAHR_TRANSITION_MATRIX="/path/to/non_nahr_transition_matrix.tsv"

# Optional: set to a truth table TSV to run the eval step; leave empty to skip.
TRUTH_TABLE=""

# ── Usage ──────────────────────────────────────────────────────────────────
usage() {
    echo "Usage: $0 <output_dir>" >&2
    exit 1
}

[[ $# -eq 1 ]] || usage

# ── Directories ────────────────────────────────────────────────────────────
WORK_DIR="$1"

PREPROCESS_DIR="${WORK_DIR}/preprocess"
INFER_DIR="${WORK_DIR}/infer"
CALL_DIR="${WORK_DIR}/call"
PLOT_DIR="${WORK_DIR}/plot"
EVAL_DIR="${WORK_DIR}/eval"

# ── Derived paths ──────────────────────────────────────────────────────────
BIN_MAPPINGS="${PREPROCESS_DIR}/bin_mappings.tsv.gz"
CN_POSTERIORS="${INFER_DIR}/cn_posteriors.tsv.gz"
GD_CALLS="${CALL_DIR}/gd_cnv_calls.tsv.gz"
PLOIDY_TABLE="${PREPROCESS_DIR}/ploidy_estimates.tsv"

# ── Step 1: preprocess ─────────────────────────────────────────────────────
echo "[1/5] preprocess"
gatk-sv-gd preprocess \
    -i "${INPUT_DEPTH}" \
    --high-res-counts "${HIGH_RESOLUTION_DEPTH}" \
    -g "${GD_TABLE}" \
    -o "${PREPROCESS_DIR}" \
    -e "${SEG_DUP_BED}" \
    -e "${CENTROMERE_BED}" \
    -e "${ACROCENTRIC_ARM_BED}" \
    --verbose

# ── Step 2: infer ──────────────────────────────────────────────────────────
echo "[2/5] infer"
gatk-sv-gd infer \
    --preprocessed-dir "${PREPROCESS_DIR}" \
    -o "${INFER_DIR}" \
    --verbose

# ── Step 3: call ───────────────────────────────────────────────────────────
echo "[3/5] call"
gatk-sv-gd call \
    --cn-posteriors "${CN_POSTERIORS}" \
    --bin-mappings "${BIN_MAPPINGS}" \
    -g "${GD_TABLE}" \
    -o "${CALL_DIR}" \
    --ploidy-table "${PLOIDY_TABLE}" \
    --verbose

# ── Step 4: plot ───────────────────────────────────────────────────────────
echo "[4/5] plot"
gatk-sv-gd plot \
    --calls "${GD_CALLS}" \
    --cn-posteriors "${CN_POSTERIORS}" \
    --raw-counts "${INPUT_DEPTH}" \
    --high-res-counts "${HIGH_RESOLUTION_DEPTH}" \
    -g "${GD_TABLE}" \
    -o "${PLOT_DIR}" \
    --gaps-bed "${GAPS_BED}" \
    --gtf "${GTF}" \
    --segdup-bed "${SEG_DUP_BED}" \
    --ploidy-table "${PLOIDY_TABLE}"

# ── Step 5: eval (optional) ────────────────────────────────────────────────
if [[ -n "${TRUTH_TABLE}" ]]; then
    echo "[5/5] eval"
    gatk-sv-gd eval \
        --calls "${GD_CALLS}" \
        --truth-table "${TRUTH_TABLE}" \
        -o "${EVAL_DIR}"
else
    echo "[5/5] eval  (skipped — TRUTH_TABLE not set)"
fi
