#!/usr/bin/env bash
set -euo pipefail

# ── Constants ──────────────────────────────────────────────────────────────
# Required inputs — set these before running.
INPUT_DEPTH="/Users/markw/Work/talkowski/sv-pipe-testing/mw_gd/gd_pyro/input/rebatch_1.10kbp_bins.rd.txt.gz"
HIGH_RESOLUTION_DEPTH="/Users/markw/Work/talkowski/sv-pipe-testing/mw_gd/gd_pyro/input/rebatch_1.RD.txt.gz"
GD_TABLE="/Users/markw/Work/talkowski/sv-pipe-testing/mw_gd/gd_pyro/input/GenomicDisorderRegions_hg38_2025-12-05.with_bp.tsv"

SEG_DUP_BED="/Users/markw/Work/talkowski/sv-pipe-testing/mw_gd/gd_pyro/input/hg38_SD.bed.gz"
CENTROMERE_BED="/Users/markw/Work/talkowski/sv-pipe-testing/mw_gd/gd_pyro/input/hg38_centromeres.bed"
ACROCENTRIC_ARM_BED="/Users/markw/Work/talkowski/sv-pipe-testing/mw_gd/gd_pyro/input/hg38_acrocentric_arms.bed"
GAPS_BED="/Users/markw/Work/talkowski/sv-pipe-testing/mw_gd/gd_pyro/input/hg38_gap.bed"
GTF="/Users/markw/Work/talkowski/sv-pipe-testing/mw_gd/gd_pyro/input/gencode.v47.basic.protein_coding.canonical.gtf"

TRANSITION_MATRIX="/Users/markw/Work/talkowski/sv-pipe-testing/mw_gd/gd_pyro/input/transition_matrix.tsv"
BREAKPOINT_TRANSITION_MATRIX="/Users/markw/Work/talkowski/sv-pipe-testing/mw_gd/gd_pyro/input/breakpoint_transition_matrix.tsv"
NON_NAHR_TRANSITION_MATRIX="/Users/markw/Work/talkowski/sv-pipe-testing/mw_gd/gd_pyro/input/transition_matrix.non_nahr.tsv"

# Optional: set to a truth table TSV to run the eval step; leave empty to skip.
TRUTH_TABLE="/Users/markw/Work/talkowski/sv-pipe-testing/mw_gd/gd_pyro/input/asd_cohort.manual_cnv_revision.GD.with_extract_non_NAHR_calls.bed"

# Optional: additional arguments
PREPROCESS_ARGS=""

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
FILTERED_GD_TABLE="${PREPROCESS_DIR}/gd_table_filtered.tsv"
CN_POSTERIORS="${INFER_DIR}/cn_posteriors.tsv.gz"
GD_CALLS="${CALL_DIR}/gd_cnv_calls.tsv.gz"
PLOIDY_TABLE="${PREPROCESS_DIR}/ploidy_estimates.tsv"

# ── Step 1: preprocess ─────────────────────────────────────────────────────
echo "[1/5] preprocess"
rm -r "${PREPROCESS_DIR}"
gatk-sv-gd preprocess \
    -i "${INPUT_DEPTH}" \
    --high-res-counts "${HIGH_RESOLUTION_DEPTH}" \
    -g "${GD_TABLE}" \
    -o "${PREPROCESS_DIR}" \
    -e "${SEG_DUP_BED}" \
    -e "${CENTROMERE_BED}" \
    -e "${ACROCENTRIC_ARM_BED}" \
    --verbose \
    ${PREPROCESS_ARGS}

# ── Step 2: infer ──────────────────────────────────────────────────────────
echo "[2/5] infer"
rm -r "${INFER_DIR}"
gatk-sv-gd infer \
    --preprocessed-dir "${PREPROCESS_DIR}" \
    -o "${INFER_DIR}" \
    --verbose

# ── Step 3: call ───────────────────────────────────────────────────────────
echo "[3/5] call"
rm -r "${CALL_DIR}"
gatk-sv-gd call \
    --cn-posteriors "${CN_POSTERIORS}" \
    --bin-mappings "${BIN_MAPPINGS}" \
    -g "${FILTERED_GD_TABLE}" \
    -o "${CALL_DIR}" \
    --ploidy-table "${PLOIDY_TABLE}" \
    --transition-matrix "${TRANSITION_MATRIX}" \
    --breakpoint-transition-matrix "${BREAKPOINT_TRANSITION_MATRIX}" \
    --non-nahr-transition-matrix "${NON_NAHR_TRANSITION_MATRIX}" \
    --verbose

# ── Step 4: plot ───────────────────────────────────────────────────────────
echo "[4/5] plot"
rm -r "${PLOT_DIR}"
gatk-sv-gd plot \
    --calls "${GD_CALLS}" \
    --cn-posteriors "${CN_POSTERIORS}" \
    --raw-counts "${INPUT_DEPTH}" \
    --high-res-counts "${HIGH_RESOLUTION_DEPTH}" \
    -g "${FILTERED_GD_TABLE}" \
    -o "${PLOT_DIR}" \
    --gaps-bed "${GAPS_BED}" \
    --gtf "${GTF}" \
    --segdup-bed "${SEG_DUP_BED}" \
    --ploidy-table "${PLOIDY_TABLE}"

# ── Step 5: eval (optional) ────────────────────────────────────────────────
if [[ -n "${TRUTH_TABLE}" ]]; then
    echo "[5/5] eval"
    rm -r "${EVAL_DIR}"
    gatk-sv-gd eval \
        --calls "${GD_CALLS}" \
        --truth-table "${TRUTH_TABLE}" \
        --ploidy-table "${PLOIDY_TABLE}" \
        -o "${EVAL_DIR}"
else
    echo "[5/5] eval  (skipped — TRUTH_TABLE not set)"
fi
