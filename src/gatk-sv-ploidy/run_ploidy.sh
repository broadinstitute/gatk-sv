#!/usr/bin/env bash
# run_ploidy.sh — run the full gatk-sv-ploidy pipeline
#
# Usage:
#   ./run_ploidy.sh <output_dir>
#
# Required constants (edit before running):
#   INPUT_DEPTH   — raw bins × samples depth TSV (may be gzipped)
#   TRUTH_JSON    — JSON mapping sample ID → true aneuploidy type
#                   leave empty ("") to skip the eval step

set -euo pipefail

# ── user-configurable constants ──────────────────────────────────────────────
INPUT_DEPTH="/path/to/depth.tsv"
TRUTH_JSON=""          # set to "/path/to/truth.json" to enable eval
# ─────────────────────────────────────────────────────────────────────────────

usage() {
    echo "Usage: $0 <output_dir>" >&2
    echo "  output_dir  Root directory for all pipeline outputs" >&2
    exit 1
}

[[ $# -eq 1 ]] || usage
WORK_DIR="$1"

# ── per-step output directories ───────────────────────────────────────────────
PREPROCESS_DIR="${WORK_DIR}/preprocess"
INFER_DIR="${WORK_DIR}/infer"
CALL_DIR="${WORK_DIR}/call"
PLOT_DIR="${WORK_DIR}/plot"
EVAL_DIR="${WORK_DIR}/eval"

# ── derived paths (outputs of earlier steps used as inputs to later steps) ───
PREPROCESSED_DEPTH="${PREPROCESS_DIR}/preprocessed_depth.tsv"
CHROM_STATS="${INFER_DIR}/chromosome_stats.tsv"
BIN_STATS="${INFER_DIR}/bin_stats.tsv.gz"
TRAINING_LOSS="${INFER_DIR}/training_loss.tsv"
PREDICTIONS="${CALL_DIR}/aneuploidy_type_predictions.tsv"

echo "=== gatk-sv-ploidy pipeline ==="
echo "  Input depth : ${INPUT_DEPTH}"
echo "  Work dir    : ${WORK_DIR}"
echo ""

# ── step 1: preprocess ───────────────────────────────────────────────────────
echo "[1/5] preprocess"
gatk-sv-ploidy preprocess \
    -i "${INPUT_DEPTH}" \
    -o "${PREPROCESS_DIR}"

# ── step 2: infer ────────────────────────────────────────────────────────────
echo "[2/5] infer"
gatk-sv-ploidy infer \
    -i "${PREPROCESSED_DEPTH}" \
    -o "${INFER_DIR}"

# ── step 3: call ─────────────────────────────────────────────────────────────
echo "[3/5] call"
gatk-sv-ploidy call \
    -c "${CHROM_STATS}" \
    -o "${CALL_DIR}"

# ── step 4: plot ─────────────────────────────────────────────────────────────
echo "[4/5] plot"
gatk-sv-ploidy plot \
    -c "${CHROM_STATS}" \
    -b "${BIN_STATS}" \
    -t "${TRAINING_LOSS}" \
    -s "${PREDICTIONS}" \
    -o "${PLOT_DIR}"

# ── step 5: eval (optional — skipped if TRUTH_JSON is unset) ─────────────────
if [[ -n "${TRUTH_JSON}" ]]; then
    echo "[5/5] eval"
    gatk-sv-ploidy eval \
        -p "${PREDICTIONS}" \
        -t "${TRUTH_JSON}" \
        -o "${EVAL_DIR}"
else
    echo "[5/5] eval  (skipped — set TRUTH_JSON to enable)"
fi

echo ""
echo "=== pipeline complete ==="
echo "  Outputs in: ${WORK_DIR}"
