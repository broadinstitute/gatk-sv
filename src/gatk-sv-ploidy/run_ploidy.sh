#!/usr/bin/env bash
# run_ploidy.sh — run the full gatk-sv-ploidy pipeline
#
# Usage:
#   ./run_ploidy.sh <output_dir>
#
# Required constants (configured via CLI args):
#   --input-depth  : raw bins × samples depth TSV (may be gzipped) (required)
#   --truth-json   : JSON mapping sample ID → true aneuploidy type (optional)
#   --work-dir     : Root directory for all pipeline outputs (required)

set -euo pipefail

INPUT_DEPTH=""
TRUTH_JSON=""
WORK_DIR=""
DRY_RUN="false"
ENABLE_PPD="false"
PREPROCESS_ARGS=""
MODEL_ARGS=""
PLOT_ARGS=""
HIGHLIGHT_SAMPLE=""
SITE_DEPTH_LIST=""
POOR_REGIONS=""
MIN_POOR_REGION_COVERAGE="0.5"

usage() {
    echo "Usage: $0 --input-depth PATH --work-dir DIR [--truth-json PATH] [--site-depth-list PATH] [--poor-regions PATH] [--ppd]" >&2
    echo "  --input-depth PATH       Raw bins×samples depth TSV (may be gzipped)" >&2
    echo "  --work-dir DIR           Root directory for all pipeline outputs" >&2
    echo "  --truth-json PATH        Optional truth JSON to enable evaluation" >&2
    echo "  --site-depth-list PATH   Optional file listing per-sample SD file paths (one per line)" >&2
    echo "  --poor-regions PATH      Optional BED of poor regions to remove during preprocess" >&2
    echo "  --min-poor-region-coverage FLOAT  Min fraction of a bin overlapped by poor regions to filter it" >&2
    echo "  --ppd                    Run posterior predictive checks and enable PPD plots" >&2
    exit 1
}

if [[ $# -eq 0 ]]; then
    usage
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-depth)
            INPUT_DEPTH="$2"; shift 2;;
        --input-depth=*)
            INPUT_DEPTH="${1#*=}"; shift;;
        --truth-json)
            TRUTH_JSON="$2"; shift 2;;
        --truth-json=*)
            TRUTH_JSON="${1#*=}"; shift;;
        --preprocess-args)
            PREPROCESS_ARGS="$2"; shift 2;;
        --preprocess-args=*)
            PREPROCESS_ARGS="${1#*=}"; shift;;
        --model-args)
            MODEL_ARGS="$2"; shift 2;;
        --model-args=*)
            MODEL_ARGS="${1#*=}"; shift;;
        --plot-args)
            PLOT_ARGS="$2"; shift 2;;
        --plot-args=*)
            PLOT_ARGS="${1#*=}"; shift;;
        --highlight-sample)
            HIGHLIGHT_SAMPLE="$2"; shift 2;;
        --highlight-sample=*)
            HIGHLIGHT_SAMPLE="${1#*=}"; shift;;
        --site-depth-list)
            SITE_DEPTH_LIST="$2"; shift 2;;
        --site-depth-list=*)
            SITE_DEPTH_LIST="${1#*=}"; shift;;
        --poor-regions)
            POOR_REGIONS="$2"; shift 2;;
        --poor-regions=*)
            POOR_REGIONS="${1#*=}"; shift;;
        --min-poor-region-coverage)
            MIN_POOR_REGION_COVERAGE="$2"; shift 2;;
        --min-poor-region-coverage=*)
            MIN_POOR_REGION_COVERAGE="${1#*=}"; shift;;
        --dry-run)
            DRY_RUN="true"; shift;;
        --dry-run=*)
            DRY_RUN="true"; shift;;
        --ppd)
            ENABLE_PPD="true"; shift;;
        --ppd=*)
            ENABLE_PPD="true"; shift;;
        --work-dir)
            WORK_DIR="$2"; shift 2;;
        --work-dir=*)
            WORK_DIR="${1#*=}"; shift;;
        -h|--help)
            usage;;
        *)
            echo "Unknown argument: $1" >&2; usage;;
    esac
done

if [[ -z "${INPUT_DEPTH}" ]]; then
    echo "Error: --input-depth is required" >&2
    usage
fi

if [[ -z "${WORK_DIR}" ]]; then
    echo "Error: --work-dir is required" >&2
    usage
fi

if [[ "${DRY_RUN}" == "true" ]]; then
    echo "DRY-RUN mode: commands will be printed but not executed"
fi

TOTAL_STEPS=5
if [[ "${ENABLE_PPD}" == "true" ]]; then
    TOTAL_STEPS=$((TOTAL_STEPS + 1))
fi

# ── per-step output directories ───────────────────────────────────────────────
PREPROCESS_DIR="${WORK_DIR}/preprocess"
INFER_DIR="${WORK_DIR}/infer"
CALL_DIR="${WORK_DIR}/call"
PLOT_DIR="${WORK_DIR}/plot"
PPD_DIR="${WORK_DIR}/ppd"
EVAL_DIR="${WORK_DIR}/eval"

# ── derived paths (outputs of earlier steps used as inputs to later steps) ───
PREPROCESSED_DEPTH="${PREPROCESS_DIR}/preprocessed_depth.tsv"
SITE_DATA="${PREPROCESS_DIR}/site_data.npz"
CHROM_STATS="${INFER_DIR}/chromosome_stats.tsv"
BIN_STATS="${INFER_DIR}/bin_stats.tsv.gz"
TRAINING_LOSS="${INFER_DIR}/training_loss.tsv"
INFERENCE_ARTIFACTS="${INFER_DIR}/inference_artifacts.npz"
PREDICTIONS="${CALL_DIR}/aneuploidy_type_predictions.tsv"
PPD_BIN_SUMMARY="${PPD_DIR}/ppd_bin_summary.tsv.gz"
PPD_CHR_SUMMARY="${PPD_DIR}/ppd_chromosome_summary.tsv"

echo "=== gatk-sv-ploidy pipeline ==="
echo "  Input depth      : ${INPUT_DEPTH}"
echo "  Site depth list  : ${SITE_DEPTH_LIST:-<none>}"
echo "  Poor regions     : ${POOR_REGIONS:-<none>}"
echo "  Poor region cov  : ${MIN_POOR_REGION_COVERAGE}"
echo "  Run PPD          : ${ENABLE_PPD}"
echo "  Work dir         : ${WORK_DIR}"
echo ""

# ── step 1: preprocess ───────────────────────────────────────────────────────
echo "[1/5] preprocess"
SD_ARGS=""
if [[ -n "${SITE_DEPTH_LIST}" ]]; then
    SD_ARGS="--site-depth-list ${SITE_DEPTH_LIST}"
fi
PR_ARGS=""
if [[ -n "${POOR_REGIONS}" ]]; then
    PR_ARGS="--poor-regions ${POOR_REGIONS} --min-poor-region-coverage ${MIN_POOR_REGION_COVERAGE}"
fi
if [[ "${DRY_RUN}" == "true" ]]; then
    echo "DRY-RUN: gatk-sv-ploidy preprocess -i \"${INPUT_DEPTH}\" -o \"${PREPROCESS_DIR}\" ${SD_ARGS} ${PR_ARGS}"
else
    gatk-sv-ploidy preprocess \
        -i "${INPUT_DEPTH}" \
        -o "${PREPROCESS_DIR}" \
        $SD_ARGS \
        $PR_ARGS \
        $PREPROCESS_ARGS
fi

# ── step 2: infer ────────────────────────────────────────────────────────────
echo "[2/${TOTAL_STEPS}] infer"
AF_ARGS=""
if [[ -f "${SITE_DATA}" ]]; then
    AF_ARGS="--site-data ${SITE_DATA}"
fi
if [[ "${DRY_RUN}" == "true" ]]; then
    echo "DRY-RUN: gatk-sv-ploidy infer -i \"${PREPROCESSED_DEPTH}\" -o \"${INFER_DIR}\" ${AF_ARGS}"
else
    gatk-sv-ploidy infer \
        -i "${PREPROCESSED_DEPTH}" \
        -o "${INFER_DIR}" \
        $AF_ARGS \
        $MODEL_ARGS
fi

# ── step 3: call ─────────────────────────────────────────────────────────────
echo "[3/${TOTAL_STEPS}] call"
if [[ "${DRY_RUN}" == "true" ]]; then
    echo "DRY-RUN: gatk-sv-ploidy call -c \"${CHROM_STATS}\" -o \"${CALL_DIR}\""
else
    gatk-sv-ploidy call \
        -c "${CHROM_STATS}" \
        -o "${CALL_DIR}"
fi

# ── step 4: posterior predictive checks (optional) ──────────────────────────
PPD_ARGS=""
if [[ -f "${SITE_DATA}" ]]; then
    PPD_ARGS="--site-data ${SITE_DATA}"
fi
if [[ "${ENABLE_PPD}" == "true" ]]; then
    echo "[4/${TOTAL_STEPS}] ppd"
    if [[ "${DRY_RUN}" == "true" ]]; then
        echo "DRY-RUN: gatk-sv-ploidy ppd -i \"${PREPROCESSED_DEPTH}\" -a \"${INFERENCE_ARTIFACTS}\" -o \"${PPD_DIR}\" ${PPD_ARGS}"
    else
        gatk-sv-ploidy ppd \
            -i "${PREPROCESSED_DEPTH}" \
            -a "${INFERENCE_ARTIFACTS}" \
            -o "${PPD_DIR}" \
            $PPD_ARGS
    fi
fi

# ── step 5: plot ─────────────────────────────────────────────────────────────
if [[ "${ENABLE_PPD}" == "true" ]]; then
    echo "[5/${TOTAL_STEPS}] plot"
else
    echo "[4/${TOTAL_STEPS}] plot"
fi
PLOT_PPD_ARGS=""
if [[ "${ENABLE_PPD}" == "true" ]]; then
    PLOT_PPD_ARGS="--ppd-bin-summary ${PPD_BIN_SUMMARY} --ppd-chr-summary ${PPD_CHR_SUMMARY}"
fi
if [[ "${DRY_RUN}" == "true" ]]; then
    echo "DRY-RUN: gatk-sv-ploidy plot -c \"${CHROM_STATS}\" -b \"${BIN_STATS}\" -t \"${TRAINING_LOSS}\" -s \"${PREDICTIONS}\" -o \"${PLOT_DIR}\" ${PLOT_PPD_ARGS}"
else
    # Include highlight sample if provided
    if [[ -n "${HIGHLIGHT_SAMPLE}" ]]; then
        gatk-sv-ploidy plot \
            -c "${CHROM_STATS}" \
            -b "${BIN_STATS}" \
            -t "${TRAINING_LOSS}" \
            -s "${PREDICTIONS}" \
            -o "${PLOT_DIR}" \
            --highlight-sample "${HIGHLIGHT_SAMPLE}" \
            $PLOT_PPD_ARGS \
            $PLOT_ARGS
    else
        gatk-sv-ploidy plot \
            -c "${CHROM_STATS}" \
            -b "${BIN_STATS}" \
            -t "${TRAINING_LOSS}" \
            -s "${PREDICTIONS}" \
            -o "${PLOT_DIR}" \
            $PLOT_PPD_ARGS \
            $PLOT_ARGS
    fi
fi

# ── step 6: eval (optional — skipped if TRUTH_JSON is unset) ─────────────────
if [[ -n "${TRUTH_JSON}" ]]; then
    if [[ "${ENABLE_PPD}" == "true" ]]; then
        echo "[6/${TOTAL_STEPS}] eval"
    else
        echo "[5/${TOTAL_STEPS}] eval"
    fi
    if [[ "${DRY_RUN}" == "true" ]]; then
        echo "DRY-RUN: gatk-sv-ploidy eval -p \"${PREDICTIONS}\" -t \"${TRUTH_JSON}\" -o \"${EVAL_DIR}\""
    else
        gatk-sv-ploidy eval \
            -p "${PREDICTIONS}" \
            -t "${TRUTH_JSON}" \
            -o "${EVAL_DIR}"
    fi
else
    if [[ "${ENABLE_PPD}" == "true" ]]; then
        echo "[6/${TOTAL_STEPS}] eval  (skipped — set TRUTH_JSON to enable)"
    else
        echo "[5/${TOTAL_STEPS}] eval  (skipped — set TRUTH_JSON to enable)"
    fi
fi

echo ""
echo "=== pipeline complete ==="
echo "  Outputs in: ${WORK_DIR}"
