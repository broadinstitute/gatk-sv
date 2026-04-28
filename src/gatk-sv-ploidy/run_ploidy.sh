#!/usr/bin/env bash
# run_ploidy.sh — run the gatk-sv-ploidy pipeline from preprocess outputs
#
# Usage:
#   ./run_ploidy.sh --work-dir DIR [options]
#
# Required constants (configured via CLI args):
#   --truth-json   : JSON mapping sample ID → true aneuploidy type (optional)
#   --work-dir     : Root directory for all pipeline outputs (required)
#
# Note: the preprocess block is intentionally commented out in this wrapper.
# It currently expects existing preprocess outputs under ${WORK_DIR}/preprocess/.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_SRC="${SCRIPT_DIR}/src"
VENV_PYTHON="${SCRIPT_DIR}/.venv/bin/python"

if [[ -x "${VENV_PYTHON}" ]]; then
    CLI_PYTHON="${VENV_PYTHON}"
elif command -v python3 >/dev/null 2>&1; then
    CLI_PYTHON="$(command -v python3)"
else
    echo "Error: could not find a Python interpreter for gatk-sv-ploidy" >&2
    exit 1
fi

run_cli() {
    PYTHONPATH="${PACKAGE_SRC}${PYTHONPATH:+:${PYTHONPATH}}" \
        "${CLI_PYTHON}" -m gatk_sv_ploidy.cli "$@"
}

describe_redacted_input() {
    local value="${1:-}"
    local present_label="${2:-provided (redacted)}"
    local absent_label="${3:-not provided}"

    if [[ -n "${value}" ]]; then
        printf '%s\n' "${present_label}"
    else
        printf '%s\n' "${absent_label}"
    fi
}

dry_run_step() {
    local step_name="$1"
    echo "DRY-RUN: ${step_name} command prepared (arguments redacted)"
}

INPUT_DEPTH=""
TRUTH_JSON=""
WORK_DIR=""
DRY_RUN="false"
ENABLE_PPD="false"
PREPROCESS_ARGS=""
MODEL_ARGS=""
CALL_ARGS=""
PLOT_ARGS=""
HIGHLIGHT_SAMPLE=""
SITE_DEPTH_LIST=""
POOR_REGIONS=""
MIN_POOR_REGION_COVERAGE="0.5"

usage() {
    echo "Usage: $0 --work-dir DIR [--input-depth PATH] [--truth-json PATH] [--site-depth-list PATH] [--poor-regions PATH] [--ppd]" >&2
    echo "  Note: preprocess is currently disabled in this wrapper; it expects existing outputs in WORK_DIR/preprocess/" >&2
    echo "  --input-depth PATH       Optional raw bins×samples depth TSV retained for provenance / future preprocess runs" >&2
    echo "  --work-dir DIR           Root directory for all pipeline outputs" >&2
    echo "  --truth-json PATH        Optional truth JSON to enable evaluation" >&2
    echo "  --site-depth-list PATH   Optional file listing per-sample SD file paths (one per line)" >&2
    echo "  --poor-regions PATH      Optional BED of poor regions to remove during preprocess" >&2
    echo "  --min-poor-region-coverage FLOAT  Min fraction of a bin overlapped by poor regions to filter it" >&2
    echo "  --call-args STRING       Extra arguments passed through to the call step" >&2
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
        --call-args)
            CALL_ARGS="$2"; shift 2;;
        --call-args=*)
            CALL_ARGS="${1#*=}"; shift;;
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

if [[ -z "${WORK_DIR}" ]]; then
    echo "Error: --work-dir is required" >&2
    usage
fi

CALL_NEEDS_PPD_FILTER="false"
if [[ "${CALL_ARGS}" == *"--min-binq"* ]]; then
    CALL_NEEDS_PPD_FILTER="true"
fi
if [[ "${CALL_NEEDS_PPD_FILTER}" == "true" && "${ENABLE_PPD}" != "true" ]]; then
    echo "Error: --call-args containing --min-binq requires --ppd so ppd_bin_quality.tsv is available" >&2
    exit 1
fi

if [[ "${DRY_RUN}" == "true" ]]; then
    echo "DRY-RUN mode: commands will be printed but not executed"
fi

TOTAL_STEPS=4
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
PPD_BIN_QUALITY="${PPD_DIR}/ppd_bin_quality.tsv"
PPD_CHR_SUMMARY="${PPD_DIR}/ppd_chromosome_summary.tsv"
IGNORED_BINS="${CALL_DIR}/ignored_bins.tsv.gz"

echo "=== gatk-sv-ploidy pipeline ==="
echo "  Input depth      : $(describe_redacted_input "${INPUT_DEPTH}")"
echo "  Site depth list  : $(describe_redacted_input "${SITE_DEPTH_LIST}")"
echo "  Poor regions     : $(describe_redacted_input "${POOR_REGIONS}")"
echo "  Poor region cov  : ${MIN_POOR_REGION_COVERAGE}"
echo "  Truth JSON       : $(describe_redacted_input "${TRUTH_JSON}")"
echo "  Highlight sample : $(describe_redacted_input "${HIGHLIGHT_SAMPLE}")"
echo "  Run PPD          : ${ENABLE_PPD}"
echo "  Work dir         : configured (redacted)"
echo "  Python           : resolved"
echo ""

CALL_STEP_ARGS="${CALL_ARGS}"
if [[ "${CALL_NEEDS_PPD_FILTER}" == "true" ]]; then
    if [[ "${CALL_STEP_ARGS}" != *"--bin-stats"* ]]; then
        CALL_STEP_ARGS="--bin-stats ${BIN_STATS} ${CALL_STEP_ARGS}"
    fi
    if [[ "${CALL_STEP_ARGS}" != *"--ppd-bin-quality"* ]]; then
        CALL_STEP_ARGS="--ppd-bin-quality ${PPD_BIN_QUALITY} ${CALL_STEP_ARGS}"
    fi
fi

# ── preprocess ─────────────────────────
# echo "[0/${TOTAL_STEPS}] preprocess"
# SD_ARGS=""
# if [[ -n "${SITE_DEPTH_LIST}" ]]; then
#     SD_ARGS="--site-depth-list ${SITE_DEPTH_LIST}"
# fi
# PR_ARGS=""
# if [[ -n "${POOR_REGIONS}" ]]; then
#     PR_ARGS="--poor-regions ${POOR_REGIONS} --min-poor-region-coverage ${MIN_POOR_REGION_COVERAGE}"
# fi
# if [[ "${DRY_RUN}" == "true" ]]; then
#     dry_run_step "preprocess"
# else
#     run_cli preprocess \
#         -i "${INPUT_DEPTH}" \
#         -o "${PREPROCESS_DIR}" \
#         $SD_ARGS \
#         $PR_ARGS \
#         $PREPROCESS_ARGS
# fi

# ── step 1: infer ────────────────────────────────────────────────────────────
echo "[1/${TOTAL_STEPS}] infer"
AF_ARGS=""
if [[ -f "${SITE_DATA}" ]]; then
    AF_ARGS="--site-data ${SITE_DATA}"
fi
if [[ "${DRY_RUN}" == "true" ]]; then
    dry_run_step "infer"
else
    run_cli infer \
        -i "${PREPROCESSED_DEPTH}" \
        -o "${INFER_DIR}" \
        $AF_ARGS \
        $MODEL_ARGS
fi

# ── step 2/3: posterior predictive checks and call ──────────────────────────
PPD_ARGS=""
if [[ -f "${SITE_DATA}" ]]; then
    PPD_ARGS="--site-data ${SITE_DATA}"
fi
CALL_WITH_TRUTH_ARGS="${CALL_STEP_ARGS}"
if [[ -n "${TRUTH_JSON}" ]]; then
    CALL_WITH_TRUTH_ARGS="--truth-json ${TRUTH_JSON} ${CALL_WITH_TRUTH_ARGS}"
fi

run_call_step() {
    if [[ "${DRY_RUN}" == "true" ]]; then
        dry_run_step "call"
    else
        run_cli call \
            -c "${CHROM_STATS}" \
            -o "${CALL_DIR}" \
            $CALL_WITH_TRUTH_ARGS
    fi
}

run_ppd_step() {
    if [[ "${DRY_RUN}" == "true" ]]; then
        dry_run_step "ppd"
    else
        run_cli ppd \
            -i "${PREPROCESSED_DEPTH}" \
            -a "${INFERENCE_ARTIFACTS}" \
            -o "${PPD_DIR}" \
            $PPD_ARGS
    fi
}

if [[ "${ENABLE_PPD}" == "true" && "${CALL_NEEDS_PPD_FILTER}" == "true" ]]; then
    echo "[2/${TOTAL_STEPS}] ppd"
    run_ppd_step
    echo "[3/${TOTAL_STEPS}] call"
    run_call_step
else
    echo "[2/${TOTAL_STEPS}] call"
    run_call_step
    if [[ "${ENABLE_PPD}" == "true" ]]; then
        echo "[3/${TOTAL_STEPS}] ppd"
        run_ppd_step
    fi
fi

# ── step 3/4: plot ───────────────────────────────────────────────────────────
if [[ "${ENABLE_PPD}" == "true" ]]; then
    echo "[4/${TOTAL_STEPS}] plot"
else
    echo "[3/${TOTAL_STEPS}] plot"
fi
PLOT_PPD_ARGS=""
if [[ "${ENABLE_PPD}" == "true" ]]; then
    PLOT_PPD_ARGS="--ppd-bin-summary ${PPD_BIN_SUMMARY} --ppd-chr-summary ${PPD_CHR_SUMMARY}"
fi
PLOT_BINQ_ARGS=""
if [[ "${ENABLE_PPD}" == "true" ]]; then
    PLOT_BINQ_ARGS="--ppd-bin-quality ${PPD_BIN_QUALITY}"
elif [[ -f "${PPD_BIN_QUALITY}" ]]; then
    PLOT_BINQ_ARGS="--ppd-bin-quality ${PPD_BIN_QUALITY}"
fi
PLOT_CHROM_STATS="${CHROM_STATS}"
PLOT_IGNORED_ARGS=""
if [[ "${CALL_NEEDS_PPD_FILTER}" == "true" ]]; then
    PLOT_CHROM_STATS="${CALL_DIR}/chromosome_stats.filtered.tsv"
    PLOT_IGNORED_ARGS="--ignored-bins ${IGNORED_BINS}"
elif [[ -f "${CALL_DIR}/chromosome_stats.filtered.tsv" ]]; then
    PLOT_CHROM_STATS="${CALL_DIR}/chromosome_stats.filtered.tsv"
fi
if [[ -f "${IGNORED_BINS}" ]]; then
    PLOT_IGNORED_ARGS="--ignored-bins ${IGNORED_BINS}"
fi
PLOT_SITE_ARGS=""
if [[ -f "${SITE_DATA}" ]]; then
    PLOT_SITE_ARGS="--site-data ${SITE_DATA}"
fi
if [[ "${DRY_RUN}" == "true" ]]; then
    dry_run_step "plot"
else
    # Include highlight sample if provided
    if [[ -n "${HIGHLIGHT_SAMPLE}" ]]; then
        run_cli plot \
            -c "${PLOT_CHROM_STATS}" \
            -b "${BIN_STATS}" \
            -t "${TRAINING_LOSS}" \
            -s "${PREDICTIONS}" \
            -o "${PLOT_DIR}" \
            $PLOT_SITE_ARGS \
            $PLOT_IGNORED_ARGS \
            $PLOT_BINQ_ARGS \
            --highlight-sample "${HIGHLIGHT_SAMPLE}" \
            $PLOT_PPD_ARGS \
            $PLOT_ARGS
    else
        run_cli plot \
            -c "${PLOT_CHROM_STATS}" \
            -b "${BIN_STATS}" \
            -t "${TRAINING_LOSS}" \
            -s "${PREDICTIONS}" \
            -o "${PLOT_DIR}" \
            $PLOT_SITE_ARGS \
            $PLOT_IGNORED_ARGS \
            $PLOT_BINQ_ARGS \
            $PLOT_PPD_ARGS \
            $PLOT_ARGS
    fi
fi

# ── step 4/5: eval (optional — skipped if TRUTH_JSON is unset) ──────────────
if [[ -n "${TRUTH_JSON}" ]]; then
    if [[ "${ENABLE_PPD}" == "true" ]]; then
        echo "[5/${TOTAL_STEPS}] eval"
    else
        echo "[4/${TOTAL_STEPS}] eval"
    fi
    if [[ "${DRY_RUN}" == "true" ]]; then
        dry_run_step "eval"
    else
        run_cli eval \
            -p "${PREDICTIONS}" \
            -t "${TRUTH_JSON}" \
            -o "${EVAL_DIR}"
    fi
else
    if [[ "${ENABLE_PPD}" == "true" ]]; then
        echo "[5/${TOTAL_STEPS}] eval  (skipped — set TRUTH_JSON to enable)"
    else
        echo "[4/${TOTAL_STEPS}] eval  (skipped — set TRUTH_JSON to enable)"
    fi
fi

echo ""
echo "=== pipeline complete ==="
echo "  Outputs written under configured work directory"
