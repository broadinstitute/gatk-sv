#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage: run_gd.sh --work-dir DIR --input-depth FILE --high-res-depth FILE --baf-table FILE --gd-table FILE \
                 --segdup-bed FILE --centromere-bed FILE --acrocentric-arm-bed FILE --gaps-bed FILE \
                 --gtf FILE --transition-matrix FILE --breakpoint-transition-matrix FILE [options]

Required arguments:
  --work-dir DIR
  --input-depth FILE
  --high-res-depth FILE
  --baf-table FILE
  --gd-table FILE
  --segdup-bed FILE
  --centromere-bed FILE
  --acrocentric-arm-bed FILE
  --gaps-bed FILE
  --gtf FILE
  --transition-matrix FILE
  --breakpoint-transition-matrix FILE

Optional arguments:
  --flank-exclusion-interval FILE
  --flank-exclusion-intervals FILE [FILE ...]
  --truth-table FILE
  --preprocess-args STRING
  --infer-args STRING
  --call-args STRING
  --eval-args STRING
  --gd-cmd CMD

For backwards compatibility, a single positional argument is accepted as --work-dir.
EOF
    exit 1
}

require_arg() {
    local value="$1"
    local flag="$2"
    if [[ -z "${value}" ]]; then
        echo "Missing required argument: ${flag}" >&2
        usage
    fi
}

WORK_DIR=""
INPUT_DEPTH=""
HIGH_RESOLUTION_DEPTH=""
BAF_TABLE=""
GD_TABLE=""
SEG_DUP_BED=""
CENTROMERE_BED=""
ACROCENTRIC_ARM_BED=""
GAPS_BED=""
GTF=""
TRANSITION_MATRIX=""
BREAKPOINT_TRANSITION_MATRIX=""
TRUTH_TABLE=""

PREPROCESS_ARGS=""
INFER_ARGS=""
CALL_ARGS=""
EVAL_ARGS=""

FLANK_EXCLUSION_INTERVALS=()
POSITIONAL_ARGS=()

GD_CMD=("gatk-sv-gd")

while [[ $# -gt 0 ]]; do
    case "$1" in
        --work-dir)
            WORK_DIR="$2"
            shift 2
            ;;
        --input-depth)
            INPUT_DEPTH="$2"
            shift 2
            ;;
        --high-res-depth|--high-res-counts)
            HIGH_RESOLUTION_DEPTH="$2"
            shift 2
            ;;
        --baf-table)
            BAF_TABLE="$2"
            shift 2
            ;;
        --gd-table)
            GD_TABLE="$2"
            shift 2
            ;;
        --segdup-bed)
            SEG_DUP_BED="$2"
            shift 2
            ;;
        --centromere-bed)
            CENTROMERE_BED="$2"
            shift 2
            ;;
        --acrocentric-arm-bed)
            ACROCENTRIC_ARM_BED="$2"
            shift 2
            ;;
        --flank-exclusion-interval)
            FLANK_EXCLUSION_INTERVALS+=("$2")
            shift 2
            ;;
        --flank-exclusion-intervals)
            shift
            while [[ $# -gt 0 && "$1" != --* ]]; do
                FLANK_EXCLUSION_INTERVALS+=("$1")
                shift
            done
            ;;
        --gaps-bed)
            GAPS_BED="$2"
            shift 2
            ;;
        --gtf)
            GTF="$2"
            shift 2
            ;;
        --transition-matrix)
            TRANSITION_MATRIX="$2"
            shift 2
            ;;
        --breakpoint-transition-matrix)
            BREAKPOINT_TRANSITION_MATRIX="$2"
            shift 2
            ;;
        --truth-table)
            TRUTH_TABLE="$2"
            shift 2
            ;;
        --preprocess-args)
            PREPROCESS_ARGS="$2"
            shift 2
            ;;
        --infer-args)
            INFER_ARGS="$2"
            shift 2
            ;;
        --call-args)
            CALL_ARGS="$2"
            shift 2
            ;;
        --eval-args)
            EVAL_ARGS="$2"
            shift 2
            ;;
        --gd-cmd)
            GD_CMD=("$2")
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        --)
            shift
            POSITIONAL_ARGS+=("$@")
            break
            ;;
        -*)
            echo "Unknown option: $1" >&2
            usage
            ;;
        *)
            POSITIONAL_ARGS+=("$1")
            shift
            ;;
    esac
done

if [[ -z "${WORK_DIR}" && ${#POSITIONAL_ARGS[@]} -eq 1 ]]; then
    WORK_DIR="${POSITIONAL_ARGS[0]}"
elif [[ ${#POSITIONAL_ARGS[@]} -gt 0 ]]; then
    echo "Unexpected positional arguments: ${POSITIONAL_ARGS[*]}" >&2
    usage
fi

require_arg "${WORK_DIR}" "--work-dir"
require_arg "${INPUT_DEPTH}" "--input-depth"
require_arg "${HIGH_RESOLUTION_DEPTH}" "--high-res-depth"
require_arg "${BAF_TABLE}" "--baf-table"
require_arg "${GD_TABLE}" "--gd-table"
require_arg "${SEG_DUP_BED}" "--segdup-bed"
require_arg "${CENTROMERE_BED}" "--centromere-bed"
require_arg "${ACROCENTRIC_ARM_BED}" "--acrocentric-arm-bed"
require_arg "${GAPS_BED}" "--gaps-bed"
require_arg "${GTF}" "--gtf"
require_arg "${TRANSITION_MATRIX}" "--transition-matrix"
require_arg "${BREAKPOINT_TRANSITION_MATRIX}" "--breakpoint-transition-matrix"

# ── Directories ────────────────────────────────────────────────────────────
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
VITERBI_PATHS="${CALL_DIR}/viterbi_paths.tsv.gz"
EVENT_MARGINALS="${CALL_DIR}/event_marginals.tsv.gz"
PLOIDY_TABLE="${PREPROCESS_DIR}/ploidy_estimates.tsv"
PREPROCESSED_BAF="${PREPROCESS_DIR}/preprocessed_baf.tsv.gz"
EVAL_REPORT="${EVAL_DIR}/truth_evaluation_report.tsv"

mkdir -p "${WORK_DIR}"

# ── Step 1: preprocess ─────────────────────────────────────────────────────
echo "[1/5] preprocess"
rm -rf "${PREPROCESS_DIR}"
PREPROCESS_CMD=(
    "${GD_CMD[@]}"
    preprocess
    -i "${INPUT_DEPTH}" \
    --high-res-counts "${HIGH_RESOLUTION_DEPTH}" \
    --baf-table "${BAF_TABLE}" \
    -g "${GD_TABLE}" \
    -o "${PREPROCESS_DIR}" \
    -e "${SEG_DUP_BED}" \
    -e "${CENTROMERE_BED}" \
    -e "${ACROCENTRIC_ARM_BED}" \
    --verbose \
)

if [[ ${#FLANK_EXCLUSION_INTERVALS[@]} -gt 0 ]]; then
    PREPROCESS_CMD+=(--flank-exclusion-intervals "${FLANK_EXCLUSION_INTERVALS[@]}")
fi

if [[ -n "${PREPROCESS_ARGS}" ]]; then
    # shellcheck disable=SC2206
    PREPROCESS_EXTRA=( ${PREPROCESS_ARGS} )
    PREPROCESS_CMD+=("${PREPROCESS_EXTRA[@]}")
fi

"${PREPROCESS_CMD[@]}"

# ── Step 2: infer ──────────────────────────────────────────────────────────
echo "[2/5] infer"
rm -rf "${INFER_DIR}"
INFER_CMD=(
    "${GD_CMD[@]}"
    infer
    --preprocessed-dir "${PREPROCESS_DIR}" \
    -o "${INFER_DIR}" \
    --verbose \
)

if [[ -n "${INFER_ARGS}" ]]; then
    # shellcheck disable=SC2206
    INFER_EXTRA=( ${INFER_ARGS} )
    INFER_CMD+=("${INFER_EXTRA[@]}")
fi

"${INFER_CMD[@]}"

# ── Step 3: call ───────────────────────────────────────────────────────────
echo "[3/5] call"
rm -rf "${CALL_DIR}"
CALL_CMD=(
    "${GD_CMD[@]}"
    call
    --cn-posteriors "${CN_POSTERIORS}" \
    --bin-mappings "${BIN_MAPPINGS}" \
    -g "${FILTERED_GD_TABLE}" \
    -o "${CALL_DIR}" \
    --ploidy-table "${PLOIDY_TABLE}" \
    --transition-matrix "${TRANSITION_MATRIX}" \
    --breakpoint-transition-matrix "${BREAKPOINT_TRANSITION_MATRIX}" \
    --verbose \
)

if [[ -n "${CALL_ARGS}" ]]; then
    # shellcheck disable=SC2206
    CALL_EXTRA=( ${CALL_ARGS} )
    CALL_CMD+=("${CALL_EXTRA[@]}")
fi

"${CALL_CMD[@]}"

# ── Step 4: eval (optional) ────────────────────────────────────────────────
if [[ -n "${TRUTH_TABLE}" ]]; then
    echo "[4/5] eval"
    rm -rf "${EVAL_DIR}"
    EVAL_CMD=(
        "${GD_CMD[@]}"
        eval
        --calls "${GD_CALLS}" \
        --truth-table "${TRUTH_TABLE}" \
        --gd-table "${FILTERED_GD_TABLE}" \
        --ploidy-table "${PLOIDY_TABLE}" \
        -o "${EVAL_DIR}" \
    )

    if [[ -n "${EVAL_ARGS}" ]]; then
        # shellcheck disable=SC2206
        EVAL_EXTRA=( ${EVAL_ARGS} )
        EVAL_CMD+=("${EVAL_EXTRA[@]}")
    fi

    "${EVAL_CMD[@]}"
else
    echo "[4/5] eval  (skipped — TRUTH_TABLE not set)"
fi

# ── Step 5: plot ───────────────────────────────────────────────────────────
echo "[5/5] plot"
rm -rf "${PLOT_DIR}"

PLOT_CMD=(
    "${GD_CMD[@]}"
    plot
    --calls "${GD_CALLS}"
    --cn-posteriors "${CN_POSTERIORS}"
    --raw-counts "${INPUT_DEPTH}"
    --high-res-counts "${HIGH_RESOLUTION_DEPTH}"
    -g "${FILTERED_GD_TABLE}"
    -o "${PLOT_DIR}"
    --gaps-bed "${GAPS_BED}"
    --gtf "${GTF}"
    --segdup-bed "${SEG_DUP_BED}"
    --ploidy-table "${PLOIDY_TABLE}"
    --viterbi-paths "${VITERBI_PATHS}"
    --event-marginals "${EVENT_MARGINALS}"
)

if [[ -n "${TRUTH_TABLE}" ]]; then
    PLOT_CMD+=(--eval-report "${EVAL_REPORT}")
fi

"${PLOT_CMD[@]}"

echo
if [[ -n "${TRUTH_TABLE}" ]]; then
    echo "Eval report: ${EVAL_REPORT}"
else
    echo "Eval report: skipped (TRUTH_TABLE not set)"
fi
echo "Plot directory: ${PLOT_DIR}"
