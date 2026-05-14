from __future__ import annotations

from pathlib import Path


def test_run_ploidy_wrapper_exposes_ppd_args_passthrough() -> None:
    script_path = Path(__file__).resolve().parents[2] / "run_ploidy.sh"
    script_text = script_path.read_text(encoding="utf-8")

    assert '# run_ploidy.sh — run the gatk-sv-ploidy wrapper pipeline' in script_text
    assert '#   preprocess -> polyploidy (if site_data.npz exists) -> infer -> call ->' in script_text
    assert '# Required arguments for a normal run:' in script_text
    assert '#   --input-depth : Raw bins x samples depth TSV consumed by preprocess' in script_text
    assert '--ppd-args STRING        Extra arguments passed through to the ppd step' in script_text
    assert '--ppd-args)' in script_text
    assert 'PPD_STEP_ARGS="$2"; shift 2;;' in script_text
    assert '--ppd-args=*)' in script_text
    assert 'PPD_STEP_ARGS="${1#*=}"; shift;;' in script_text
    assert 'PPD extra args   :' in script_text
    assert '$PPD_INPUT_ARGS \\' in script_text
    assert '$PPD_STEP_ARGS' in script_text


def test_run_ploidy_wrapper_uses_step_passthroughs_for_infer_call_and_preprocess() -> None:
    script_path = Path(__file__).resolve().parents[2] / "run_ploidy.sh"
    script_text = script_path.read_text(encoding="utf-8")

    assert 'Usage: $0 --input-depth PATH --work-dir DIR [options]' in script_text
    assert 'Runs preprocess, infer, call, plot, and optional ppd/eval under WORK_DIR' in script_text
    assert '--preprocess-args STRING Extra arguments passed through to the preprocess step' in script_text
    assert '--infer-args STRING      Extra arguments passed through to the infer step' in script_text
    assert '--infer-args)' in script_text
    assert 'INFER_ARGS="$2"; shift 2;;' in script_text
    assert '--infer-args=*)' in script_text
    assert 'INFER_ARGS="${1#*=}"; shift;;' in script_text
    assert '--polyploidy-args)' in script_text
    assert '--polyploidy-args=*)' in script_text
    assert 'Baseline CN tool : polyploidy' in script_text
    assert 'run_cli polyploidy \\' in script_text
    assert '--min-poor-region-coverage)' not in script_text
    assert '--use-callq20)' in script_text
    assert 'USE_CALLQ20="false"' in script_text
    assert 'CALL_STEP_ARGS="--use-callq20 ${CALL_STEP_ARGS}"' in script_text
    assert 'PR_ARGS="--poor-regions ${POOR_REGIONS}"' in script_text
    assert 'Error: --input-depth is required' in script_text


def test_run_ploidy_wrapper_recomputes_polyploidy_after_preprocess() -> None:
    script_path = Path(__file__).resolve().parents[2] / "run_ploidy.sh"
    script_text = script_path.read_text(encoding="utf-8")

    assert 'RUN_POLYPLOIDY_STATUS="pending site_data check"' in script_text
    assert 'echo "  Run polyploidy   : ${RUN_POLYPLOIDY_STATUS}"' in script_text
    assert 'RUN_POLYPLOIDY="false"' in script_text
    assert 'if [[ -f "${SITE_DATA}" ]]; then' in script_text
    assert 'elif [[ "${DRY_RUN}" == "true" && -n "${SITE_DEPTH_LIST}" ]]; then' in script_text
    assert 'Warning: site-depth list was provided but preprocess did not create ${SITE_DATA}; AF-enabled polyploidy and raw AF plots will be skipped.' in script_text
