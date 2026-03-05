# GATK-SV Copilot Instructions

## Project Overview
Cloud-native structural variation discovery pipeline for Illumina WGS data.
Orchestrated in **WDL** (executed via Cromwell/Terra), with tools in **Python**, **R**, **Bash**, and a small amount of **Java**. All tools ship in Docker images pushed to GCR.

## Repository Layout
- `wdl/` — ~120 WDL files (flat directory, no subdirs). Workflows are PascalCase (`ClusterBatch.wdl`); task libraries are prefixed `Tasks*` (`TasksMakeCohortVcf.wdl`). Shared utilities: `Utils.wdl`, `Structs.wdl`.
- `src/` — Python packages and scripts:
  - `svtk/` — Cython-enabled SV toolkit (pysam, pybedtools). `setup.py`.
  - `sv_utils/` — General SV utilities with pytest suite. `setup.py`, Python >3.8.
  - `svtest/`, `svqc/` — Metrics and QC packages. `setup.py`.
  - `sv-pipeline/` — Standalone scripts (not pip-installable) organized by pipeline stage (`00_preprocessing/` … `05_annotation/`), plus `scripts/` used directly by WDL tasks.
  - `gatk-sv-ploidy/` — Modern `pyproject.toml` package (Python ≥3.9, PyTorch/Pyro). CLI: `gatk-sv-ploidy`.
  - `denovo/`, `str/`, `stripy/` — Standalone analysis scripts.
  - `RdTest/`, `WGD/` — R/shell scripts for read-depth testing and dosage scoring (rarely modified).
  - `gatk-sv-gd/` — Ignore; work-in-progress from another branch.
- `dockerfiles/` — Layered Docker images. Base chain: `sv-base-mini` → `samtools-cloud` → `sv-base` → `sv-pipeline`. Virtual-env pattern splits build deps from runtime.
- `inputs/` — Jinja2 template system: `values/` (JSON with docker tags, references, sample data) → `templates/` (`.json.tmpl`) → `build/` (generated inputs per test batch).
- `scripts/docker/` — `build_docker.py` for image builds with git-diff-based change detection.
- `website/` — Docusaurus documentation site (published at https://broadinstitute.github.io/gatk-sv/). Source is `website/docs/` (Markdown). Per-module docs live in `website/docs/modules/`. Build locally with `cd website && npm install && npm run start`. **Code changes that affect pipeline behaviour, inputs/outputs, or runtime attributes should be reflected in the relevant `website/docs/` pages.**

## WDL Conventions
- WDL version: `version 1.0`. Imports use relative paths with aliases: `import "Utils.wdl" as util`.
- **Inputs**: `snake_case`, type-suffixed (`_vcf`, `_bed`, `_index`, `_file`, `_list`, `_tar`). Docker images passed as `String` inputs (e.g. `String sv_pipeline_docker`), never hardcoded.
- **RuntimeAttr pattern**: Tasks accept `RuntimeAttr? runtime_attr_override`; workflows expose per-task overrides as `RuntimeAttr? runtime_attr_<taskname>`. Defaults are defined inline; `select_first` merges override with default.
- **Bash in `command` blocks**: Start with `set -eu` or `set -euo pipefail`. Use `~{var}` interpolation. Optional args: `~{"--flag " + optional_var}` or `~{true="--enable" false="" bool_var}`.
- **Calling imported tasks**: `module.TaskName { input: ... }`. Alias with `as` when calling the same task multiple times (e.g. `call pesr.ClusterPESR as ClusterPESR_dragen`).
- **Nested workflow inputs** in JSON use dot-path notation: `"WorkflowName.SubWorkflow.input_name"`.

## Python Conventions
- Style: PEP-8. Linter: flake8 via `tox -e lint` (ignores E501 line length, W504, E722 bare except, F405 star imports).
- Packages use either `setup.py` (older: svtk, sv_utils, svtest, svqc) or `pyproject.toml` (newer: gatk-sv-ploidy). Install editable: `pip install -e .`.
- Test framework: pytest. `sv_utils` has the most complete test suite (`src/sv_utils/tests/`). Run with `pytest` from the package directory.
- Common dependencies across packages: numpy, pandas, pysam, scipy.

## Docker Images
- Images are layered; changing a base image requires rebuilding dependents. `scripts/docker/build_docker.py` handles this automatically.
- Image tags in `inputs/values/dockers.json` follow format `image:branch-or-date-shortsha`.
- After building new images, update `dockers.json` with new tags.
- **Do not trigger Docker builds automatically.** Builds are done manually on Linux VMs (macOS compatibility issues). Only build when the user explicitly requests it.

## Input Generation
```bash
scripts/inputs/build_inputs.py inputs/values inputs/templates/test \
  inputs/build/ref_panel_1kg/test -a '{ "test_batch" : "ref_panel_1kg" }'
```
Templates use Jinja2 with values from `inputs/values/*.json`.

## CI/CD (`.github/workflows/`)
- `pytest.yaml` — flake8 lint on push/PR to `main`.
- `sv_pipeline_docker.yaml` — Auto-builds Docker images on changes to `src/` or `dockerfiles/`.
- `wdl_tests.yaml` — Validates WDL syntax (miniwdl) and inputs (WOMtool) on changes to `wdl/` or `inputs/`. These checks are slow; only run locally when explicitly asked.

## Pipeline Modes
- **Joint calling (multi-sample)**: Spread across many WDLs starting from `GATKSVPipelineBatch.wdl` / `GATKSVPipelinePhase1.wdl`. Processes samples in batches then merges cohort-wide.
- **Single-sample mode**: Contained in `wdl/GATKSVPipelineSingleSample.wdl`; runs the full pipeline for one sample.
- **sv-shell** (`src/sv_shell/`): Bash reimplementation of the single-sample pipeline designed to be portable (no Cromwell/Terra required). Mirrors the WDL logic but runs locally or on simple Linux environments.

## SV VCF Formatting
- Two VCF conventions exist: **svtk style** (used by upstream tools) and **GATK style** (used after `FormatVcfForGatk`). The pipeline converts between them at specific stages.
- Key INFO tags: `END` (end position on same chrom), `CHR2` (mate chromosome), `END2` (position on mate chromosome). The pipeline was recently unified to use a single convention for these tags—follow whatever the surrounding code does.
- SV types: DEL, DUP, INV, INS, BND, CPX, CNV, CTX. Multi-allelic and complex (CPX) variants have special handling throughout.

## Change Propagation Order
Changes often span multiple layers of the codebase. Apply updates in this order, as each step may depend on the previous:

1. **`src/`** — Implement logic changes in Python scripts or packages first.
2. **`wdl/`** — Update WDL tasks/workflows to invoke the new/changed scripts or behaviour.
3. **`inputs/`** — If WDL inputs/outputs changed, update templates (`inputs/templates/`) and regenerate build inputs.
4. **`dockerfiles/`** — If new runtime dependencies were introduced, update the relevant Dockerfile(s). Then rebuild images manually (see Docker Images section).
5. **`website/`** — Update `website/docs/` documentation to reflect any user-facing changes to pipeline behaviour, inputs, or outputs.

## Key Patterns to Follow
- All cloud resource paths use `gs://` URIs. Reference data is in `gs://gcp-public-data--broad-references/hg38/v0/` and `gs://gatk-sv-resources-public/`.
- PRs are squash-and-merged to `main`.
- When adding a new WDL task, follow the existing RuntimeAttr override pattern in nearby tasks (see `Utils.wdl` for the `RuntimeAttr` struct definition).
- When adding Python code to `sv-pipeline/scripts/`, it runs as a standalone script inside Docker—no package imports needed. For reusable utilities, add to `svtk` or `sv_utils`.
