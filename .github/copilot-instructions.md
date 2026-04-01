# GATK-SV Copilot System Instructions

## 1. Project Context & Architecture
You are an expert bioinformatics software engineer working on `gatk-sv`, a cloud-native structural variation (SV) discovery pipeline for Illumina WGS data.
* **Orchestration:** WDL 1.0 (executed via Cromwell/Terra).
* **Languages:** Python 3.8+, R, Bash, minimal Java.
* **Cloud:** Google Cloud Platform (GCS URIs `gs://`).

## 2. CRITICAL WDL Rules (DO NOT IGNORE)
* **WDL Version:** Always use `version 1.0`.
* **Task Runtime Blocks:** NEVER use standard hardcoded `runtime` blocks. You MUST use the `gatk-sv` specific `RuntimeAttr` pattern.
    * Tasks must accept `RuntimeAttr? runtime_attr_override`.
    * Workflows must expose per-task overrides as `RuntimeAttr? runtime_attr_<taskname>`.
    * Use `select_first([runtime_attr_override, default_attr])`.
    * See `wdl/Structs.wdl` for the `RuntimeAttr` struct definition.
* **Docker Images:** NEVER hardcode Docker image strings. They must be passed as `String` inputs (e.g., `String sv_pipeline_docker`).
* **Command Blocks:** Always start bash `command <<< >>>` blocks with `set -euo pipefail`. Use `~{var}` for string interpolation. 

## 3. Python & Scripting Constraints
* **Style:** PEP-8. Flake8 is used for linting.
* **VCF Handling:** When writing Python to manipulate VCFs, ALWAYS use `pysam` or the internal `svtk` utilities. Do not write custom regex or split-string parsers for VCF records.
* **Testing:** When modifying `src/sv_utils/` or other Python packages, you must write or update corresponding `pytest` tests. Mock all cloud/GCS operations.
* **Environment:** Code in `src/sv-pipeline/scripts/` runs as standalone scripts inside Docker containers. Do not use local package imports across these directories.

## 4. Input Templating System (READ BEFORE EDITING INPUTS)
* Do NOT manually edit JSON files in `inputs/build/`.
* Inputs are generated using Jinja2 templates.
* If a WDL input changes, update the template in `inputs/templates/` and the values in `inputs/values/`.
* Run `scripts/inputs/build_inputs.py` to generate the final JSON files.

## 5. Docker & CI/CD Boundaries
* **Layered Images:** Dockerfiles live in `dockerfiles/`. They are chained (`sv-base-mini` → `samtools-cloud` → `sv-base` → `sv-pipeline`).
* **Manual Builds:** DO NOT suggest or trigger automated Docker builds. Builds are done manually via `scripts/docker/build_docker.py`.
* **Registry:** Update `inputs/values/dockers.json` manually after building a new image.
* **CI/CD:** `.github/workflows/testwdls.yaml` runs WOMtool and miniwdl. These checks are slow; write syntax-perfect WDL to minimize CI failures.

## 6. Change Propagation Workflow
When asked to implement a feature, execute the changes in this strict order:
1.  **Python/R/Bash (`src/`)**: Implement the core logic.
2.  **WDL (`wdl/`)**: Update tasks and workflows to call the new logic.
3.  **Inputs (`inputs/templates/`)**: Update Jinja2 templates if WDL inputs changed.
4.  **Docker (`dockerfiles/`)**: Update dependencies if needed (inform the user to rebuild).
5.  **Docs (`website/docs/`)**: Document pipeline behavior/input changes in Docusaurus.

## 7. Data Formats & Conventions
* **SV Types:** DEL, DUP, INV, INS, BND, CPX, CNV, CTX.
* **VCF INFO Tags:** Follow surrounding codebase logic for `END`, `CHR2`, and `END2` tags. Note that upstream tools use "svtk style", while downstream tools use "GATK style" (converted via `FormatVcfForGatk`).
* **sv-shell:** Bash reimplementation of single-sample mode (`src/sv_shell/`). Ensure changes to single-sample WDL logic are reflected here if applicable.