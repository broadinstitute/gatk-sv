#!/usr/bin/env python3

"""
Terra pipeline test runner for GATK-SV.

Sets up a Terra workspace for a feature branch, configures workflows from
Dockstore, uploads input/output configurations, and runs the pipeline
sequentially, polling for completion between modules.

Requirements:
    pip install firecloud

Usage:
    python scripts/test/run_terra_workflows.py \
        --branch my-feature-branch \
        --credentials creds.json \
        --start-module 1 \
        [--poll-interval 30] \
        [--repo-root /path/to/gatk-sv]

Credentials JSON format:
    {
        "terra_namespace": "my-billing-project",
        "terra_workspace": "my-workspace-name",
        "method_namespace": "gatk-sv",
        "entity_name": "ref_panel_1kg",
        "entity_type": "sample_set"
    }
"""

import argparse
import json
import logging
import os
import re
import subprocess
import sys
import time
import copy

import firecloud.api as fapi

logging.basicConfig(
    format="%(asctime)s %(levelname)s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
)
log = logging.getLogger(__name__)

# Ordered list of pipeline modules.  The number prefix matches the Terra
# workspace convention (e.g. "1-GatherSampleEvidence").  Each entry is
# (step_number, workflow_name, entity_type).
PIPELINE_MODULES = [
    (1, "GatherSampleEvidence", "sample"),
    (2, "EvidenceQC", "sample_set"),
    (3, "TrainGCNV", "sample_set"),
    (4, "GatherBatchEvidence", "sample_set"),
    (5, "ClusterBatch", "sample_set"),
    (6, "GenerateBatchMetrics", "sample_set"),
    (7, "FilterBatchSites", "sample_set"),
    (8, "PlotSVCountsPerSample", "sample_set"),
    (9, "FilterBatchSamples", "sample_set"),
    (10, "FilterOutlierSamples", "sample_set"),
    (11, "MergeBatchSites", "sample_set_set"),
    (12, "GenotypeBatch", "sample_set"),
    (13, "RegenotypeCNVs", "sample_set_set"),
    (14, "CombineBatches", "sample_set_set"),
    (15, "ResolveComplexVariants", "sample_set_set"),
    (16, "GenotypeComplexVariants", "sample_set_set"),
    (17, "CleanVcf", "sample_set_set"),
    (18, "RefineComplexVariants", "sample_set_set"),
    (19, "JoinRawCalls", "sample_set_set"),
    (20, "SVConcordance", "sample_set_set"),
    (21, "ScoreGenotypes", "sample_set_set"),
    (22, "FilterGenotypes", "sample_set_set"),
    (23, "MainVcfQc", "sample_set_set"),
    (24, "AnnotateVcf", "sample_set_set"),
    (25, "VisualizeCnvs", "sample_set_set"),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _run(cmd, cwd=None, check=True):
    """Run a shell command, streaming output.  Crash on failure."""
    log.info("Running: %s", cmd)
    result = subprocess.run(
        cmd, shell=True, cwd=cwd,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    if result.stdout:
        for line in result.stdout.splitlines():
            log.info("  | %s", line)
    if check and result.returncode != 0:
        raise RuntimeError(
            f"Command failed (rc={result.returncode}): {cmd}\n{result.stdout}"
        )
    return result


def _check_response(response, expected_code, context=""):
    """Assert that a FireCloud API response has the expected status code."""
    if isinstance(expected_code, int):
        expected_code = [expected_code]
    if response.status_code not in expected_code:
        raise RuntimeError(
            f"Terra API error {context}: "
            f"status={response.status_code}, body={response.text}"
        )


def _terra_config_name(step_number, workflow_name):
    """Return the Terra method-configuration name with zero-padded numeric prefix."""
    return f"{step_number:02d}-{workflow_name}"


# Matches numbers like 42, 3.14, -1e-6 etc.
_NUMBER_RE = re.compile(r'^-?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$')


def _stringify_values(d):
    """Convert a dict of workflow-config values to the string format Terra
    expects.

    Terra's Rawls API requires every input/output value to be a *string* that
    is a valid Terra expression:
      - Workspace/entity references: ``workspace.attr`` or ``this.attr``
      - String literals:             ``"quoted value"`` (quotes included)
      - Number literals:             ``42`` / ``3.14``
      - Boolean literals:            ``true`` / ``false``
      - Array literals:              ``["a","b"]``

    The JSON config files produced by the build system use a Jinja-style
    ``${workspace.attr}`` / ``${this.attr}`` wrapper that must be stripped.
    Plain string literals (e.g. ``EXACT``) must be wrapped in quotes so that
    Terra's expression parser sees them as string literals.
    """
    result = {}
    for k, v in d.items():
        if v is None:
            continue  # skip optional inputs that were not set
        if isinstance(v, bool):
            result[k] = str(v).lower()  # True/False → "true"/"false"
        elif isinstance(v, (int, float)):
            result[k] = str(v)
        elif isinstance(v, (list, dict)):
            # Compact JSON: Terra accepts array/object literals directly.
            result[k] = json.dumps(v, separators=(',', ':'))
        else:  # str
            # 1. Strip ${...} Jinja wrapper → bare Terra expression
            m = re.match(r'^\$\{(.+)\}$', v.strip())
            bare = m.group(1) if m else v

            # 2. Classify and format the bare value
            if re.match(r'^(this|workspace)\.', bare):
                # Entity or workspace attribute reference – keep as-is
                result[k] = bare
            elif bare in ('true', 'false'):
                result[k] = bare
            elif _NUMBER_RE.match(bare):
                result[k] = bare
            elif (bare.startswith('"') and bare.endswith('"')) or \
                 bare.startswith('[') or bare.startswith('{'):
                # Already a quoted string literal or JSON array/object – keep
                result[k] = bare
            else:
                # Plain string literal (e.g. "EXACT", "gs://…") – wrap in
                # double quotes so Terra sees it as a string expression.
                result[k] = f'"{bare}"'
    return result


def _load_credentials(creds_path):
    """Load and validate the credentials/config JSON."""
    with open(creds_path) as f:
        creds = json.load(f)
    required = ["terra_namespace", "terra_workspace", "sample_set_set_id"]
    missing = [k for k in required if k not in creds]
    if missing:
        raise ValueError(f"Credentials JSON missing keys: {missing}")
    creds.setdefault("method_namespace", "gatk-sv")
    if "samples_tsv" not in creds:
        raise ValueError(
            "Credentials JSON must contain 'samples_tsv' "
            "(e.g. \"samples_1kgp_156.tsv\").  "
            "Remove the old 'entity_name'/'entity_type' keys."
        )
    return creds


def _detect_batches(terra_dir, samples_tsv_name):
    """Determine which batches are covered by the given samples TSV.

    Reads the sample IDs from ``samples_tsv_name``, then cross-references
    ``sample_set_membership_1kgp.tsv`` to find every batch (sample_set) that
    has at least one sample present in the TSV.

    Returns:
        cohort_name (str)  – stem of the samples TSV filename, used as the
                             ``sample_set_set`` entity name on Terra.
        batches (list)     – sorted batch names (sample_set IDs) found.
    """
    samples_tsv_path = os.path.join(terra_dir, samples_tsv_name)
    if not os.path.isfile(samples_tsv_path):
        raise FileNotFoundError(
            f"samples_tsv '{samples_tsv_name}' not found in {terra_dir}"
        )

    sample_ids = set()
    with open(samples_tsv_path) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if parts and parts[0]:
                sample_ids.add(parts[0])

    membership_tsv = os.path.join(terra_dir, "sample_set_membership_1kgp.tsv")
    batch_samples: dict = {}
    with open(membership_tsv) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                batch, sample = parts[0], parts[1]
                batch_samples.setdefault(batch, set()).add(sample)

    batches = sorted(
        batch for batch, members in batch_samples.items()
        if members & sample_ids
    )
    if not batches:
        raise RuntimeError(
            f"No batches found in membership TSV for samples in '{samples_tsv_name}'"
        )

    cohort_name = os.path.splitext(samples_tsv_name)[0]  # e.g. "samples_1kgp_156"
    log.info(
        "Detected %d batch(es) for cohort '%s': %s",
        len(batches), cohort_name, batches,
    )
    return cohort_name, batches


# ---------------------------------------------------------------------------
# Step 1: checkout branch
# ---------------------------------------------------------------------------

def checkout_branch(repo_root, branch):
    """Ensure the repo is on the requested branch."""
    result = _run("git rev-parse --abbrev-ref HEAD", cwd=repo_root, check=True)
    current = result.stdout.strip()
    if current == branch:
        log.info("Already on branch '%s'", branch)
    else:
        log.info("Checking out branch '%s' (currently on '%s')", branch, current)
        _run(f"git checkout {branch}", cwd=repo_root)
    # Log the current SHA for traceability
    sha_result = _run("git rev-parse --short HEAD", cwd=repo_root)
    log.info("HEAD is at %s", sha_result.stdout.strip())


# ---------------------------------------------------------------------------
# Step 2: build default inputs
# ---------------------------------------------------------------------------

def build_default_inputs(repo_root):
    """Run build_default_inputs.sh to regenerate all Terra input JSONs."""
    script = os.path.join(repo_root, "scripts", "inputs", "build_default_inputs.sh")
    if not os.path.isfile(script):
        raise FileNotFoundError(f"Cannot find {script}")
    _run(f"bash {script} -d {repo_root}", cwd=repo_root)


# ---------------------------------------------------------------------------
# Step 3: update .dockstore.yml for modified workflows
# ---------------------------------------------------------------------------

def get_modified_workflows(repo_root, branch):
    """Diff current branch against main to find modified WDL workflow names.

    Returns a set of workflow names (without .wdl extension) that have been
    modified on this branch relative to main.
    """
    result = _run(
        "git diff --name-only main...HEAD -- wdl/",
        cwd=repo_root, check=True,
    )
    modified_wdls = set()
    for line in result.stdout.strip().splitlines():
        line = line.strip()
        if line.endswith(".wdl"):
            wdl_name = os.path.basename(line).replace(".wdl", "")
            modified_wdls.add(wdl_name)
    log.info("Modified WDLs: %s", sorted(modified_wdls))
    return modified_wdls


def update_dockstore_yml(repo_root, branch, modified_wdls):
    """Add the branch to .dockstore.yml for each modified workflow.

    Returns the set of workflow names that were actually updated (i.e. the
    branch was added where it wasn't already present).
    """
    dockstore_path = os.path.join(repo_root, ".github", ".dockstore.yml")
    with open(dockstore_path) as f:
        # We parse manually to preserve YAML formatting and avoid introducing
        # a ruamel/pyyaml dependency; the structure is simple and regular.
        lines = f.readlines()

    updated_workflows = set()
    # Walk through the file looking for workflow name lines, then their
    # branches lists.
    i = 0
    while i < len(lines):
        line = lines[i]
        # Detect a workflow name line
        if line.strip().startswith("name:"):
            wf_name = line.strip().split("name:")[1].strip()
            if wf_name in modified_wdls:
                # Scan forward to find the branches list
                j = i + 1
                in_branches = False
                branch_already_present = False
                last_branch_line = None
                while j < len(lines):
                    sline = lines[j].strip()
                    if sline == "branches:":
                        in_branches = True
                        j += 1
                        continue
                    if in_branches:
                        if sline.startswith("- "):
                            branch_val = sline[2:].strip()
                            if branch_val == branch:
                                branch_already_present = True
                            last_branch_line = j
                            j += 1
                            continue
                        else:
                            # End of branches list
                            break
                    # If we hit another top-level key, stop
                    if sline.startswith("- subclass:") or (
                        sline.startswith("tags:") and not in_branches
                    ):
                        break
                    j += 1

                if not branch_already_present and last_branch_line is not None:
                    # Determine indentation from the last branch entry
                    indent = ""
                    for ch in lines[last_branch_line]:
                        if ch in (" ", "\t"):
                            indent += ch
                        else:
                            break
                    new_line = f"{indent}- {branch}\n"
                    lines.insert(last_branch_line + 1, new_line)
                    updated_workflows.add(wf_name)
                    log.info("Added branch '%s' to .dockstore.yml for %s",
                             branch, wf_name)
        i += 1

    if updated_workflows:
        with open(dockstore_path, "w") as f:
            f.writelines(lines)

        # Also determine the set of ALL dockstore workflows that include
        # this branch (not just the newly added ones) since we need to
        # configure all of them in Terra.
        _run("git add .github/.dockstore.yml", cwd=repo_root)
        commit_msg = (
            f"Update .dockstore.yml: add branch '{branch}' to workflows: "
            + ", ".join(sorted(updated_workflows))
        )
        _run(f'git commit -m "{commit_msg}"', cwd=repo_root)
        log.info("Committed .dockstore.yml changes")
    else:
        log.info("No .dockstore.yml updates needed")

    return updated_workflows


def get_all_branch_workflows(repo_root, branch):
    """Return the set of workflow names in .dockstore.yml that list this branch."""
    dockstore_path = os.path.join(repo_root, ".github", ".dockstore.yml")
    with open(dockstore_path) as f:
        lines = f.readlines()

    workflows = set()
    current_wf = None
    in_branches = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("name:"):
            current_wf = stripped.split("name:")[1].strip()
            in_branches = False
        elif stripped == "branches:":
            in_branches = True
        elif in_branches and stripped.startswith("- "):
            branch_val = stripped[2:].strip()
            if branch_val == branch and current_wf:
                workflows.add(current_wf)
        elif in_branches and not stripped.startswith("- "):
            in_branches = False

    return workflows


# ---------------------------------------------------------------------------
# Step 4: update Terra workspace data (dashboard & workspace attributes)
# ---------------------------------------------------------------------------

def _upload_tsv_string(namespace, workspace, tsv_text, context):
    """Upload an in-memory TSV string as a Terra entity table."""
    import tempfile
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", delete=False
    ) as tmp:
        tmp.write(tsv_text)
        tmp_path = tmp.name
    try:
        r = fapi.upload_entities_tsv(namespace, workspace, tmp_path,
                                     model="flexible")
        _check_response(r, [200, 201], context)
    finally:
        os.unlink(tmp_path)


def update_terra_workspace(repo_root, namespace, workspace, samples_tsv_name,
                           cohort_name):
    """Upload workspace data files (dashboard markdown, workspace TSV,
    the selected samples TSV, and sample-set membership) to the Terra workspace.
    The sample_set_set entity identified by ``cohort_name`` is assumed to
    already exist in the workspace."""
    terra_dir = os.path.join(
        repo_root, "inputs", "build", "ref_panel_1kg", "terra"
    )

    # --- Update workspace description (dashboard) ---
    dashboard_path = os.path.join(terra_dir, "cohort_mode_workspace_dashboard.md")
    if os.path.isfile(dashboard_path):
        with open(dashboard_path) as f:
            dashboard_md = f.read()
        updates = [fapi._attr_set("description", dashboard_md)]
        r = fapi.update_workspace_attributes(namespace, workspace, updates)
        _check_response(r, 200, "update workspace dashboard")
        log.info("Updated workspace dashboard description")

    # --- Update workspace attributes from workspace.tsv ---
    workspace_tsv = os.path.join(terra_dir, "workspace.tsv")
    if os.path.isfile(workspace_tsv):
        with open(workspace_tsv) as f:
            header = f.readline().strip().split("\t")
            values = f.readline().strip().split("\t")
        # Strip "workspace:" prefix from header if present
        header = [h.replace("workspace:", "") for h in header]
        if len(header) != len(values):
            raise RuntimeError(
                f"workspace.tsv header/value length mismatch: "
                f"{len(header)} vs {len(values)}"
            )
        updates = [fapi._attr_set(k, v) for k, v in zip(header, values)]
        r = fapi.update_workspace_attributes(namespace, workspace, updates)
        _check_response(r, 200, "update workspace attributes")
        log.info("Updated %d workspace attributes", len(updates))

    # --- Detect batches from the chosen samples TSV ---
    _, batches = _detect_batches(terra_dir, samples_tsv_name)

    # --- Upload entity data: samples first, then membership TSVs ---
    # 1. Sample entities from the chosen TSV only.
    samples_tsv_path = os.path.join(terra_dir, samples_tsv_name)
    log.info("Uploading sample entities from: %s", samples_tsv_name)
    r = fapi.upload_entities_tsv(namespace, workspace, samples_tsv_path,
                                 model="flexible")
    _check_response(r, [200, 201], f"upload entities {samples_tsv_name}")

    # 2. Sample-set (batch) membership.
    membership_tsv_path = os.path.join(terra_dir, "sample_set_membership_1kgp.tsv")
    log.info("Uploading batch membership TSV")
    r = fapi.upload_entities_tsv(namespace, workspace, membership_tsv_path,
                                 model="flexible")
    _check_response(r, [200, 201], "upload batch membership")

    log.info("Using existing sample_set_set '%s' with batches: %s",
             cohort_name, batches)
    return batches


# ---------------------------------------------------------------------------
# Step 5: point workflows at current branch via Dockstore
# ---------------------------------------------------------------------------

def configure_workflow_branch(namespace, workspace, method_namespace,
                              branch, branch_workflows):
    """For each workflow that this branch updates, ensure the Terra method
    configuration points to the current branch's snapshot on Dockstore.

    Terra/Dockstore convention: the method source is 'dockstore' and the
    methodPath is github.com/broadinstitute/gatk-sv/<WorkflowName>.
    """
    r = fapi.list_workspace_configs(namespace, workspace, allRepos=True)
    _check_response(r, 200, "list workspace configs")
    existing_configs = {c["name"]: c for c in r.json()}

    for step_num, wf_name, entity_type in PIPELINE_MODULES:
        if wf_name not in branch_workflows:
            continue

        terra_name = _terra_config_name(step_num, wf_name)
        if terra_name not in existing_configs:
            log.warning(
                "Workflow config '%s' not found in workspace; skipping branch "
                "update (will be created during config upload)", terra_name,
            )
            continue

        cfg = existing_configs[terra_name]
        cfg_namespace = cfg.get("namespace", method_namespace)

        # Fetch the full configuration
        r = fapi.get_workspace_config(
            namespace, workspace, cfg_namespace, terra_name
        )
        _check_response(r, 200, f"get config {terra_name}")
        config_body = r.json()

        # Update the method reference to point at the correct branch
        method_ref = config_body.get("methodRepoMethod", {})
        current_version = method_ref.get("methodVersion", "")

        if current_version == branch:
            log.info("'%s' already on branch '%s'", terra_name, branch)
            continue

        log.info(
            "Changing '%s' from version '%s' to '%s'",
            terra_name, current_version, branch,
        )
        method_ref["methodVersion"] = branch
        # Remove methodUri so Rawls regenerates it from
        # methodPath + methodVersion; leaving a stale URI causes
        # Rawls to silently keep the old version.
        method_ref.pop("methodUri", None)
        config_body["methodRepoMethod"] = method_ref

        r = fapi.overwrite_workspace_config(
            namespace, workspace, cfg_namespace, terra_name, config_body
        )
        _check_response(r, 200, f"update branch for {terra_name}")

        # Verify the change persisted
        r2 = fapi.get_workspace_config(
            namespace, workspace, cfg_namespace, terra_name
        )
        _check_response(r2, 200, f"verify branch for {terra_name}")
        actual = r2.json().get("methodRepoMethod", {}).get("methodVersion")
        if actual != branch:
            raise RuntimeError(
                f"Failed to set branch on '{terra_name}': "
                f"expected '{branch}', got '{actual}'"
            )
        log.info("Verified '%s' is now on branch '%s'", terra_name, branch)


# ---------------------------------------------------------------------------
# Step 6: upload workflow input/output configurations
# ---------------------------------------------------------------------------

def upload_workflow_configs(repo_root, namespace, workspace, method_namespace,
                           branch_workflows, branch):
    """Upload (overwrite) workflow input and output configs from the built
    JSONs in inputs/build/ref_panel_1kg/terra/workflow_configurations/."""
    config_dir = os.path.join(
        repo_root, "inputs", "build", "ref_panel_1kg", "terra",
        "workflow_configurations",
    )
    output_config_dir = os.path.join(config_dir, "output_configurations")

    # First, get all existing configs so we can merge
    r = fapi.list_workspace_configs(namespace, workspace, allRepos=True)
    _check_response(r, 200, "list configs for upload")
    existing_configs = {c["name"]: c for c in r.json()}

    for step_num, wf_name, entity_type in PIPELINE_MODULES:
        if wf_name not in branch_workflows:
            continue

        terra_name = _terra_config_name(step_num, wf_name)

        # Load input configuration
        input_json_path = os.path.join(config_dir, f"{wf_name}.json")
        if not os.path.isfile(input_json_path):
            log.warning("No input config found: %s", input_json_path)
            continue

        with open(input_json_path) as f:
            inputs_raw = json.load(f)

        inputs = _stringify_values(inputs_raw)

        # Load output configuration if it exists
        outputs = {}
        output_json_path = os.path.join(
            output_config_dir, f"{wf_name}Outputs.json"
        )
        if os.path.isfile(output_json_path):
            with open(output_json_path) as f:
                outputs = _stringify_values(json.load(f))

        # Build or update the method configuration
        if terra_name in existing_configs:
            cfg_ns = existing_configs[terra_name].get(
                "namespace", method_namespace
            )
            r = fapi.get_workspace_config(
                namespace, workspace, cfg_ns, terra_name
            )
            _check_response(r, 200, f"get config {terra_name}")
            config_body = r.json()

            # Update inputs, outputs, and branch
            config_body["inputs"] = inputs
            if outputs:
                config_body["outputs"] = outputs
            config_body["rootEntityType"] = entity_type
            method_ref = config_body.setdefault("methodRepoMethod", {})
            method_ref["methodVersion"] = branch
            # Remove methodUri so Rawls regenerates it from
            # methodPath + methodVersion.
            method_ref.pop("methodUri", None)

            r = fapi.overwrite_workspace_config(
                namespace, workspace, cfg_ns, terra_name, config_body
            )
            _check_response(r, 200, f"overwrite config {terra_name}")
            log.info("Updated config '%s' inputs/outputs", terra_name)
        else:
            # Create a new method configuration.
            # The Rawls API expects the config nested under "methodConfiguration".
            method_config = {
                "namespace": method_namespace,
                "name": terra_name,
                "rootEntityType": entity_type,
                "inputs": inputs,
                "outputs": outputs,
                "prerequisites": {},
                "methodRepoMethod": {
                    "sourceRepo": "dockstore",
                    "methodPath": (
                        f"github.com/broadinstitute/gatk-sv/{wf_name}"
                    ),
                    "methodVersion": branch,
                },
                "methodConfigVersion": 1,
                "deleted": False,
            }
            config_body = {"methodConfiguration": method_config}
            r = fapi.create_workspace_config(namespace, workspace, config_body)
            _check_response(r, [200, 201], f"create config {terra_name}")
            log.info("Created new config '%s'", terra_name)


# ---------------------------------------------------------------------------
# Step 7: run the pipeline
# ---------------------------------------------------------------------------

def _submission_is_done(status):
    """Return True if a submission status indicates completion."""
    return status in ("Done", "Aborted")


def _submission_succeeded(submission_json):
    """Return True if a completed submission had no failures."""
    statuses = submission_json.get("workflowStatuses", {})
    failed = statuses.get("Failed", 0)
    aborted = statuses.get("Aborted", 0)
    return failed == 0 and aborted == 0


def submit_workflow(namespace, workspace, config_namespace, terra_name,
                    entity_name, entity_type, expression=None,
                    use_callcache=True):
    """Submit a single workflow and return the submission ID.

    ``expression`` is passed straight through to Terra's create_submission and
    is used for per-sample workflows submitted against a sample_set entity
    (e.g. expression="this.samples" scatters over all samples in the set).

    Raises RuntimeError on failure.  If the failure is due to 'Extra inputs',
    the offending keys are extracted from the error and reported clearly so
    the developer can remove them from the relevant template in
    inputs/templates/terra_workspaces/cohort_mode/workflow_configurations/.
    """
    r = fapi.create_submission(
        namespace, workspace,
        config_namespace, terra_name,
        entity=entity_name,
        etype=entity_type,
        expression=expression,
        use_callcache=use_callcache,
    )

    if r.status_code == 400 and "Extra inputs:" in r.text:
        m = re.search(r'Extra inputs:\s*([^"]+?)(?:[,.]\s*(?:Invalid|$)|\.?")', r.text)
        extra_keys = m.group(1).strip() if m else "(see error body above)"
        tmpl = (
            "inputs/templates/terra_workspaces/cohort_mode/"
            f"workflow_configurations/{terra_name.split('-', 1)[-1]}.json.tmpl"
        )
        raise RuntimeError(
            f"Terra rejected submission of '{terra_name}' due to extra inputs "
            f"that the WDL no longer declares.\n"
            f"  Extra keys: {extra_keys}\n"
            f"  Fix: remove those keys from {tmpl}\n"
            f"  Then re-run build_default_inputs.sh and retry."
        )

    _check_response(r, 201, f"submit {terra_name}")
    submission_id = r.json()["submissionId"]
    log.info("Submitted '%s' (ns=%s) → submission %s",
             terra_name, config_namespace, submission_id)
    return submission_id


def poll_submission(namespace, workspace, submission_id, poll_interval):
    """Poll a submission until it completes.  Returns the final submission
    JSON.  Crashes if the submission fails."""
    while True:
        r = fapi.get_submission(namespace, workspace, submission_id)
        _check_response(r, 200, f"poll submission {submission_id}")
        sub = r.json()
        status = sub.get("status", "Unknown")
        wf_statuses = sub.get("workflowStatuses", {})

        status_str = ", ".join(f"{k}={v}" for k, v in wf_statuses.items())
        log.info(
            "Submission %s: status=%s  workflows=[%s]",
            submission_id, status, status_str,
        )

        if _submission_is_done(status):
            return sub

        time.sleep(poll_interval)


def run_pipeline(namespace, workspace, method_namespace,
                 cohort_name, batches,
                 start_module, poll_interval, branch_workflows):
    """Run the pipeline sequentially from start_module, polling between
    modules.

    Entity routing per module type:
      ``sample``        – submitted against each batch as a sample_set with
                          expression ``this.samples`` (one sub per batch).
      ``sample_set``    – submitted once per batch.
      ``sample_set_set``– submitted once against the cohort sample_set_set.
    """
    modules_to_run = [
        (num, name, etype)
        for num, name, etype in PIPELINE_MODULES
        if num >= start_module and name in branch_workflows
    ]

    if not modules_to_run:
        log.warning(
            "No modules to run from step %d with branch workflows: %s",
            start_module, sorted(branch_workflows),
        )
        return

    log.info(
        "Will run %d modules: %s",
        len(modules_to_run),
        [f"{n}-{name}" for n, name, _ in modules_to_run],
    )

    # Fetch the actual config namespace for each workflow from the workspace.
    r = fapi.list_workspace_configs(namespace, workspace, allRepos=True)
    _check_response(r, 200, "list workspace configs for submission")
    workspace_configs = {c["name"]: c for c in r.json()}

    def _run_module(step_num, terra_name, config_namespace,
                    entity_name, entity_type, expression=None):
        """Submit + poll one module run; raise on failure."""
        log.info("=" * 70)
        label = f"{terra_name} on {entity_type}:{entity_name}"
        if expression:
            label += f" expr={expression}"
        log.info("Submitting: %s (ns=%s)", label, config_namespace)
        log.info("=" * 70)

        submission_id = submit_workflow(
            namespace, workspace, config_namespace, terra_name,
            entity_name, entity_type, expression=expression,
        )
        sub = poll_submission(namespace, workspace, submission_id, poll_interval)
        if not _submission_succeeded(sub):
            wf_statuses = sub.get("workflowStatuses", {})
            raise RuntimeError(
                f"Module {terra_name} FAILED on {entity_type}:{entity_name}.  "
                f"Submission {submission_id} statuses: {wf_statuses}"
            )
        log.info("Completed: %s on %s:%s", terra_name, entity_type, entity_name)

    for step_num, wf_name, etype in modules_to_run:
        terra_name = _terra_config_name(step_num, wf_name)

        if terra_name not in workspace_configs:
            raise RuntimeError(
                f"Cannot submit '{terra_name}': config not found in workspace "
                f"{namespace}/{workspace}"
            )
        config_namespace = workspace_configs[terra_name]["namespace"]

        if etype == "sample_set_set":
            # Cohort-level module: one submission against the sample_set_set.
            _run_module(step_num, terra_name, config_namespace,
                        cohort_name, "sample_set_set")

        elif etype == "sample_set":
            # Batch-level module: one submission per batch.
            for batch in batches:
                _run_module(step_num, terra_name, config_namespace,
                            batch, "sample_set")

        elif etype == "sample":
            # Per-sample module: scatter over samples in each batch via
            # expression so we get one workflow per sample without creating
            # hundreds of separate Terra submissions.
            for batch in batches:
                _run_module(step_num, terra_name, config_namespace,
                            batch, "sample_set",
                            expression="this.samples")
        else:
            raise RuntimeError(
                f"Unknown entity type '{etype}' for module {terra_name}"
            )

    log.info("=" * 70)
    log.info("Pipeline run complete!")
    log.info("=" * 70)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Run GATK-SV pipeline on Terra for a feature branch.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--branch", required=True,
        help="Git branch to test.",
    )
    p.add_argument(
        "--credentials", required=True,
        help="Path to JSON file with Terra workspace credentials/config.",
    )
    p.add_argument(
        "--start-module", type=int, required=True,
        help="Module number to start the pipeline from (1-based).",
    )
    p.add_argument(
        "--poll-interval", type=int, default=300,
        help="Seconds between status polls (default: 300).",
    )
    p.add_argument(
        "--repo-root", default=None,
        help="Path to the gatk-sv repo root. Auto-detected if not given.",
    )
    p.add_argument(
        "--skip-build-inputs", action="store_true",
        help="Skip running build_default_inputs.sh.",
    )
    p.add_argument(
        "--skip-dockstore-update", action="store_true",
        help="Skip updating .dockstore.yml.",
    )
    p.add_argument(
        "--skip-workspace-update", action="store_true",
        help="Skip uploading workspace data to Terra.",
    )
    p.add_argument(
        "--skip-config-upload", action="store_true",
        help="Skip uploading workflow configurations.",
    )
    p.add_argument(
        "--skip-run", action="store_true",
        help="Set up everything but do not submit workflows.",
    )
    p.add_argument(
        "--run-all-modules", action="store_true",
        help="Run all modules from start-module, not just branch-modified ones.",
    )
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    # Determine repo root
    if args.repo_root:
        repo_root = os.path.abspath(args.repo_root)
    else:
        # Auto-detect: this script lives at scripts/test/run_terra_workflows.py
        repo_root = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "..")
        )
    if not os.path.isdir(os.path.join(repo_root, "wdl")):
        raise RuntimeError(f"Cannot find wdl/ in repo root: {repo_root}")
    log.info("Repo root: %s", repo_root)

    branch = args.branch
    creds = _load_credentials(args.credentials)
    namespace = creds["terra_namespace"]
    workspace = creds["terra_workspace"]
    method_namespace = creds["method_namespace"]
    samples_tsv = creds["samples_tsv"]
    cohort_name = creds["sample_set_set_id"]

    # Detect batches from the chosen samples TSV up-front so any path/content
    # errors fail before we touch the remote workspace.
    terra_dir = os.path.join(
        repo_root, "inputs", "build", "ref_panel_1kg", "terra"
    )
    _, batches = _detect_batches(terra_dir, samples_tsv)
    log.info("Cohort: %s  Batches: %s", cohort_name, batches)

    # Step 1: checkout branch
    checkout_branch(repo_root, branch)

    # Step 2: build default inputs
    if not args.skip_build_inputs:
        build_default_inputs(repo_root)
    else:
        log.info("Skipping build_default_inputs.sh")

    # Step 3: determine modified workflows and update .dockstore.yml
    modified_wdls = get_modified_workflows(repo_root, branch)
    if not args.skip_dockstore_update:
        update_dockstore_yml(repo_root, branch, modified_wdls)
    else:
        log.info("Skipping .dockstore.yml update")

    # Determine the full set of workflows this branch covers in Dockstore
    branch_workflows = get_all_branch_workflows(repo_root, branch)
    log.info(
        "Workflows in .dockstore.yml for branch '%s': %s",
        branch, sorted(branch_workflows),
    )

    if args.run_all_modules:
        all_wf_names = {name for _, name, _ in PIPELINE_MODULES}
        branch_workflows = all_wf_names
        log.info("--run-all-modules: will run all pipeline modules")

    if not branch_workflows:
        raise RuntimeError(
            f"No workflows in .dockstore.yml reference branch '{branch}'. "
            f"Modified WDLs: {sorted(modified_wdls)}"
        )

    # Step 4: update Terra workspace data
    if not args.skip_workspace_update:
        batches = update_terra_workspace(
            repo_root, namespace, workspace, samples_tsv, cohort_name
        )
    else:
        log.info("Skipping Terra workspace data update")

    # Step 5: configure workflows to use the current branch
    configure_workflow_branch(
        namespace, workspace, method_namespace, branch, branch_workflows
    )

    # Step 6: upload workflow input/output configurations
    if not args.skip_config_upload:
        upload_workflow_configs(
            repo_root, namespace, workspace, method_namespace,
            branch_workflows, branch
        )
    else:
        log.info("Skipping workflow config upload")

    # Step 7: run the pipeline
    if not args.skip_run:
        run_pipeline(
            namespace, workspace, method_namespace,
            cohort_name, batches,
            args.start_module, args.poll_interval, branch_workflows,
        )
    else:
        log.info("Skipping pipeline run (--skip-run)")

    log.info("Done.")


if __name__ == "__main__":
    main()
