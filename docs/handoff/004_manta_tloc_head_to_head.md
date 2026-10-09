# 004 — Production head-to-head for the single-tloc flag (2026-10-09)

Follows 003 (rebase onto main). Records the validation campaign for PR #968
(`mw_manta_tloc_autoresolve_pr`), which until today had only ever been tested locally.

## Question

`test_single_tloc.py` (26 cases) and one local A/B were the entire evidence base, and three of
its assertions skip locally because local pysam is 0.24.1 while the image pins **0.15.4**
(`dockerfiles/sv-pipeline-virtual-env/Dockerfile:34-40`). The END/END2 write could therefore not
be observed at all in the environment that ships. CI cannot close this: the only Python job is
`tox -e lint` (flake8) and `tox.ini` excludes both changed files, so adding tests to the PR would
not make CI execute them either.

## Leg

Step **04** (`GatherBatchEvidence` → `TinyResolve` → `ResolveManta`) is the only leg that runs this
code. `mantatloccheck.sh:30` is the line that calls `svtk resolve ... --resolve-single-tlocs`, and
`main`'s copy of the same script has no flag. Step 06 (`GenerateBatchMetrics`) does not reference
`TinyResolve` at all — an earlier plan froze 21.55 GiB of step-06 inputs before this was caught, and
no compute was spent on it. The testkit's rerun tool models steps 06-10 only, so the leg ran through
the validation harness WDL `wdl/TlocResolveOnly.wdl` (registered on Dockstore for the
`mw_manta_tloc_autoresolve` branch only; marked "VALIDATION HARNESS ONLY" in `.github/.dockstore.yml`).

## Method (reproducible)

- Images built on a throwaway GCE VM (no local Docker), same Dockerfile both sides:
  `sv-pipeline:main-743dd9d4` (control, `origin/main` at rebase time) and
  `sv-pipeline:mw-manta-tloc-autoresolve-pr-a6d61b` (PR head). Control is required: PR #966 put the
  paired-path tloc fix on main nine days earlier, so any older baseline would mix two variables.
- Inputs: September's frozen harness inputs (317 objects, 8.65 GiB) copied into the submitting
  workspace's own bucket — Cromwell reads inputs with that workspace's pet SA, so referencing another
  workspace's bucket fails even when `gsutil stat` succeeds.
- Two method configs, identical except `TlocResolveOnly.sv_pipeline_docker` (asserted programmatically
  before POST, and read back after). Root entity `sample_set/all_samples` is an anchor only.
  `useCallCache=False` stated on the submission: a run served from cache did not run the code it came
  to measure, and its output still passes.
- Arms: PR `725888ac…` (workflow `ee471849…`), control `038d1696…` (workflow `98d80456…`), both
  `Succeeded`, 54 min, $0.33.

## Result

| | control `main-743dd9d4` | PR `a6d61b` |
|---|---|---|
| records | 961,550 | 988,817 (+27,267) |
| `<CTX>` | 7 | 27,274 |
| normalised records present in control, absent from PR | — | **0** |
| records added | — | 10,519 distinct, all `<CTX>` |
| non-CTX records whose ID differs | — | 2,814 / 6,399 (44.0%) |
| CTX per sample | ≤1 | min 142, median 172, mean 174.8, max 222 |

CTX END in the shipped image (27,274 records): `END < POS` 19,203 (70.4%), `END > POS` 8,071 (29.6%),
`END == POS` **0**, `END == END2` 100%. The local 0.24.1 proxy reads these as `END == POS` because
newer pysam clamps — the proxy was wrong in the direction it was suspected of being wrong, and the
comment that quoted it was corrected in `25b8a98b` (comment-only; code identical with comments
stripped). Local and production otherwise agree: ID churn 44.8% → 44.0%, CTX/sample 172 → median 172.

Not shown by this A/B and not claimable from it: the strand/arm heuristic's cost. It is a
flag-on/flag-off effect on the same code (146 of 318 candidates on one sample locally), invisible in
a main-vs-PR diff because those records never resolved in `main` either. A production flag-off arm
would mean editing the script inside the image.

## Gotchas worth not re-learning

- `gcloud compute instances set-metadata --items=k=v` **fails on base64 `=` padding**, printing a help
  hint instead of applying. Use the REST `setMetadata` (needs `fingerprint`, is **asynchronous** —
  re-read before believing it).
- `gcloud compute instances get-serial-port-output --start N` is a **byte offset into retained
  history**, not "from the beginning"; omit it for the current window. Serial output is lost on delete,
  so save it first (`…/testkit/logs/gsv-control-build-serial.txt`).
- The builder writes its spec/markers back to instance metadata on completion, clobbering a metadata
  edit made while it runs. This builder leaves **no** `GATK_SV_BUILD_*` metadata keys — provenance lives
  in serial.
- Entity TSV import dialect is `entity:<type>_id<TAB>attr…`; the `#type` export form is rejected as
  "Unknown firecloud model entity type".
- Rawls here answers `[]` for method-config **listings** in both the REST and fiss clients while
  reading a config **by name** works. `POST /methodconfigs` requires `methodConfigVersion`.
  `/api/workspaces/v1/...` and `workflow-runner/*` are 405 (runner API not enabled) — a workflow must
  reach Terra via Dockstore + method config.
- Two arms writing to the same output attribute (`this.tloc_vcf`) **overwrite each other**; recover
  files from `…/submissions/<sub>/TlocResolveOnly/<wf>/call-ResolveManta/shard-N/attempt-M/…` and pick
  the highest attempt per shard (preemptibles retried shard-1 in one arm, shard-2 in the other).
- `jarprobe-6d795a` in broad-dsde-methods is **not ours** (parallel workstream). My builder
  `gsv-mw-tloc-pr-300g` was deleted after both images pushed; no `gsv-*` VMs or disks remain.

## Still open

- PR #968 `blocked` on human review; PR #970 too. No reviewers requested.
- PR body still carries the pre-correction claims; `/tmp/tloc_fix/pr_body_addendum.md` is drafted,
  not posted — changing claims on an open review is the owner's call.
- `test_single_tloc.py` lives on the working branch only, per the decision recorded in 002 §1a.
- Cleanup after merge: retire `wdl/TlocResolveOnly.wdl` + unpublish its Dockstore version, delete the
  two `04t-*` configs, the `tloc_frozen/` copy (8.65 GiB) and the harness WDL's `.dockstore.yml` entry,
  drop the safety refs `tmp/pre-{rebase,fix}-2026-10-0*`, remove `wt/gitignore-wt` once #970 merges.
