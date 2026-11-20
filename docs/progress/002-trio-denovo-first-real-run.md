# 002 — Trio mode ran end-to-end on real data and produced MOI; 3 blocking defects fixed, 1 defect open; the `Failed` verdict is a stale output name in the Terra config

**Session date:** system clock `2026-09-25 12:29 EDT` at handoff (repo declares no date source).
Git author dates for this session's commits: `b0756590` 2026-09-22, `6b666685` 2026-09-23,
`4ba78d4b` 2026-09-24. Terra's own timestamps for the five submissions run this session span
`2026-09-22T22:20Z` … `2026-09-25T15:40Z`. Both clocks are recorded because they disagree with
each other and with doc 001 (`2026-09-21`); neither was corrected.

Predecessor doc: `001-trio-denovo-single-sample.md` (implementation + two review rounds, "never
run"). This session is the first real-data execution.

---

## 0. Rebase 2026-10-08: every SHA below moved

Everything after this section was written before `main` moved. On 2026-10-08 the branch was rebased
onto `origin/main` `743dd9d4` (`a55c498f` plus the bot's image-list update) and force-pushed with
`--force-with-lease`: head `34dd0010` -> 10 commits on `main`, in two steps (the second one
removed the two Dockstore commits, see the readonly rule below). Pre-rebase history is kept at
`tmp/pre-rebase-34dd0010`, so the SHAs quoted below still resolve and none was rewritten in place. A
SHA in this doc is a fact about a Terra run, not a pointer to current code, and §13 already says to
re-verify rather than trust it.

| pre-rebase | new | commit |
|---|---|---|
| `4419315c` | dropped | publish `mw_fix_single_sample_blocking` on Dockstore |
| `b495befe` | `25d18997` | add trio de novo mode |
| `056e88da` | `c809468d` | adversarial review fixes |
| `df2ef88c` | `a8069c97` | subagent review fixes |
| `5b1c0612` | `286f3044` | doc 001 handoff |
| `e138b169` | `6583c570` | remove `validated_samples.list` |
| `91ea18cf` | dropped | Dockstore branch filter |
| `b0756590` | `859a6539` | trio list inside `ValidateTrioInputs` |
| `6b666685` | `29eff125` | stop localizing parents' CRAMs |
| `c7161978` | `08eb2a00` | this doc |
| `34dd0010` | `07b017fe` | the §13 correction commit |
| (this section) | this commit | the §0 record |

Six commits left the branch:

| dropped | why |
|---|---|
| `144a3fae` ConcatBaf | all three hunks already on `main` via `a55c498f`: the `File reference_dict` input, the call-site pass, and `java -jar ${GATK_JAR}` with `--sequence-dictionary` |
| `9f2ebe25` workspace table | patch already upstream |
| `4ba78d4b` CondenseReadCounts SIGPIPE | `main` fixed the same defect; see §14 #1 and #2 |
| `8c70779b` drop `pkg_resources` | owner decision to keep `main`'s `pkg_resources` import and the `setuptools<81` cap. The work (5 files + the Dockerfile import gate) is kept at `tmp/pkg-resources-drop-b9139b5e`, not deleted |
| `4419315c` publish `mw_fix_single_sample_blocking` on Dockstore | readonly file, see below. That branch was deleted from `origin` on 2026-10-08 anyway |
| `91ea18cf` Dockstore branch filter | readonly file, see below |

**`.github/.dockstore.yml` and `inputs/values/dockers.json` are readonly for PRs**, enforced by
`.github/workflows/readonly_check.yaml` (`##[error]Readonly file modified: ...`, per #922; an earlier
workstream hit it as `857419a0`). It red on the PR at `cf826e01` even though that commit was a
*revert* moving the file back toward `main`, because the check names any file in the diff rather than
asking whether the net change is empty. Two consequences for the next session: a restore commit does
not fix it, so the commits touching those paths have to come out of history; and with `91ea18cf`
gone, this branch no longer publishes its single-sample WDL to Dockstore, so Terra testing driven by
a Dockstore-published descriptor has to come from `main` after merge or from a maintainer edit.

Two claims below are wrong as written now. The headline's "3 blocking defects fixed" counts the
CondenseReadCounts SIGPIPE fix, which this branch no longer carries, so 2 of the 3 are this branch's.
And the `wdl/CollectCoverage.wdl` row in §13 is no longer a deliverable of this branch. Doc 001's
SHAs (`8995b2d0`, `1763a774`, `a51ae8d6` on base `857419a0`) were already orphaned before this
rebase, from the branch's earlier base; left as written.

Static checks were re-run against the rebased tree (gatk-sv-testkit `0a56383` +dirty, measured
`main`@`743dd9d4` vs `trio_denovo_single_sample`@`c49e2b9e`): no new WDL launchability findings, the
`sv_shell` contract findings are identical to `main` (the same 7 unsupplied reads), 0 new jq nulls
against `main`, `annotate_moi` tests 4/4 passed, and `wdl_reach` confirms `AnnotateModeOfInheritance`
is reached from the single-sample workflow at `GATKSVPipelineSingleSample.wdl:1798` and reaches
`annotate_moi.py`. Two gaps, both honest: the womtool required-input half is skipped on *both* refs
(no `WOMTOOL_JAR`, local java is 11), and `svtk` was never imported locally (no local install; the
cover is byte-identity with `main`).

---

## 1. Headline: trio MOI annotation works

Full-genome trio run (`ffe3c189-5997-4c66-8de2-ba3ec6d77d03`, workflow `d4456607-ad4e-4c13-8de2-cc877a104681`)
completed every task. Cromwell's own root status is **`Succeeded`**. Terra marks the workflow
`Failed` for one reason only, recorded in §7.

`SM-GN4BI.moi.moi_summary.tsv` (135 bytes, verified by `gsutil cat` of the workflow output):

| MOI | COUNT |
|---|---|
| DE_NOVO | 749 |
| INHERITED_FROM_MOTHER | 2506 |
| INHERITED_FROM_FATHER | 2542 |
| INHERITED_FROM_BOTH | 6713 |
| PARENT_ONLY | 4365 |
| UNASSESSABLE | 914 |

Independently re-counted from the final VCF itself (not from the summary): 17 789 variant records,
`MOI` distribution identical to the table above (6713 / 4365 / 2542 / 2506 / 914 / 749).
`MOI_CONFIDENCE`: `CONFIRMED` 16758, `UNCONFIRMED` 1031. Verified by
`python3` over `gzip.open()` of the 10 972 534-byte `final_vcf`.

Header facts from the same file:

- sample columns, in order: `SM-GN4BI` (case), `SM-GZQLZ` (mother), `SM-GZQMC` (father)
- `##INFO=<ID=MOI,Number=1,Type=String,Description="Mode of inheritance relative to the case sample: DE_NOVO|INHERITED_FROM…">`

Other workflow outputs collected by Cromwell (17 total): `final_vcf` + `_idx`, `pre_cleanup_vcf` +
`_idx`, `metrics_file`, `qc_file`, `ploidy_matrix`, `ploidy_plots`, `working_ped`,
`non_genotyped_unique_depth_calls` + `_idx`, `stripy_{html,json,tsv,vcf}_output`,
`gd_output_tarball = None`. That null is **correct**: the GD branch is skipped in trio mode.

Durable copies of five artifacts (workspace scratch buckets are deleted with the workspace):

```
/Users/markw/Work/genotypebatch_debug/testkit/staging/trio-results-ffe3c189/
  SM-GN4BI.moi.moi_summary.tsv        135 bytes  sha256 aeb15729fbd8…
  SM-GN4BI.moi.vcf.gz              10972534 bytes  sha256 f76cf2fa6fc1…
  single_sample.SM-GN4BI.metrics.tsv  16102 bytes  sha256 7789430a446b…
  sv_qc.SM-GN4BI.tsv                   1450 bytes  sha256 f1c149213024…
  trio_combined_ped_file.ped           3601 bytes  sha256 902a894ae5e1…
```

`sha256sum` prefixes recorded above are 12 chars; recompute locally with
`cd …/trio-results-ffe3c189 && sha256sum *`.

### QC verdict, recorded as observed, not explained away

`sv_qc.SM-GN4BI.tsv`: **10 PASS / 17 FAIL**. The failures are all `final_vcf_*` count/size-band
assertions (e.g. `final_vcf_DEL_pass_count => FAIL`, `final_vcf_DUP_pass_size_500_5000 => PASS`).
Thresholds in force: workspace attribute
`ref_panel_qc_definitions = gs://gatk-sv-resources-public/hg38/v0/sv-resources/ref-panel/1KG/v2/single_sample.qc_def`
— the shipped **single-sample** QC definitions. *Inferred, not proven:* those bands were calibrated
for a 1-sample callset and a trio callset has different composition, so most FAILs are expected
rather than regressions. **Not verified either way** — nobody has diffed this run's QCDF against a
v1.1.1 single-sample baseline QCDF. That is open item §14 #5.

`SingleSampleMetrics` contains **0 keys** matching `moi|trio|inherit` (413 metric rows total,
verified by `cut -f1 | grep -icE`). MOI coverage exists only in `moi_summary` + VCF `INFO`, so
nothing in metrics would catch MOI regressions. See open item §14 #4.

---

## 2. The five submissions, with root cause per failure

| # | submission | Terra date | outcome | cost | root cause |
|---|---|---|---|---|---|
| 1 | `6da799b0-3805-4158-b868-94d026d4566a` | 09-22 22:20Z | Failed ~10 s | $0.00 | defect A: workflow-scope `write_lines` |
| 2 | `79b10bb0-64f9-410e-b176-95687131ebc5` | 09-22 22:39Z | Failed | $6.47 | defect B: `ValidateTrioInputs` localized both parental CRAMs |
| 3 | `57b1b186-6d27-4d9c-afdc-06a9d1f46026` | 09-23 15:12Z | Failed 11 h 47 m | $1.16 | defect C: SIGPIPE rc=141 in `CondenseReadCounts` |
| 4 | `ffe3c189-5997-4c66-8de2-ba3ec6d77d03` | 09-24 15:29Z | **all work Succeeded**; Terra `Failed` | $18.87 | defect E: stale `final_bed` in config outputs |
| 5 | `3d7f4e75-4eba-4d70-9f93-f7f8454f2ebe` | 09-25 15:40Z | **`Running`, cost ~$0.00 at handoff** | ~$0.00 | — (cache-hit rerun of #4 with the outputs map fixed) |

Total spend at handoff: **$26.50** (owner-approved envelope was ~$15–25 *per run*; run 4 alone
exceeded it). No build VMs and no `gsv-*` compute instances existed at handoff
(`gcloud compute instances list --project broad-dsde-methods --filter="name~gsv-"` returned zero
rows; the build VM from §3 was deleted after reading its serial log).

Preemption context for run 3 and 4 (measured from `?expandSubWorkflows=true` metadata):
`GCP Batch task exited with VMPreemption(50001)` hit **hundreds** of task attempts
(`CNMOPS.NormalR1/R2`, `SVCluster`, `CalcMedCov`, `CreatePloidyTableFromPed`, `ConcatBafCase`,
`MergeSREvidence`, `SDtoBAF`, `ResolveManta`, `CondenseReadCounts:1`). All recovered via
`preemptible_tries: 3` except the ones named in defect C. Run 3's 12 h at $1.16 was **spot capacity
queueing, not cheapness** — run 4 had real capacity and spent $8 in its first 50 minutes.

---

## 3. Image build (sv-pipeline only)

`./docker/gatk-sv-build.sh trio_denovo_single_sample sv-pipeline` from the **testkit** checkout
(that script lives in `gatk-sv-testkit`, not in `gatk-sv`); VM `gsv-trio-denovo-single-sample-91ea18cf`,
e2-standard-8, 150 GB pd-ssd, `us-central1-a`; result `### GATK_SV_BUILD_RESULT=SUCCESS`, 0 errors.

```
us.gcr.io/broad-dsde-methods/markw/gatk-sv/sv-pipeline@sha256:dcd40e2a01492351b9d93c2f9c1e45e1b8d0bfdc1d2bf2497ca1f5ecc97cf980
us.gcr.io/broad-dsde-methods/markw/gatk-sv/sv-shell@sha256:ef40df98…   (dependency, same build)
```
Tag `trio-denovo-single-sample-91ea18` built from `91ea18cf`, i.e. **before** the three defect fixes
`b0756590`, `6b666685`, `4ba78d4b`. Those fixes are WDL-only, so the image stayed valid; the
workspace pins the **digest**, not the tag, so the run is byte-reproducible.

Build attempt #1 failed with `ModuleNotFoundError: No module named 'pkg_resources'`; root cause
(`pip install -e /opt/gatk-sv-gd` floats conda-pinned setuptools 62.1 → 84.x, `pkg_resources` removed
at 82.0.0) was independently fixed upstream as `8c70779b`. VM deleted after reading the serial log.

---

## 4. Dockstore: how the branch got published, and how to *prove* which commit is served

`.github/.dockstore.yml` filters `branches` — a feature branch is never indexed unless it is listed.
Commit `91ea18cf` added `- trio_denovo_single_sample` to the `SingleSamplePipeline` entry
(`name: SingleSamplePipeline`, `primaryDescriptorPath: /wdl/GATKSVPipelineSingleSample.wdl`,
`version: 1.2`).

Working verification, all of it unauthenticated:

```bash
ID='%23workflow%2Fgithub.com%2Fbroadinstitute%2Fgatk-sv%2FSingleSamplePipeline'
B=https://dockstore.org/api/ga4gh/trs/v2       # NOT /api/ga4gh/v2 — that 404s
curl -s "$B/tools/$ID/versions/trio_denovo_single_sample/WDL/descriptor" \
  | python3 -c "import json,sys,hashlib; print(hashlib.sha256(json.load(sys.stdin)['content'].encode()).hexdigest()[:16])"
git -C <worktree> cat-file -p <sha>:wdl/GATKSVPipelineSingleSample.wdl | sha256sum | cut -c1-16
```
The two hex prefixes must match. Observed re-index latency after `git push`: **~200 s** (poll
`seq 1 8` × 40 s). Measured matches: `dd972614…` = `91ea18cf`, `62331f05…` = `b0756590`,
`270997c6…` = `6b666685`.

**This proves only the primary descriptor.** `wdl/CollectCoverage.wdl` (defect C) is an *import*, so
its sha is invisible to the descriptor check, and TRS on Dockstore has **no** bundle/files endpoint
(`/WDL/uri` and `/versions/{v}/files/{path}` both `{"code":404}`). Ground truth for imports is the
**rendered script in the run's scratch bucket**:

```bash
gsutil cat '<bucket>/…/call-CondenseReadCounts/shard-0/script' | grep -nE '!found|found = 1'
```
Run 4 printed both lines ⇒ Rawls re-resolved the branch including imports. Rawls did **not** serve a
stale branch snapshot on any of the five resubmits.

---

## 5. Coordinates: Terra sandbox and config assembly

| what | value |
|---|---|
| workspace | `gsv-trio-denovo-a0e10b99` (cloned from `help-gatk/GATK-Structural-Variants-Single-Sample`, so the cohort sandbox was never touched) |
| config | `single-sample-trio-a0e10b99`, 116 inputs, 11 outputs, `methodVersion = trio_denovo_single_sample` |
| entity | `sample` / `SM-GN4BI`, expression **`this`** |
| case | demo CRAM `SM-GN4BI`; parents = set_5 `SM-GZQLZ` (mother), `SM-GZQMC` (father), CRAM + `.crai` |
| durable payload | `/Users/markw/Work/genotypebatch_debug/testkit/staging/single-sample-trio-a0e10b99.config.json` |
| helper | `/Users/markw/Work/genotypebatch_debug/testkit/staging/gsvpeek.sh` (bounded call-level peek) |

Two 400s that cost a round trip each and are **not** documented anywhere:

1. Config `inputs` keys must be **workflow-qualified** (`GATKSVPipelineSingleSample.foo`). Stripped
   keys make every input look extra (117 extra / 87 missing). Rebuilt from the branch's flat bundle
   + the official `help-gatk/gatk-sv-single-sample` @v1.1 call-scoped overrides + the 6 trio inputs.
2. Submission entity expression must be `this`, not `this.sample.SM-GN4BI`.

Also: `create_workspace_config`/`overwrite_workspace_config` act as validators and return
`extraInputs` / `invalidInputs` / `invalidOutputs` — but **they do not validate output names**
(see §7). `overwrite_workspace_config` requires the full body including `namespace`.
`sourceRepo` must be `dockstore`; `github` → 400 `Illegal method repo 'github'`.

---

## 6. Defects

### A — workflow-scope `write_lines()` cannot run on PAPIv2. FIXED `b0756590`

`wdl/GATKSVPipelineSingleSample.wdl:637` had `File raw_trio_samples_list = write_lines(trio_samples)`.
Exact Cromwell text:

> Failed to evaluate 'raw_trio_samples_list': Evaluating write_lines(trio_samples) failed: Could not
> build the path "write_lines_677072024d572b4785763bd65898e6f6.tmp". It may refer to a filesystem not
> supported by this instance of Cromwell. Supported filesystems are: DRS, Google Cloud Storage, HTTP.

It was the only workflow-scope `write_*` in any production WDL (`wdl/TestUtils.wdl` has them and is
local-CI only), and the repo's own workaround task `WriteLines` (`wdl/Utils.wdl:698`) was never
called. Fix: pass `Array[String] trio_samples` into `ValidateTrioInputs` and materialize inside the
command. Proven by the artifact `write_lines_677072024d572b4785763bd65898e6f6.tmp` in the run-2 call
directory containing exactly `SM-GN4BI / SM-GZQLZ / SM-GZQMC`.

### B — `ValidateTrioInputs` downloaded both whole-genome parental CRAMs. FIXED `6b666685`

`File? mother_cram` / `File? father_cram` existed only to test `defined()`, so Cromwell localized
them into a task asking for `local-disk 10 HDD`. From `attempt-2/gcs_localization.sh`:

```
gs://…/SM-GZQLZ/v1/NA19238_NA19238_A_SM-GZQLZ_v1.cram
gs://…/SM-GZQMC/v1/NA19239_NA19239_A_SM-GZQMC_v1.cram
```

Both attempts ran ~20 min, wrote **no stdout/stderr**, and died with
`The job was stopped before the command finished. Check GCP Batch job logs for details.` — the
command never started. Fixed by passing `has_mother_cram`/`has_father_cram` Booleans from the call
site. Note the same trap still applies to `dragen_vcf` / `case_*_vcf` on that task — only when
supplied, far smaller, and rejected in trio mode; deliberately left alone, open item §14 #3.

### C — SIGPIPE rc=141 in `CondenseReadCounts` on whole-genome runs. FIXED on this branch `4ba78d4b`; still open on main

`wdl/CollectCoverage.wdl:134`, `set -euxo pipefail` + `existing_sample_id=$(zcat ~{counts} | awk '/^@RG/ { …; exit }')`.
Failure text: `Task GatherBatchEvidence.CondenseReadCounts:0:2 failed. Job exit code 141.` (141 =
128+SIGPIPE). Measured from this run's own input: `SM-GN4BI.counts.tsv.gz` = 131 MiB with the first
`@RG` at **header line 3368** (dictionary inlined as thousands of `@SQ` lines) — awk closed the pipe
with >130 MB still pending, far beyond the 64 KB pipe buffer. Deterministic for that header, not a race.
Shards 0 and 2 exhausted retries (fatal); sibling probes at `:113` and `:164` already carry `|| true`.

Reproduced locally: old block → `+ id=SM-GN4BI` then `exit=141`; new block → same value, `exit=0`.

**DIVERGENCE TO RECONCILE.** Another agent independently fixed the same bug on
`mw_fix_single_sample_blocking` as `04fa5142` ("Guard zcat header probes against SIGPIPE under
pipefail (exit 141)"), using **`|| true`** in six sites (`CollectCoverage`, `AnnotateExternalAFPerShard`,
`RefineComplexVariants` ×2, `TasksMakeCohortVcf`, `Vapor`, `WGD`) plus `scripts/test/test_sigpipe.sh`
and a Cloud Build run (`abe1c6d0-11a2-40e9-9828-68e4559a1e67`, SUCCESS). Their commit independently
reaches the same conclusion (135 MB input, "not a flake"). Their fix is **broader** (5 sites mine
doesn't touch); mine is **stricter** (drops the early `exit` so no error is masked — the owner chose
that approach for this branch). Both are currently pushed and neither is merged. Whichever branch
merges must decide once. Open item §14 #1.

### D — `AddTrioSamplesToPed` writes an invalid pedigree. OPEN, not fixed

`wdl/GatherBatchEvidence.wdl`, task `AddTrioSamplesToPed`. Its own comment says "Parents carry `0`
for their own parent columns" but the code zeroes only the column matching the sample:

```sh
PED_FATHER="$FATHER_ID"; if [ "$sample" = "$FATHER_ID" ]; then PED_FATHER="0"; fi
PED_MOTHER="$MOTHER_ID"; if [ "$sample" = "$MOTHER_ID" ]; then PED_MOTHER="0"; fi
printf 'trio_denovo\t%s\t%s\t%s\t%s\t1\n' "$sample" "$PED_FATHER" "$PED_MOTHER" "$SEX"
```

Observed in the run's own `working_ped` (159 rows; trio rows appended to the 156-row ref-panel ped):

```
trio_denovo  SM-GN4BI  SM-GZQMC  SM-GZQLZ  2 1     (case: PAT=father, MAT=mother  — correct)
trio_denovo  SM-GZQLZ  SM-GZQMC  0         2 1     (mother's FATHER is the case's father)
trio_denovo  SM-GZQMC  0         SM-GZQLZ  1 1     (father's MOTHER is the case's mother)
```

The parents are not founders and the pedigree is a 2-cycle. `SEX` is correct in all three rows, so
ploidy/sex assignment was unaffected in this run — that is an observation about *this* run, not proof
of harmlessness for any consumer that traverses the pedigree. One-line fix: parents emit `0` for both
columns. Open item §14 #2.

### E — stale `final_bed` output name marks a successful run Failed. OPEN (config-side)

Run 4 message: `output named GATKSVPipelineSingleSample.final_bed does not exist` while Cromwell's
root status was `Succeeded`. `final_bed` is declared by **no** single-sample WDL — 0 occurrences in
this branch's `wdl/GATKSVPipelineSingleSample.wdl` **and** 0 in `origin/main`'s copy; repo-wide the
string exists only in `wdl/RunDeNovoSVs.wdl`, a different workflow. It reached the config because the
config was cloned from the published `help-gatk/gatk-sv-single-sample` @v1.1 outputs map. Rawls does
not reject unknown output names at config time (`invalidOutputs=0`), so the cost is discovered only
after the whole run. Consequence: **any** current single-sample run driven by that published config
gets a `Failed` verdict after full success. Open item §14 #6.

---

## 7. Corrections — these earlier claims were wrong and one already shaped a decision

1. **"Cannot determine which commit Dockstore serves for a branch version."** FALSE. TRS is mounted
   at `/api/ga4gh/trs/v2` (earlier probes used `/api/ga4gh/v2`, which 404s) and the descriptor
   endpoint returns the WDL inside a JSON `content` field. Method in §4, verified three times.
2. **"`main` does not have the `gatk-sv-gd` lines, so the `pkg_resources` break is not a main problem."**
   FALSE — local `main` was stale (`4bc70a69` vs `origin/main e1909d2f`). Current main **does** have
   the lines; the breakage is a main-line problem, subsequently fixed as `8c70779b`.
3. **"`a0e10b99` (my own ConcatBaf jar fix) is the fix in force."** SUPERSEDED — dropped in the rebase
   in favour of `144a3fae`, which adds the required `--sequence-dictionary`. `RunCNVNonGenotyper`
   still runs `/gatk/gatk` under `sv_base_mini_docker`; that one is **still not exercised** (GD is
   skipped in trio mode), so it remains an untested claim, not a fixed one.
4. **Standing instructions live in `docs/CLAUDE.md`.** FALSE for both repos: neither `gatk-sv`
   (main or this worktree) nor `gatk-sv-testkit` has `CLAUDE.md` or `AGENTS.md` anywhere
   (`ls` verified). Durable facts therefore went into `gatk-sv-testkit/docs/terra-head-to-head.md`
   §8 ("API edges that bite"), chosen because it is the doc you are already reading when driving a
   Terra submission — there is no standing-instructions file to hold them.
5. **"A resubmit will be near-free: every task is cached, so run 5 costs ~$0–1 in tens of minutes."**
   WRONG, and the option the owner chose was described with that number. Measured 3 h 15 m into run 5:
   33 calls `Done` at $0.00 (cache replay) and **`MakeCohortVcf` re-executing**, so everything from
   there downstream — `RefineComplexVariants`, `AnnotateModeOfInheritance`, `SingleSampleMetrics`,
   `SingleSampleQC` — pays again. *Inferred mechanism, not proven:* run 4's own log showed
   `call-MakeBincovMatrixColumns/shard-0/cacheCopy/…`, i.e. caching works when inputs are external
   (reference, demo CRAM); tasks consuming large intermediates living in each submission's own
   `gs://fc-…/submissions/<id>/…` bucket path cannot key-match across submissions, so the cache chain
   breaks at the first task whose input is a big per-submission file. Practical rule: **cache reuse
   covers evidence gathering, not the back half.** Correct any estimate that assumes otherwise.

---

## 8. Resume here

```bash
SID=3d7f4e75-4eba-4d70-9f93-f7f8454f2ebe; WF=41a320cb-0929-4372-a8cf-e621d4d28e47   # run 5
 cd /Users/markw/IdeaProjects/gatk-sv-testkit
 NS=$(./kit/gsvtk-config get GSVTK_TERRA_NAMESPACE); WS=gsv-trio-denovo-a0e10b99
 # one cheap snapshot (exit 0 terminal / 3 still running / 4 terminal-with-failure):
 python3 /Users/markw/.pi/agent/skills/terra-monitor/scripts/twatch.py status $SID -w "$NS/$WS"
 # or block until terminal (never pass --diagnose; it 405s):
 python3 /Users/markw/.pi/agent/skills/terra-monitor/scripts/twatch.py watch $SID -w "$NS/$WS" --max-wait 1500
 # call-level progress (reads metadata into /tmp, prints only a summary):
 /Users/markw/Work/genotypebatch_debug/testkit/staging/gsvpeek.sh $SID $WF $WS
 # expected outputs of a good run 5, straight from Cromwell:
 TOK=$(gcloud auth print-access-token)
 curl -s -H "Authorization: Bearer $TOK" \
   "https://api.firecloud.org/api/workspaces/$NS/$WS/submissions/$SID/workflows/$WF" -o /tmp/m.json
 python3 -c "import json;d=json.load(open('/tmp/m.json'));[print(k.split('.')[-1],v) for k,v in (d.get('outputs') or {}).items()]"
```

Branch state must stay: local head == `origin/trio_denovo_single_sample` == `4ba78d4b`, worktree
clean (`wt/trio-denovo`). Re-verify, do not trust this line:

```bash
 git -C /Users/markw/IdeaProjects/gatk-sv/wt/trio-denovo status -sb
 git -C /Users/markw/IdeaProjects/gatk-sv/wt/trio-denovo log --oneline -1
 git -C /Users/markw/IdeaProjects/gatk-sv/wt/trio-denovo log --oneline origin/trio_denovo_single_sample -1
```

## 9. What good looks like — expected values for the next check

| check | command | expected | why |
|---|---|---|---|
| run 5 verdict | `twatch.py status` above | `Succeeded`, not `Failed` | §6 E: `final_bed` removed from the config outputs map at handoff |
| run 5 cost | same | **not** ≤$1. Observed at handoff: $0.00 after 3 h, rising once `MakeCohortVcf` started. Budget single-digit $, expect the tail to pay. | call caching covered 33 calls but **not** `MakeCohortVcf` — see the correction in §7 item 5 |
| run 5 duration | same | hours, not tens of minutes | same |
| MOI registered as a Terra output | outputs listing above | `moi_summary` present **and** the config now registers it | added `moi_summary` to the outputs map so the trio result is visible in the UI |
| `moi_summary` values unchanged | `gsutil cat` the new run's `moi_summary` | DE_NOVO 749 / MOTHER 2506 / FATHER 2542 / BOTH 6713 / PARENT_ONLY 4365 / UNASSESSABLE 914 | no code change between run 4 and run 5, so identical counts prove the resubmit was a cache replay, not a new run |
| import freshness, if WDL changes again | `gsutil cat …/call-CondenseReadCounts/shard-0/script \| grep '!found'` | grep hits | §4: descriptor sha cannot prove import freshness |
| static gate after any WDL edit | `testkit/scripts/gsvtk gate trio_denovo_single_sample --repo /Users/markw/IdeaProjects/gatk-sv` | `GATKSVPipelineSingleSample IC=2 stale=0`, `SVShell IC=1 stale=0`, unchanged vs base `4419315c` | measured identical at `b0756590` **and** `6b666685`; IC/stale counts are the zero-delta baseline |
| miniwdl | `for f in wdl/*.wdl; do /tmp/wdl-venv/bin/miniwdl check "$f"; done` | 119 OK, 0 FAIL (note: `miniwdl check --quiet` is invalid) | measured 119/0 after each of the three fixes |
| MOI unit test | `/tmp/moi-venv/bin/python -m pytest src/sv-pipeline/05_annotation/scripts/test_annotate_moi.py` | 4 passed | measured this session |

## 10. Gotchas actually hit this session (each with the text that produced it)

- **Per-workflow metadata can be a stale cached snapshot.** Two consecutive fetches 30+ min apart
  returned byte-identical JSON while the scratch bucket proved the run had advanced (new
  `call-SampleFilterMetrics`, `call-SampleFilterQC` dirs). Use the bucket + submission `cost` for
  liveness; treat metadata as lagging.
- **Metadata omits sub-workflow internals by default.** Top-level `calls` stayed at 9 entries for
  hours while ~200 tasks ran inside `GatherBatchEvidence`. `?expandSubWorkflows=true` returns
  **44 838 261 bytes** — always `curl -o` to a file and parse a summary; never inline.
- `list_workspace_configs` **and** `list_submissions` return `[]` for this workspace even with 1
  config and 5 submissions present. Use `curl` on `/api/workspaces/{ns}/{ws}/submissions` instead —
  that is how a double-submit was ruled out after a `terra.submit` raised mid-call.
- **Output expressions must not be entity-qualified.** `400 … Invalid outputs:
  GATKSVPipelineSingleSample.moi_summary -> Entity references not permitted in the middle of output
  expressions`. The value is `this.moi_summary` (workflow prefix belongs to the key, not the value).
- **Config input keys must be workflow-qualified** (§5 bullet 1) and the submission expression must be
  exactly `this` (§5 bullet 2).
- `File` inputs are **downloaded** into the task's working directory before the command runs, even if
  the command only tests `defined()` — that is defect B, and it produced no logs at all, which looks
  like an image-pull failure until you read `gcs_localization.sh`.
- `141` is SIGPIPE, and it can be deterministic. `set -o pipefail` + an `awk`/`head` that exits early
  against a multi-MB `zcat` (§6 C). `bgzip`/`gzip -dc` on a truncated stream can yield **0** bytes, so
  `gsutil cat … | head -c N | gzip -dc` reads a header as empty; download the file (11 MB) instead.
- `twatch.py --diagnose` fails with HTTP **405** (~405 s of retry) — skill defect, reported, not
  patched. Cromwell VM hostnames are unreachable from here.
- The Browser tool is unusable until the vendor binary is repaired: `Managed BetterChromium is
  outdated … Run betterwright setup`, and the wrapper itself fails `env: bun: No such file or
  directory`. `bun` was deliberately **not** installed for a single lookup.
- No PyYAML in any local venv; `.dockstore.yml` checks were done with text inspection, so a YAML
  syntax error in that file would not have been caught locally.

## 11. Mutation ledger — verify / undo each entry

Reversible / verifiable, all made by this session:

| mutation | verify / undo |
|---|---|
| pushed `91ea18cf` (force, `--force-with-lease`, rebased trio onto `4419315c`) | `git ls-remote origin refs/heads/trio_denovo_single_sample` |
| pushed `b0756590`, `6b666685`, `4ba78d4b` normal fast-forwards | same; `git log origin/trio_denovo_single_sample -4` |
| branch `trio-rebase-rehearsal` + scratch worktree `/tmp/trio-rebase` (branch tip `5b1c12`… kept) | `git worktree list`; delete with `git worktree remove /tmp/trio-rebase` |
| GCR images `markw/gatk-sv/{sv-pipeline,sv-shell}` at the §3 digests, path no pipeline reads | `gcloud container images list-tags …` |
| GCE VM `gsv-trio-denovo-single-sample-91ea18cf` created, serial log read, **deleted** | `gcloud compute instances list --filter=name~gsv-` → zero rows |
| workspace `gsv-trio-denovo-a0e10b99` created (clone), 95 attributes, `sv_pipeline_docker` + `sv_pipeline_qc_docker` pinned to the sv-pipeline **digest** | `terra.workspace(ns,'gsv-trio-denovo-a0e10b99')` |
| config `single-sample-trio-a0e10b99` created, then `overwrite_workspace_config` × ~6 (last two removed `final_bed`, added `moi_summary`) | `terra.config_payload(ns,ws,ns,name)` |
| 5 workflow submissions (§2), total **$26.50** | Job Manager / the §8 REST call |
| durable files under `…/testkit/staging/` (payload, `workspace.tsv`, trio results dir, `gsvpeek.sh`) | `ls -l` that dir |
| ponytail-review of PR #971, then the 2026-10-08 review-fix round: 3 commits (`886d93c7` pedigree + parent roles + trio guard, `5bffc79a` filter dedupe + header read, `ca635e74` CNV→BND abort + MOI header + new tests) plus the `single.md` pass, each cherry-picked from a private lane worktree (`pr971fix/lane{1,2,3}-*`) onto `6499b294` | `git log --oneline 6499b294..HEAD`; `git show <sha>`; undo with `git reset --hard 6499b294` |
| pushed the review-fix round + this doc to `trio_denovo_single_sample` (plain fast-forward, no force) | `git ls-remote origin refs/heads/trio_denovo_single_sample` equals local `git rev-parse HEAD` |
| `gatk-sv-testkit/docs/terra-head-to-head.md` §8 added (docs-only, path-scoped commit — the repo's other uncommitted files were left untouched) | `git -C …/gatk-sv-testkit log -1 --stat` |

Not touched, deliberately: the shared baseline workspace
`broad-firecloud-dsde-methods/GATK-Structural-Variants-Joint-Calling`, the cohort sandbox, `origin/main`.

## 12. Deliverables

| file | change | commit |
|---|---|---|
| `wdl/GATKSVPipelineSingleSample.wdl` | trio list materialized inside `ValidateTrioInputs`; CRAM inputs → Booleans | `b0756590`, `6b666685` |
| `wdl/CollectCoverage.wdl` | not a deliverable after the 2026-10-08 rebase: `main` `a55c498f` carries a SIGPIPE guard, this branch's variant was dropped | none |
| `.github/.dockstore.yml` | `trio_denovo_single_sample` added to `SingleSamplePipeline` branch filter | `91ea18cf` |
| repo root | removed accidentally committed `validated_samples.list` | `e138b169` |
| `docs/progress/002-trio-denovo-first-real-run.md` | this doc | this commit |
| `gatk-sv-testkit/docs/terra-head-to-head.md` | new §8: Terra/Dockstore/Cromwell API facts with the exact failure text for each | separate testkit commit |

## 13. State at handoff, re-queried not carried forward

- Run 5 `3d7f4e75-4eba-4d70-9f93-f7f8454f2ebe` / wf `41a320cb-0929-4372-a8cf-e621d4d28e47`, submitted
  `2026-09-25T15:40:54Z`: **`Running`**, 33 calls `Done` + `MakeCohortVcf` `Running`, cost **$0.00**,
  ~3 h 30 m elapsed at handoff. Submission-level status reads `Submitted` while its workflow reads
  `Running` — that lag is normal here, trust the workflow.
  **Liveness confirmed on the bucket, not on metadata** (which had served the same 641 532-byte
  snapshot twice, 10 min apart): `call-MakeCohortVcf/` holds 1 object (`script`) and **no `stdout`**,
  i.e. the task is executing now. Expect the back half to run for real — §7 item 5.
  **No outcome is assumed.** Re-run the §8 command; expect §9's corrected cost/duration rows.
- Workspace submissions total **5**, all mine (`6da799b0`, `79b10bb0`, `57b1b186`, `ffe3c189` Done;
  `3d7f4e75` active). Re-listed with `curl` at handoff precisely because `list_submissions` returns
  `[]` (§10) and because I needed to rule out a parallel agent submitting into this sandbox — there is
  none. Note each finished submission stores the config under a per-submission snapshot name
  (`single-sample-trio-a0e10b99_B0EJFlC5SLk`), so do not mistake those for extra configs.
- `wt/trio-denovo`: clean. Code head `4ba78d4b` (last WDL change); docs commits on top — `c7161978`
  (this doc) then the correction commit that records this line. Local head == `origin/trio_denovo_single_sample`
  at every check made this session; re-verify per §8 rather than trusting the sha here.
- `wt/fix-ss-blocking`: clean, head `04fa5142`, == remote (another agent's branch, untouched here).
- `wt/jrc_args`: **ahead=2, nothing pushed**; `wt/tloc-pr`: **no upstream, never pushed** — neither is
  mine; both are other agents' work and were left alone.
- `gatk-sv-testkit`: was **DIRTY** at freeze time with changes this session did not make —
  `docs/setup.md`, `scripts/probe_fixes.py`, `scripts/selftest.sh`. This session committed **only**
  `docs/terra-head-to-head.md` (path-scoped, `d63c771`) and never staged, reverted or reviewed the
  foreign files. **They were committed by their author mid-session** as `d958e8b` ("selftest: a
  skipped probe no longer satisfies the pinned count, and says what to install"), so at handoff the
  worktree is **clean**, `ahead=0`, local head == `origin/main` == `d63c771`. §14 #9 closed on that
  evidence rather than being carried forward as a guess.
- `gatk-sv` main worktree: dirty only with untracked tool dirs (`.claude/`, `.serena/`, `.tokensave/`,
  `wt/`), head `e1909d2f` == `origin/main`.
- Compute: zero `gsv-*` instances in `broad-dsde-methods`.

## 14. Open items / next steps

- [x] **1. Reconcile the two SIGPIPE fixes** — **closed by the 2026-10-08 rebase, but not as the
      union asked for.** `main` `a55c498f` took the `|| true` form (3 sites in
      `wdl/CollectCoverage.wdl`) and the trio variant `4ba78d4b` (drop `exit`) was dropped, so there
      is one fix in one place. The strictness trade went the other way: `|| true` swallows every
      failure of that pipeline, not just rc 141. See #13 for the test that did not come with it.
- [x] **2. Port the SIGPIPE fix to main** — **done, and not by this branch.** `wdl/CollectCoverage.wdl`
      on `origin/main` `a55c498f` now guards the header probes; verify with
      `git log -1 -- wdl/CollectCoverage.wdl` on `main`.
- [x] **3. Fix defect D** — **fixed in the 2026-10-08 review-fix round** (`886d93c7`).
      `AddTrioSamplesToPed` now writes both provided parents as founders (`0` in both parent
      columns), so the full-trio PED is a pedigree instead of an `F`↔`M` 2-cycle. The task also
      gained a self-check over the finished file that fails the task if a declared parent row is
      not a founder, and a role check: a declared mother must show chrX CN 2, a declared father
      CN 1 (field 2 of `sample_sex_assignments.txt.gz` is the chrX copy number per
      `src/WGD/bin/estimatePloidy.R:643`, not a sex code, so the old 1-or-2 test could not tell
      the parents apart — transposing the two ids used to pass and invert every `INHERITED_FROM_*`
      label downstream). Verified by extracting the real command block out of the WDL and running
      it against a synthetic panel PED + ploidy tarball for all three shapes (`CASE F M` /
      `M 0 0` / `F 0 0`, mother-only `CASE 0 M` + `M 0 0`, father-only `CASE F 0` + `F 0 0`),
      plus two negative runs (swapped ids → exit 1; founder logic reverted → the self-check exits
      1). Half-trio rows are unchanged from before. `miniwdl check` exits 0.
      Remaining gap: the harness lived in `/tmp`, so nothing in the repo pins these three row
      sets — the in-task self-check is the durable guard, not a test.
- [ ] **4. Same localization trap on the remaining `File?` inputs** of `ValidateTrioInputs`
      (`dragen_vcf`, `case_*_vcf`) — Booleans instead, as done for the CRAMs. **Still open after the
      rebase, and now counted:** `wdl_semantics.py` reports 5 `DEFINED-ONLY` findings on this branch
      against 2 on `main`, exactly these inputs (`dragen_vcf`, `case_manta_vcf`, `case_melt_vcf`,
      `case_scramble_vcf`, `case_wham_vcf`), each localized and never opened.
- [ ] **5. Judge run 4's QC** — 10 PASS / 17 FAIL against the shipped `single_sample.qc_def`.
      Diff this run's QCDF against a v1.1.1 single-sample baseline QCDF to separate expected trio-mode
      drift from regression. Not done: no baseline QCDF was ever pulled.
- [ ] **6. `final_bed`** — decide whether the single-sample WDL should re-declare it or the published
      `help-gatk/gatk-sv-single-sample` @v1.1 config should drop it. Until then every run driven by
      that config reports `Failed` after success.
- [ ] **7. No MOI metrics coverage** — `SingleSampleMetrics` has 0 trio/MOI keys; nothing would catch
      an MOI regression except reading `moi_summary` by hand.
- [ ] **8. Confirm run 5 landed `Succeeded` with `moi_summary` registered** (§8/§9) and record the
      verdict here; if it is still `Running`, do not report it as done. Its MOI table must match
      §1 exactly — if `RefineComplexVariants`/`AnnotateModeOfInheritance` re-ran (which §7 item 5 says
      it did), any difference in the counts is a real finding about run-to-run reproducibility, not noise.
- [x] **9. Account for the dirty `gatk-sv-testkit` files** in §13 — **resolved without this session's
      involvement**: their author committed them as `d958e8b` while this session was writing the doc;
      the worktree is clean at handoff. Left checked-out here so the next session can confirm nobody
      is still editing them.
- [ ] **10. RD cutoffs are still not quotable** — the only full-interval run is this trio run, whose
      QC is not validated (item 5), and cutoffs are meaningless without a baseline comparison.
- [ ] **11. Report the tool defects** — `terra-monitor/scripts/twatch.py --diagnose` → HTTP 405;
      Browser tool's BetterChromium wrapper needs `bun` (not installed); **`session-handoff` ships
      `SKILL.md` text that references `scripts/handoff_scan.sh`, which does not exist** (only
      `repo_state.sh` and `references/HANDOFF_TEMPLATE.md` are installed), so the skill's own
      "prove it before declaring the handoff done" step cannot be run and its parser properties had to
      be verified by hand instead.
- [ ] **12. Nothing on the `RunCNVNonGenotyper` `/gatk/gatk` claim** — that task has still never
      executed (GD skipped in trio mode). It is *untested*, not *fixed*; see §7 item 3.
- [ ] **13. The SIGPIPE regression test has no branch of its own** — `scripts/test/test_sigpipe.sh`
      came from `mw_fix_single_sample_blocking` (`04fa5142`, `083e9956`). `main` took the guard but
      not the test (`git ls-tree -r main | grep -i sigpipe` is empty), and that branch was deleted
      from `origin` on 2026-10-08. The commits are still reachable, through `mw_genotype_scale`
      history and the `tmp/pre-rebase-62d7f80d` backup ref, but nothing owns them. Salvage if wanted:
      `git branch tmp/salvage-test-sigpipe 083e9956`.
- [ ] **14. CI runs no tests at all** — `.github/workflows/pytest.yaml` installs pytest and then runs
      only `tox -e lint`, which is bare flake8 (`tox.ini [testenv:lint]`). So `test_annotate_moi.py`
      and the new `test_convert_cnvs_without_depth_support_to_bnds.py` are never collected: PR #971's
      `Linting` check passed in 22 s without touching either. Both files were run locally, in a
      throwaway venv with `pysam` + `pytest`, to get the 8 passing results. Fixing the runner is a
      3-line change, but it gates every file in the repo and so is deliberately NOT part of this PR;
      it needs its own branch and a decision about which test trees to enable.
- [ ] **15. `sv_shell` silently ignores trio inputs** — `src/sv_shell/single_sample_pipeline.sh` has
      no `mother_*`/`father_*` handling (its README scopes it to case mode), so an input JSON carrying
      `mother_cram` + `mother_sample_id` is accepted, never read, and yields a case-only VCF with no
      MOI and no error. Three lines next to the other input reads would turn that into a refusal. Out
      of scope for PR #971 by provenance (`src/sv_shell/` is untouched by it); recorded so it is not
      lost. Its copy of the case genotype filter is independent of the WDL dedupe in `5bffc79a`, so it
      still names its output `.<base>.filter_by_<sample_id>_gt.vcf.gz` where the WDL now writes
      `.<base>.filter_by_trio_gt.vcf.gz`.
