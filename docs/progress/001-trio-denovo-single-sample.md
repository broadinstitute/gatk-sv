# 001 — Trio de novo mode for the single-sample pipeline: implemented, two review rounds, awaiting first real-data run

**Session date:** system clock `2026-09-21 15:18 EDT` at handoff (repo declares no date source).
Work spanned three days per git: commits dated 2026-09-11, 2026-09-13, 2026-09-21 (author dates
equal committer dates; verified by `git log --format='%h ad=%ad cd=%cd'`).

No prior progress doc exists for this repo (this is 001; the repo has no `CLAUDE.md`/`AGENTS.md`
and no numbered docs series — see Open items). Session = implementation + owner adversarial
review + delegated (subagent) adversarial review + a no-op review of `.dockstore.yml`.

---

## 1. Trio de novo feature (workstream A) — implemented, validated statically, never run

Branch `trio_denovo_single_sample`, 3 commits on base `857419a0` (local `mw_gd_external`):

- `8995b2d0` feature: optional mother/father CRAM+`sample_id`(+index) inputs to
  `wdl/GATKSVPipelineSingleSample.wdl`; per-parent `GatherSampleEvidence` calls; trio arrays feed
  GatherBatchEvidence/ClusterBatch/genotyping; all trio members' alt calls retained in final VCF;
  MOI annotation as final step (`wdl/AnnotateModeOfInheritance.wdl` +
  `src/sv-pipeline/05_annotation/scripts/annotate_moi.py` + pytest; INFO `MOI`, `MOI_CONFIDENCE`);
  multi-sample filter tasks in `wdl/SingleSampleFiltering.wdl`; `AddTrioSamplesToPed` +
  trio-aware CNMOPS ped in `wdl/GatherBatchEvidence.wdl`; `website/docs/execution/single.md` section.
- `1763a774` owner-review fixes: half-trio `select_first` crash in MOI inputs; metrics
  `select_first` on gated task; new `ValidateTrioInputs` task (fail-fast; its output file feeds
  every consumer of `trio_samples_list` so it cannot be pruned); GBE `combined_ped_file` batch
  null-semantics restored; sex-code ∈ {1,2} check.
- `a51ae8d6` subagent-review fixes: metrics `extra_samples` input (svtest asserts VCF sample set
  == metrics list — would kill *every* trio run); scramble-requires-manta validation; GD gated off
  in trio mode (`gd_output_tarball` now `File?`); exact `-F'\t'` ploidy lookups (grep `-w`
  false-matches `PROBAND` vs `PROBAND-M`); PED parents get `0 0` (no self-parent rows);
  `while read … || [ -n "$sid" ]` hardening; `pysam.tabix_index(preset="vcf", force=True)`
  (production pins pysam **0.15.4**, `dockerfiles/sv-pipeline-virtual-env/Dockerfile:40` —
  verified by read); output types restored to original (`File` for
  `non_genotyped_unique_depth_calls[_idx]`); docs corrections.

Validation evidence (all re-run at/after final commit, exact strings reproducible):
`pytest … test_annotate_moi.py -q` → `4 passed`; flake8 clean;
`miniwdl_validation.py --imports-dir wdl wdl/*.wdl` → "All the WDLs successfully passed the
miniwdl validation."; `scripts/test/validate.sh -d . -j womtool.jar` (Cromwell 84 womtool) →
"48 TESTS PASSED SUCCESSFULLY!"; womtool vs a synthetic trio input JSON → "Success!"; the
`ValidateTrioInputs` bash rendered from real JSON values + 6 synthetic configs (scratch under
/tmp, ephemeral) behaved correctly on all.

**Not exercised: anything requiring containers or real data** (no docker on this machine).
`annotate_moi.py` ran only on pysam 0.24 locally — pysam **0.15.4** compatibility is inferred from
repo precedent (`sv_utils/fix_vcf.py` etc. call `tabix_index` in that image), not executed.
End-to-end trio run remains the whole untested surface.

## 2. Coordinates of anything created or changed outside this repo

| thing | identifier | how to verify | how to undo |
|---|---|---|---|
| Remote branch `trio_denovo_single_sample` on `origin` | `https://github.com/broadinstitute/gatk-sv.git` | `git ls-remote origin refs/heads/trio_denovo_single_sample` | `git push origin --delete trio_denovo_single_sample` (reversible; pushed by this handoff) |
| Base branch state | `857419a0` = local `mw_gd_external`; **`origin/mw_gd_external` no longer exists** (`git ls-remote origin refs/heads/mw_gd_external` → empty, verified 2026-09-21; tracking ref shows `[gone]`) | `git merge-base trio_denovo_single_sample 857419a0` | n/a — record only; do not "restore" without owner |
| Subagent adversarial review | run `105529fe-8574-48b6-9562-5cd4535b9bf2`; session file `~/.pi/agent/sessions/--Users-markw-IdeaProjects-gatk-sv--/2026-09-11T19-50-29-003Z_01a09205-7f8a-77c5-b0b0-92e5ebb6f257/47e8db6f-3fa9-492f-9291-ad47d3c157fa/run-0/session.jsonl` | read file | n/a |
| Scratch venvs (OS purges /tmp periodically) | `/tmp/moi-venv` (pysam 0.24, pytest, flake8), `/tmp/wdl-venv` (miniwdl, jinja2, +pyyaml installed this session) | `ls /tmp/moi-venv/bin/python` | disposable; rebuild recipe in §3 |
| womtool jar | Cromwell **84**, ~167 MB, downloaded from GitHub releases **and deleted after each validation run** — not in tree | `git status --short` clean | n/a |

**External mutations ledger (this session):** local commits `8995b2d0`, `1763a774`, `a51ae8d6` +
this handoff commit (undo: `git reset --hard 857419a0`, reversible, nothing was amended); push of
the branch (listed above, reversible); `/tmp` fixtures/scripts (ephemeral, self-deleting). No
cloud/Terra/gcs writes. No irreversible mutations.

## 3. Resume here (paste-able)

```
cd /Users/markw/IdeaProjects/gatk-sv/wt/trio-denovo
 git status -sb && git --no-pager log --oneline -4
 git --no-pager log --oneline origin/trio_denovo_single_sample -1   # must equal local head
# rebuild scratch venvs if /tmp was purged:
 python3 -m venv /tmp/moi-venv && /tmp/moi-venv/bin/pip -q install pysam pytest flake8
 python3 -m venv /tmp/wdl-venv && /tmp/wdl-venv/bin/pip -q install miniwdl jinja2
/tmp/moi-venv/bin/python -m pytest src/sv-pipeline/05_annotation/scripts/test_annotate_moi.py -q
/tmp/moi-venv/bin/flake8 --config tox.ini src/sv-pipeline/05_annotation/scripts/annotate_moi.py \
  src/sv-pipeline/scripts/single_sample/convert_cnvs_without_depth_support_to_bnds.py
/tmp/wdl-venv/bin/python scripts/test/miniwdl_validation.py --imports-dir wdl wdl/*.wdl
# womtool matrix (re-download; delete jar before committing):
 curl -sL -o womtool.jar https://github.com/broadinstitute/cromwell/releases/download/84/womtool-84.jar
bash scripts/test/validate.sh -d . -j womtool.jar ; rm -f womtool.jar
# trio-path static check (rebuilds the synthetic trio JSON):
 python3 - <<'PY'
import json
d=json.load(open("inputs/build/NA19240/test/GATKSVPipelineSingleSample.json"))
p="GATKSVPipelineSingleSample."
d.pop(p+"melt_docker",None); d[p+"use_melt"]=False
d[p+"mother_sample_id"]="SAMPLE_MOM"; d[p+"mother_cram"]="gs://t/m.cram"; d[p+"mother_cram_index"]="gs://t/m.cram.crai"
d[p+"father_sample_id"]="SAMPLE_DAD"; d[p+"father_cram"]="gs://t/f.cram"; d[p+"father_cram_index"]="gs://t/f.cram.crai"
open("/tmp/trio_inputs.json","w").write(json.dumps(d,indent=2))
PY
 java -jar womtool.jar validate wdl/GATKSVPipelineSingleSample.wdl -i /tmp/trio_inputs.json
```

## 4. What "good" looks like on the next check

| check | expected | what a mismatch means |
|---|---|---|
| `git status -sb` in worktree | clean, head == `origin/trio_denovo_single_sample` head | uncommitted or unpushed work |
| pytest command above | `4 passed in 0.0Xs` | MOI script/tests regressed |
| miniwdl validation | `All the WDLs successfully passed the miniwdl validation.` | WDL type error introduced (one real catch this session: gated call output needs `select_all([...])` to reach `Array[File]`, miniwdl: "Expected Array[File] instead of Array[File?]+") |
| `validate.sh` matrix | `48 TESTS PASSED SUCCESSFULLY!` | womtool-visible API break vs existing input JSONs |
| womtool + trio JSON | `Success!` | trio inputs broke static contract |
| first real trio run, final metrics step | completes; `metrics_file` counts include parent calls | svtest sample-set assertion (`ValueError: One or more sample(s) found in VCF header but not samples list`) means `extra_samples` wiring regressed |

**Cost / time / size:** not measured — no cloud run occurred. Expect ~3× `GatherSampleEvidence`
cost vs single-case (case + up to 2 parents) plus parents in CNMOPS/gCNV/genotyping; read real
figures from the Terra/Cromwell cost report of the first trio run.

## 5. Gotchas found (actually hit, with the error text)

1. `/tmp/moi-venv` silently purged mid-project → `No module named pytest` / `flake8: No such file
   or directory` → macOS reaps /tmp; rebuild recipe in §3; never assume /tmp survives across days.
2. `svtest vcf` hard-asserts VCF column set == samples list (`src/svtest/svtest/utils/TestUtils.py:8`
   raises `ValueError: One or more sample(s) found in VCF header but not samples list: …`) — any
   mode that adds VCF columns must also extend metrics `extra_samples`.
3. `grep -w "^$sample"` matches `PROBAND` inside `PROBAND-M` (`-`/`.` are word boundaries) →
   reproduced as `Error: ploidy-derived sex code '1\n2\n1' for sample PROBAND is not 1 (male) or
   2 (female)` → use exact `awk -F'\t' '$1==s'` lookups.
4. Conditional-call outputs are optionized outside the `if` block → `[GD.gd_output_tarball]` typed
   `Array[File?]`; fix is `select_all([...])` inside the `defined()` gate.
5. The repo's own `inputs/build/NA19240/test/GATKSVPipelineSingleSample.json` is
   runtime-invalid (`use_melt: true` with no `melt_docker` — pre-change code would crash on it in
   `select_first`; CI `validate.sh` is **static womtool only**). The new `ValidateTrioInputs` now
   rejects it with `ERROR: use_melt requires melt_docker or case_melt_vcf` — that rejection on the
   test JSON is *correct*, not a false positive.
6. BSD/macOS tool quirks in this shell: `sed -n 'N,+40p'` invalid; plain `grep` here may be a
   wrapper that truncates/annotates output — quote `/usr/bin/grep` for exact matching.

## 6. Corrections to earlier documents

- **Was:** STRipy-merged records can be `DE_NOVO` in a half-trio. **Now:** merge clears GT for all
  samples (`merge_stripy_vcfs.py:204` `record.samples[sample]['GT'] = (None, None)`) so they are
  always `UNASSESSABLE`/`UNCONFIRMED` — docs and script updated in `a51ae8d6`. Mattered: it shaped
  the MOI confidence design claim.
- **Was (commit 8995b2d0):** gating `GetUniqueNonGenotypedDepthCalls` to non-trio is safe;
  widened output type to `File?`. **Now:** its consumer `SingleSampleMetrics` evaluates
  `select_first([…])` unconditionally → trio runs would abort; task unconditional again, output
  type restored to original `File`. Evidence: original lines in `git show mw_gd_external:…`
  (`File non_genotyped_unique_depth_calls = GetUniqueNonGenotypedDepthCalls.out`).
- **Was (owner round-1 review conclusion):** "non-trio behavior unchanged in semantics".
  **Now:** `select_all` generalization had silently swallowed missing-caller-docker
  misconfigurations the original crashed on; fixed by making `ValidateTrioInputs` enforce
  `use_X ⇒ X_docker || case_X_vcf` in **all** modes (contract encoded from the original
  `select_first` crash sites).
- **Was (inherited session note):** "Cromwell 84 lacks `sep()`". Proven wrong this session —
  production WDL uses `~{sep=" " raw_vcfs}` (`wdl/PESRPreprocessing.wdl:78`). The `tr '\n' ' ' <
  write_lines(...)` renderings added earlier work regardless; claim retired, origin unverifiable
  (predates compaction).

## 7. Deliverables

| file | change |
|---|---|
| `wdl/GATKSVPipelineSingleSample.wdl` | trio inputs/evidence/arrays, `ValidateTrioInputs`, MOI wiring, GD gated off in trio, metrics `extra_samples`, output-type fixes |
| `wdl/SingleSampleFiltering.wdl` | multi-sample genotype filters, proband-generalized CNV→BND task, loop hardening |
| `wdl/GatherBatchEvidence.wdl` | `AddTrioSamplesToPed` (exact matching, valid PED rows), trio-aware CNMOPS ped, batch-mode output parity |
| `wdl/AnnotateModeOfInheritance.wdl` | new MOI task (no unused `vcf_idx`) |
| `wdl/GATKSVPipelineSingleSampleMetrics.wdl` | `extra_samples` input (CRASH-1 fix) |
| `src/sv-pipeline/05_annotation/scripts/annotate_moi.py` (+ test) | MOI classification, summary TSV, pysam-0.15-safe indexing |
| `src/sv-pipeline/scripts/single_sample/convert_cnvs_without_depth_support_to_bnds.py` | multi-proband support |
| `website/docs/execution/single.md` | trio section + limitations + corrected MOI/STRipy/output-table facts |
| `docs/progress/001-trio-denovo-single-sample.md` | this document |

**Commits / pushes:** gatk-sv (worktree `wt/trio-denovo`) → `8995b2d0`, `1763a774`, `a51ae8d6`, +
handoff commit; pushed to `origin/trio_denovo_single_sample`, remote head re-read with
`git ls-remote` and confirmed equal to local head (see §Resume; result recorded in the handoff
message, and re-verifiable there).
**Other repos in the inventory (untouched, for completeness):** main checkout `~/IdeaProjects/gatk-sv`
on `mw_gd_external` @ `857419a0` (untracked personal dirs only); `wt/manta_tloc_autoresolve`
(clean vs `origin`, one untracked scratch file `test_single_tloc.py`); two worktrees under
`~/Work/genotypebatch_debug/` belong to a different workstream.

## Open items / next steps

- [ ] **Run the first end-to-end trio test on real data** (Terra or Cromwell on local CRAM trio) —
      everything above is static/unit validation only. Success = final VCF carries parent columns,
      every record has `MOI`/`MOI_CONFIDENCE`, `moi_summary` written, metrics step completes.
- [ ] Verify `annotate_moi.py` inside the real image (pysam **0.15.4**): `tabix_index`,
      `header.info` membership, INFO String writes were not executed there.
- [ ] Owner decision: GD-in-trio — keep the skip, or implement a case+panel RD-matrix subset task
      so GD can run in trio mode (would need bincov column subsetting + tabix regeneration).
- [ ] Owner decision: rename branch to the `mw_*` convention before PR (origin currently hosts
      `mw_gd*`, `mw_manta_*`; also note `origin/mw_gd_external` was deleted upstream — confirm
      the intended PR base).
- [ ] Before PR: delete/relocate `docs/progress/` (session artifact) and decide on a trio example
      input JSON for `inputs/` (none added).
- [ ] Upstream hygiene (out of scope, worth PRs to `main`): fix the runtime-invalid NA19240 test
      JSON (`use_melt` without `melt_docker`); add the missing final newline in
      `.github/.dockstore.yml` (branch copy verified byte-identical to `main`, blob
      `70eb0d249a1de9affe507410676392716e610e05`; file is readonly-for-PRs per #922 — no edit was
      made here).
- [ ] Consider trio-mode QC threshold guidance: `metrics_file`/`qc_file` variant counts and `AF`
      now include parents (documented, thresholds untouched).
- [ ] Note for future sessions: `handoff_scan.sh` referenced by the session-handoff skill is not
      present in the install (`scripts/` contains only `repo_state.sh`); structure verified by eye
      against heading anchors instead.
