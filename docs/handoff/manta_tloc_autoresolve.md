# 001 — Manta tloc auto-resolve: implemented, adversarially reviewed, opt-in scoped; PR + docker validation remain

**Session date:** 2026-09-21 (system `date -u` = `2026-09-21T19:18Z`; consistent with the pi
workflow receipt timestamps `2026-09-21T14:31Z` for the same session's work — clocks agree).
No prior handoff doc exists for this repo; numbering starts here.

> **SUPERSEDED IN PART (2026-09-25) — read `002_manta_tloc_autoresolve_pr_and_terra_ab.md`.** The title's remaining work is done: **PR #968 is open** (onto #966's branch, 3 files), and the flag-on production-stack path that this doc listed as untested has now been run twice and PASSes (3 samples on the real `sv-pipeline` image; 156 samples on Terra). Also stale: merge-base with `main` is no longer `4bc70a69` — the branch was rebased onto `mw_fix_single_sample_blocking` (`4419315c`), and `test_single_tloc.py` was decided *out* of the PR, not in.

Session arc: implemented the "Manta tloc insertions should be automatically resolved" TODO on
branch `mw_manta_tloc_autoresolve`, then ran a three-lane adversarial subagent review that found
two majors (strand-class semantics) and one blocker (blast radius), all fixed and pushed.

---

## 1. Feature workstream — `mw_manta_tloc_autoresolve` (3 commits on merge-base `4bc70a69`)

**b6c6e8de** — `CANDIDATE_SINGLE_TLOC`: a single interchromosomal manta/dragen BND carrying
`CHR2`/`END2` (standardization keeps one record per MATEID pair) resolves to `SVTYPE=CTX` inside
`svtk resolve` instead of being dropped as `SINGLE_ENDER` → `UNRESOLVED` by `mantatloccheck.sh`.
File: `src/svtk/svtk/cxsv/complex_sv.py`.

**b0cd91e4** — review fix, the label/cross-check contract:
- CPX_TYPE label from cytoband arms (same arms → `CTX_PP/QQ`, different → `CTX_PQ/QP`).
- Strand **class** cross-check: `+-`/`-+` must be same-arm, `++`/`--` must be cross-arm;
  contradictions demote to `UNRESOLVED_TYPE=<label>_MISMATCH` (suffixed, stamped on constituent
  records), mirroring `resolve_translocation` + `cpx_tloc.classify_simple_translocation` +
  `reformat_CPX_bed_and_generate_script.py:426-431` PE-orientation queries. First cut accepted
  only antiparallel records → would have silently dropped every pericentric p↔q tloc.
- Resolved record rebuilt via `rec.copy()` from the surviving record (pass-2 SR-only shrink
  provenance; `pysam 0.15.4` `copy()` = `bcf_dup`, no shared `bcf1_t` — reviewer source-dived).
- Missing `STRANDS` → `STRAND_MISMATCH_TLOC` demotion (no KeyError);
  cytoband fetch guarded by `except (StopIteration, ValueError, IndexError)` → `CTX_UNR`
  (a bad band row previously escaped the `resolve_complex_sv` generator as PEP-479 RuntimeError).

**4d673b7f** — blocker fix, blast radius: auto-resolve is now **opt-in**.
`ComplexSV(..., resolve_single_tlocs=False)`; new CLI flag `svtk resolve --resolve-single-tlocs`;
only `src/sv-pipeline/00_preprocessing/scripts/mantatloccheck.sh` passes it. Verified *by execution*
(real `svtk resolve` CLI on manta/dragen single-BND fixture): flag-off → 0 resolved, all 3 records
`UNRESOLVED_TYPE=SINGLE_ENDER`, identical to `main`. Without this, cohort `ResolveComplexVariants`/
`ResolveCpxAll` + sv-shell + single-sample would have flipped every interchromosomal manta/dragen
BND site to `FILTER=PASS CTX` (per-contig scatter ⇒ those mates never cluster; see §5 gotcha 5).

**Review panel dispositions** (full findings lived in this session's transcript; substance folded
here — /tmp scratch is gone):
- Logic lane: majors 1–2 fixed by b0cd91e4; minors 3–5 fixed (rec.copy, KeyError/ValueError,
  bool+comma-split). Its "paired path stamps unsuffixed UNRESOLVED_TYPE" line-ref claim was
  **wrong** — the suffixed stamp via resolve.py defaults was verified by trace and adopted.
- Blast lane: blocker fixed by 4d673b7f. Its other majors are consequences of the blocker
  (cohort RD-QC exemption flood, `update_variant_representations` CTX→4×BND path) → dormant
  behind the flag, ticketed as open items. eph `ClusterTloc.wdl`
  (`git show origin/eph_cluster_tloc_review_pe:wdl/ClusterTloc.wdl`) verified to accept emitted
  records as-is (`INFO END=` parse, SVTYPE filter, arm relabel consistency).
- Runtime lane: no blocker on pinned `pysam 0.15.4` (copy semantics, `END<POS` write parity,
  INFO type matrices). Adopted: `IndexError` catch; production `Number=.` ALGORITHMS header +
  fresh fixture loads in the harness.
- *Not exercised locally:* full flag-on CLI end-to-end (local pysam 0.24 cannot write resolved
  records through resolve.py's `bcftools sort` stdin pipe; see §5 gotcha 3) — **requires the
  sv-pipeline docker rebuild**. Flag-on record-level semantics are covered by the harness instead.

## 2. Coordinates of anything created or changed outside this repo

| thing | identifier | how to verify | how to undo |
|---|---|---|---|
| Remote branch (reversible) | `origin/mw_manta_tloc_autoresolve` @ `4d673b7f` at last check (b6c6e8de, b0cd91e4, 4d673b7f + this handoff commit) | `git ls-remote origin refs/heads/mw_manta_tloc_autoresolve` | `git push origin :mw_manta_tloc_autoresolve` |
| Worktree | `/Users/markw/IdeaProjects/gatk-sv/wt/manta_tloc_autoresolve` | `git worktree list` | `git worktree remove` |
| Test venv (untracked) | `wt/manta_tloc_autoresolve/.venv-tloc` (Python 3.12.13, pysam 0.24.1 + setuptools<71 shim, editable svtk) | `.venv-tloc/bin/python --version` | `rm -rf .venv-tloc` |
| Harness (UNTRACKED) | `wt/manta_tloc_autoresolve/test_single_tloc.py` (19 checks) | `git status` in worktree | delete file |
| Broken async workflow | pi run `2df9d1ba-…` (3 failed children, infra-only) | `subagent status` (dir in receipt) | n/a (already terminated) |

Untouched other-agent trees (informational, last observed at handoff): main tree
`~/IdeaProjects/gatk-sv` [mw_gd_external] `857419a0` clean of tracked changes (local-only branch,
not on origin); `wt/trio-denovo` [trio_denovo_single_sample] `a51ae8d6`;
`~/Work/genotypebatch_debug/wt/gatk-sv-scale` `1fe2d87e`; `…/gatk-sv-v111` `b5d5049c`.

## 3. Resume here (paste-able)

```
cd /Users/markw/IdeaProjects/gatk-sv/wt/manta_tloc_autoresolve
git fetch origin && git status -sb && git --no-pager log --oneline main..HEAD
 .venv-tloc/bin/python test_single_tloc.py        # expect: 19 PASS lines + ALL PASS
git ls-remote origin refs/heads/mw_manta_tloc_autoresolve   # must equal local HEAD
```

## 4. What "good" looks like on the next check

| check | expected | what a mismatch means |
|---|---|---|
| `test_single_tloc.py` | exactly 19 `PASS` lines, zero `FAIL`, final `ALL PASS` (two `NOTE …` lines about pysam ≥0.22 `stop<POS` clamp are normal on the local venv) | regression in resolve_single_tloc semantics or harness drift |
| remote vs local head | `ls-remote` sha == `git rev-parse HEAD` | a push was lost / branch moved elsewhere |
| flag-off e2e (real CLI, manta single BNDs) | 0 resolved records; every input `UNRESOLVED_TYPE=SINGLE_ENDER` | opt-in gate leaked; blast-radius bug reintroduced |
| flag-on docker e2e (**RUN 2026-09-23, PASS** — see 002 §1b/§1c: 3 samples on the built image, 156 on Terra; `CTX` 172/160/155 docker, +27,267 cohort; zero other label moved) | same-arm `+-`→`CTX_PP/QQ`; cross-arm `++`→`CTX_PQ/QP`; `+-`-cross-arm → `CTX_PQ/QP_MISMATCH`; wham stays `SINGLE_ENDER`; no POSTHOC duplicates | production-stack divergence (pysam 0.15.4 path untested locally) |
| `git merge-base HEAD origin/main` | currently `4bc70a69`; `origin/main` is 8 commits ahead at handoff | if merged past, PR diff changes — re-review touched files |

**Cost / time / size:** pipeline runtime impact — **now measured (2026-09-23), was not measured when this line was written**: +2.4% cpu-VM-minutes on the tloc task (392.5 -> 402.1 for 156 samples) and +2.8% cohort complex output; full step-04 rerun is 31,206 VM-minutes, its `TinyResolve` alone 339. See 002 §1c. (Relevant because dockerfiles/sv-pipeline-virtual-env pins pysam==0.15.4.)

## 5. Gotchas found (hit this session, with the error text)

1. Async subagent lanes dead: `Error: Cannot find module
   '/Users/markw/.pi/agent/npm/node_modules/pi-subagents/src/runs/background/subagent-runner.ts'`
   (Node v26.4.0) → pi-subagents 0.70.1 npm package ships only `subagent-runner.js`, but spawn
   forces `PI_ASYNC_NATIVE_RUNNER=1` (`async-execution.js:503`) which needs the `.ts` → foreground
   sequential subagent lanes worked fine; background/async needs an upstream/package fix or
   non-Node-26 runner.
2. `ComplexSV` record mutation bleed: sharing one parsed-records dict across cases leaked
   `UNRESOLVED_TYPE` stamps between tests → reload the fixture per case.
3. Flag-on full-CLI e2e impossible on local pysam 0.24: `[E::vcf_format] Invalid BCF, the INFO
   tag id=16 is too large` when writing any record carrying a header-`add_line`-added tag through
   `resolve.py:444`'s `pysam.VariantFile(pipe.stdin,'w')`; wrapped as `TypeError: expected str,
   bytes … not BufferedWriter`. Pre-existing env limitation (production 0.15.4 ships this exact
   flow). Path-based writes of the same records work.
   ALSO: never `proc.wait()` on the pipe consumer while the pysam file object is alive —
   that deadlocked a throwaway probe for ~2 h (harmless, killed; not pipeline behavior).
4. `svtk resolve -q` expects a value (argparse has no `store_true`) — `-q true`; `bgzip -c`
   writes to stdout (use `bgzip -f` for files); tabix rejects empty beds and `/dev/null`.
5. ResolveCpxSv scatters per contig (`ResolveComplexVariants.wdl:66-90` + PullVcfShard) ⇒
   interchromosomal mate pairs can never cluster in the cohort path; treat that path as
   single-record-only when reasoning about tlocs.
6. `git ls-remote` needs auth here (interactive auth already available; fine).

## 6. Corrections to earlier documents/claims

- **Was:** (commit b6c6e8de message + code docstring) "STRANDS on a standardized record is not
  reliable for the strand-based PP/QQ-vs-PQ/QP split". **Now:** the strand *class*
  (`+-`/`-+` vs `++`/`--`) is mate-order- and swap-invariant and IS used as a consistency gate;
  arms supply the label. Evidence: `cpx_tloc.py:75-96`, `resolve_translocation`,
  `reformat_CPX_bed_and_generate_script.py:426-431`; harness cases `++ cross-arm → CTX_PQ/QP`
  (b6c6e8de's gate would have DROPPED pericentric tlocs — the wrong claim had already shaped the
  first implementation). Correction is named in b0cd91e4's commit message.
- **Was:** (prior session summary) "fix applies to both manta and dragen via ALGORITHMS gate".
  **Still true**, but the gate is intersection-based: merged sites like `ALGORITHMS=manta,wham`
  DO resolve (documented + harness case `merged manta+wham`). Deliberate; revisit if undesired.
- Reviewer-claim correction (no repo text): blast/runtime reviewers each asserted something the
  trace refuted (paired-path unsuffixed stamp; `update_best_genotypes` scope) — resolved in-session
  before acting on them.

## 7. Deliverables

| file | change |
|---|---|
| `src/svtk/svtk/cxsv/complex_sv.py` | CANDIDATE_SINGLE_TLOC machinery (+opt-in param, +arm/strand contract, +robustness) |
| `src/svtk/svtk/cli/resolve.py` | `--resolve-single-tlocs` flag threaded through `resolve_complex_sv`/`_v2` and all 4 `ComplexSV(...)` sites |
| `src/sv-pipeline/00_preprocessing/scripts/mantatloccheck.sh` | sole caller passing the new flag |
| `test_single_tloc.py` (worktree, UNTRACKED) | 19-case regression harness incl. flag-off regression + CLI merge/sanity simulation |
| `docs/handoff/manta_tloc_autoresolve.md` (this file) | handoff doc; drop/exclude this commit before PR merge if undesired |

**Commits / pushes:** repo `~/IdeaProjects/gatk-sv` (origin = github.com/broadinstitute/gatk-sv) →
branch `mw_manta_tloc_autoresolve`, feature heads b6c6e8de → b0cd91e4 → 4d673b7f; remote head was
re-read at each push (`4d673b7f…` confirmed; re-confirm after this handoff commit).
**External mutations made:** created+pushed `origin/mw_manta_tloc_autoresolve` (reversible: delete
branch); created worktree + `.venv-tloc` (reversible: rm); `/tmp` scratch dirs
(`manta_tloc_test_*`, `mwa_review`, `tloc_e2e`) — all deleted, reported as deleted. No pipeline
jobs, no external services touched.

## Open items / next steps

- [x] DECIDED 2026-09-23: **no test file in the PR** (owner decision; this repo has no python test directory to follow — `find` for `test_*.py` hits only my own untracked file), evidence goes in the description; harness stays untracked. Offered as `src/svtk/svtk/test/` if reviewers ask. Original recommendation, now overridden: yes — `git add` in the
      worktree, amend or new commit, push; verifies: file listed in `git show --stat HEAD`).
- [x] DONE 2026-09-23 as **PR #968**, but not this shape: opened from a separate branch `mw_manta_tloc_autoresolve_pr` stacked onto `mw_fix_single_sample_blocking` (#966), because this branch is rebased onto that unmerged branch and a PR to main would have carried its 8 files too; harness WDL + docs/handoff excluded (002 §1a). Original text: Open PR from `mw_manta_tloc_autoresolve` (base `main` @ `e1909d2f`, 8 ahead of merge-base
      `4bc70a69` — decide rebase vs leave; if rebasing, re-run the harness afterwards).
- [x] DONE 2026-09-23 — Rebuilt sv-pipeline docker from this branch and ran the whole TinyResolve/ResolveManta set end-to-end; flag-on behavior confirmed as §4 row 4 predicted (002 §1b/§1c). Original text: Rebuild sv-pipeline docker from this branch and run one TinyResolve/ResolveManta shard
      end-to-end; expect flag-on behavior per §4 row 4 (this is the only untested production path).
- [ ] Follow-up ticket (pre-existing, now dormant behind the flag):
      `src/sv-pipeline/scripts/single_sample/update_variant_representations.py` CTX→4×BND expansion
      — wrong REF anchor for M4 (~line 115) and self-referential `MATEID` (~line 157); blast lane
      reproduced both. Do not enable `--resolve-single-tlocs` in the single-sample path until fixed.
- [ ] Confirm with owner: within the manta-tloc workflow the `EVIDENCE=PE` gate is inert
      (mantatloccheck stamps PE on every record, incl. SR-only) → `manta_tloc` evidence volume now
      equals all interchromosomal manta BNDs. Assumed intended for the eph tloc workflow
      (it does its own AF/filtering) — needs the eph owner's ack.
- [ ] Fix pi-subagents background lanes (upstream: native-runner resolution wants `subagent-runner.ts`;
      local workaround: none found this session — foreground lanes fine).
