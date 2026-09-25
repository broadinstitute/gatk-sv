# 002 — Manta tloc auto-resolve: PR #968 open, flag validated in docker (3 samples) and Terra (156 samples), both PASS

**Handoff written:** `2026-09-25T16:29Z`. Clock agreement verified, not assumed:
`date -u` = `terra (api.firecloud.org HTTP Date)` = `github` = `Fri, 25 Sep 2026 16:28:49 GMT`.
The A/B evidence below was produced **2026-09-23** (Terra submission timestamps) — the session spans
two days of real wall clock, there is no clock disagreement to reconcile.

Predecessor: `001` = `docs/handoff/manta_tloc_autoresolve.md` on this branch. Read it for the design
contract; this doc records what happened after it and **corrects three of its claims** (§6).

---

## 1. Workstreams

### 1a. PR is open — #968 (the thing to review)

`https://github.com/broadinstitute/gatk-sv/pull/968` — base `mw_fix_single_sample_blocking`
(#966's branch), head `mw_manta_tloc_autoresolve_pr` @ `ad80380f`, `OPEN`, not draft,
**3 files, +136/−11** (`gh pr view 968 --json changedFiles,additions,deletions`).

Stacked deliberately: this branch was rebased onto #966, and #966 is still `OPEN`/`MERGEABLE`, so a
PR to `main` would have carried 10 commits including #966's 8 unrelated files. Owner decision this
session (asked explicitly): stacked now, retarget to `main` when #966 merges.

Verified before pushing (`wt/tloc-pr`): `git diff --numstat origin/mw_fix_single_sample_blocking...HEAD`
= exactly the 3 code files; `git diff mw_manta_tloc_autoresolve HEAD -- <those 3 files>` **empty**, so
the PR tree is byte-identical to the tree both A/Bs validated; `python3 -m py_compile` rc=0 on both
python files; `bash -n` rc=0 on a 33-line (non-empty) `mantatloccheck.sh`.

Excluded from the PR by owner decision: the harness WDL + its Dockstore entry, `docs/handoff/`
(content moved into the PR description), and any test file.

`gh pr checks 968` → `no checks reported` — **expected, not a failure**: `docs.yml` triggers on
`pull_request` to `main` only, and the image workflow builds on `push: main` only. No automated
signal exists on this base; the manual build + A/B in the description is the evidence.

### 1b. Docker A/B on the rebased base — 3/3 PASS

Image built from `bc89a7c9` (branch tip at the time, base `4419315c`):
`us.gcr.io/broad-dsde-methods/markw/gatk-sv/sv-pipeline:mw-manta-tloc-autoresolve-bc89a7` =
`sha256:d3ab3af0c04ced200d540b79f24679cd7f2bbffca50055e698a1632e4f273244` (verified by `docker pull`
+ `docker image inspect`). Build **25.4 VM-minutes**. Baseline arm image
`us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-10-02-v1.1-483973d6` is byte-identical to what
the production run used. `pysam 0.15.4` in both arms (the pinned production stack, which is what
`001` §5 gotcha 3 could not reach locally).

| sample | COMPLEX | `CTX` (PP/QQ + PQ/QP) | `SINGLE_ENDER` | mismatches | closure |
|---|---|---|---|---|---|
| HG00096 | 6007 → 6179 | 0 → 172 (98 + 74) | 360 → 42 | 69 + 77 | `172+146 = 360-42` |
| HG00129 | 5458 → 5618 | 0 → 160 (90 + 70) | 321 → 46 | 44 + 71 | `160+115 = 321-46` |
| HG00140 | 6066 → 6221 | 0 → 155 (81 + 74) | 337 → 46 | 70 + 66 | `155+136 = 337-46` |

`compare_newbase.py`: 3/3 PASS on CONTROL (baseline == frozen production artefact, 0/18 labels
differing) / CLOSURE / PURITY. Same numbers as the pre-rebase base (`compare_ab.py`).

### 1c. Terra head-to-head, whole cohort — PASS (this is new evidence, `001` never had it)

Trimmed workflow `TlocResolveOnly` = step 04's tloc path only (`UntarFiles` → `GetShardInputs`
inlined as pure WDL → `ResolveManta`). `build_tloc_wdl.py` extracts `struct RuntimeAttr` +
`task UntarFiles` + `task ResolveManta` **verbatim** from the repo and re-reads its own output to
assert byte-identity; `miniwdl check` rc=0. Both arms differ only in `sv_pipeline_docker`,
`useCallCache=False`, `deleteIntermediateOutputFiles=False`.

| set | complex records | `CTX_PP/QQ` | `CTX_PQ/QP` |
|---|---|---|---|
| production run's own outputs | 961,550 | 6 | 1 |
| baseline arm (prod image) | 961,550 | 6 | 1 |
| branch arm | 988,817 | 15,117 | 12,157 |

`compare_terra.py`: **CONTROL 156/156** samples identical to production on records + every label;
**EFFECT** all 156 samples change with `+27,267` records == added `CTX` sample-by-sample;
**PURITY** zero non-`CTX` label moved. Verdict PASS.

Runtime/size impact now measured (`001` had it as "not measured"): the tloc task
**392.5 → 402.1 cpu-VM-minutes (+2.4%)** for the same 156 samples; cohort complex output **+2.8%**.
Both arms also had exactly one `ResolveManta` shard fail attempt 1 and succeed on retry — shard 0 in
baseline, shard 6 in branch, i.e. unrelated shards ⇒ preemption, not image-specific.

**Scope limit, do not overread:** the trimmed workflow outputs only the **complex** VCF, so
`SINGLE_ENDER` counts and the closure identity are checkable only in the Docker leg (1b). The trim is
sound because `TinyResolve` consumes none of the dropped tasks' outputs (they are siblings, not
ancestors) — but it cannot detect an interaction with them, and it is not production step 04.

## 2. Coordinates of anything created or changed outside this repo

| thing | identifier | how to verify | how to undo |
|---|---|---|---|
| PR | #968 base `mw_fix_single_sample_blocking` | `gh pr view 968 --repo broadinstitute/gatk-sv --json state,changedFiles` | `gh pr close 968` (leaves branch) |
| Branch (pushed) | `mw_manta_tloc_autoresolve_pr` @ `ad80380f` | `git ls-remote origin refs/heads/mw_manta_tloc_autoresolve_pr` | `git push origin :mw_manta_tloc_autoresolve_pr` |
| Branch (pushed) | `mw_manta_tloc_autoresolve` @ `78dd16ab` (harness WDL + Dockstore entry + 001/002 docs) | same command, other branch | force-push back to `42b9877c`, or delete branch |
| Local safety ref | `tmp/pre-rebase-cb97242a` → `cb97242a` (pre-rebase state, local only) | `git rev-parse tmp/pre-rebase-cb97242a` | `git branch -D` |
| Worktrees | `wt/manta_tloc_autoresolve` (78dd16ab), `wt/tloc-pr` (ad80380f) | `git worktree list` | `git worktree remove` |
| Terra workspace | `broad-firecloud-dsde-methods/GATK-SV-manta-tloc-mw-2026-09-21`, bucket `gs://fc-5fb1cc39-5691-424e-8a52-622c1cc7beee` | `terra_probe.py workspaces` | delete workspace (orphan bucket needs separate delete) |
| Method configs | `04t-baseline`, `04t-branch` (root type `tloc_run`) | `GET .../methodconfigs/{ns}/04t-branch` | `DELETE` same path (204) |
| Root entity | `tloc_run/run1` | `GET .../entities/tloc_run/run1` | `DELETE .../entities/tloc_run/run1` |
| Submissions | baseline `48ca3774-dcf9-4c74-9a22-ba8df52f9b0b` (wf `443d43a0-…`), branch `e7d25c1e-11e0-4092-9e7e-972be686e631` (wf `3e74a180-…`) — both `Succeeded` | terra-monitor `status <id>` | abort is moot (terminal); outputs stay in the bucket |
| Frozen inputs | 317 objects / 8.65 GiB under `gs://fc-5fb1cc39-…/tloc_frozen/`, crc32c+size **317/317** (`freeze_verified.json`) | `python3 freeze_copy.py verify` | `gcloud storage rm -r gs://fc-5fb1cc39-…/tloc_frozen` |
| Dockstore | `github.com/broadinstitute/gatk-sv/TlocResolveOnly`, branch `mw_manta_tloc_autoresolve`; **owner published it manually** (set default version + Publish) | authoritative check = POST a throwaway config with the `dockstore://` ref (`201` resolves / `404 Cannot get` not published), then `DELETE` it | unpublish / delete entry |
| Registry image | `markw/gatk-sv/sv-pipeline:mw-manta-tloc-autoresolve-bc89a7` (`sha256:d3ab3af0…`) | `gcloud container images describe` | `gcloud container images delete` (branch-only tag, never a shared tag) |
| GCE | **none** — `gsv-mw-manta*` instances and disks both empty | `gcloud compute instances list --filter name~'gsv-mw'` | n/a (already deleted, verified rc=1 "was not found") |
| Local installs | `bun` 1.4.2 (brew), BetterChromium (`~/.betterwright/chromium`) | `bun --version` | `brew uninstall bun`; `rm -rf ~/.betterwright/chromium` |
| Harness dir | `/Users/markw/Work/manta_tloc_testkit` — **NOT under version control** | `git -C . rev-parse` → `fatal: not a git repository` | n/a; see Open items (backup risk) |

## 3. Resume here (paste-able)

```
 cd /Users/markw/IdeaProjects/gatk-sv && git fetch origin
 git worktree list
 for b in mw_manta_tloc_autoresolve mw_manta_tloc_autoresolve_pr mw_fix_single_sample_blocking; do
   printf '%-32s local=%-9s remote=%s\n' "$b" "$(git rev-parse --short $b)" \
     "$(git ls-remote origin refs/heads/$b | cut -c1-8)"; done
 gh pr view 968 --repo broadinstitute/gatk-sv --json state,baseRefName,changedFiles,additions,deletions
 gh pr view 966 --repo broadinstitute/gatk-sv --json state,mergeable
 cd /Users/markw/Work/manta_tloc_testkit
 /Users/markw/IdeaProjects/gatk-sv-testkit/.venv/bin/python compare_terra.py   # cached VCFs, no cloud
 cd /Users/markw/.pi/agent/skills/terra-monitor
 TERRA_WORKSPACE=broad-firecloud-dsde-methods/GATK-SV-manta-tloc-mw-2026-09-21 \
   python3 scripts/twatch.py list --fresh                                      # expect: no rows
 cd /Users/markw/IdeaProjects/gatk-sv/wt/manta_tloc_autoresolve
 .venv-tloc/bin/python test_single_tloc.py                                     # 19 PASS + 2 NOTE
```

## 4. What "good" looks like on the next check

| check | expected | what a mismatch means |
|---|---|---|
| `gh pr view 968 --json changedFiles,additions,deletions` | `3`, `136`, `11` | branch moved, or #966 merged and GitHub re-based the diff → re-verify touched files before trusting review |
| `gh pr view 966 --json state` | `OPEN` | if `MERGED`, retarget #968 to `main` and re-read the diff (should stay 3 files) |
| `gh pr checks 968` | `no checks reported` | checks appearing means base changed to `main` (docs.yml) — read them, they are free signal |
| `compare_terra.py` (offline) | `CONTROL 156/156 … MATCH`, `EFFECT 156/156`, `total CTX added 27267`, `records +27267`, `PURITY … 0`, `VERDICT: PASS` | `tloc_vcf/` cache or `tloc_compare_sets.json` drifted, or someone reran an arm with different inputs |
| `compare_newbase.py` | 3/3 `PASS`, CTX 172 / 160 / 155 | evidence files renamed (they are `results/ab2nb_*.txt`) |
| `twatch.py list --fresh` | no rows | a submission was started since; check cost before trusting the ledger |
| `gcloud compute instances/disks list --filter name~'gsv-mw'` | both empty | someone provisioned compute; VM-minutes are the unit, never currency |
| `test_single_tloc.py` | 19 `PASS`, 0 `FAIL`, `ALL PASS`, 2 `NOTE` (pysam ≥0.22 `stop<POS` clamp) | resolve semantics regressed, or harness drift |
| `git ls-remote` vs local for both my branches | equal (78dd16ab / ad80380f) | a push was lost or the branch moved elsewhere |
| `freeze_copy.py verify` | `crc32c+size match: 317/317`, `FROZEN VERIFIED` | freeze mutated → re-runs would not be byte-comparable to production |
| Dockstore probe | `201` | `404 Cannot get dockstore://…` = unpublished again, not "missing" |

## 5. Gotchas found (only ones actually hit, with the error text)

1. `POST /methodconfigs` → `400 "The request content was malformed:\nunexpected json type"` from
   **two** distinct causes: `methodVersion: "master"` on an agora ref (must be Int), and real JSON
   list/int values in `inputs` (must be strings). The message names neither field.
2. `404 "Cannot get dockstore://github.com%2Fbroadinstitute%2Fgatk-sv%2FTlocResolveOnly/… from
   method repo."` looked like a missing registration; it meant **unpublished** (§6 correction 1).
3. Rawls **persists a config whose method resolution failed**, so the retry says
   `409 04t-baseline already exists`. `valid` reads `None` forever after, and
   `POST .../methodconfigs/validate` → `405 supported methods: OPTIONS` (unknown path). Free
   pre-flight instead: `miniwdl input_template` locally vs the config's keys.
4. Inverted `outputs` map was accepted at creation and only failed at submission:
   `Validation errors: Invalid outputs: this.tloc_baseline_vcf -> Error while parsing the expr`.
5. `rootEntityType: workspace` is a dead end here: `400 Your method config defines a root entity but
   you haven't passed one to the submission.` → with an entity:
   `500 AttributeEntityReference(workspace,GATK-SV-manta-tloc-mw-2026-09-21) not found` → creating it:
   `400 Entity type workspace is reserved and cannot be overwritten`. Fix: custom type
   `tloc_run/run1`.
6. Both workflows `Succeeded` yet the root entity's attributes stayed **empty** — read Cromwell
   `.../workflows/v1/{id}/metadata` `outputs` instead. Also: call objects report
   **`backendStatus`**; `jobStatus` is `None`.
7. `gcloud storage describe` → `ERROR: (gcloud.storage) Invalid choice: 'describe'`. My checksum
   verifier turned that rc≠0 into `None`, then into "317 objects MISSING" for a copy that had
   completed. `gcloud storage -m cp -I` is also a usage error (rc=2) here. Fix: per-object
   `gcloud storage cp` + JSON API metadata (404 distinguishable from other failures).
8. miniwdl: `invalid choice: 'inputs'` → use `input_template`; and its required-only list is **not**
   "all declared inputs" (`samples_per_shard` has a default, so it looks "unknown to the WDL").
9. `rc=0` printed after a pipe was `head`'s status, which hid the miniwdl usage error above.
10. `bash -n` passes on an empty file and `grep -c '@@'` returns 0 for one, so both guards approved a
    **0-byte startup script** (sed typo `-e "@@MEI@@|…"` missing the `s|`); ~20 idle VM-minutes.
11. The docker A/B runner printed `MISSING LOCAL IMAGE` and carried on; both arms died `rc=1` with
    `records=0` — which reads exactly like "the flag changed nothing". Fresh VMs also lacked
    `/ab/in/cytobands.bed.gz` (`cp: cannot stat`), because the previous VM had leftovers that masked
    the missing download step.
12. Browser tool: `Managed BetterChromium is outdated or its verified installation receipt is
    missing. Run betterwright setup` — needed `brew install bun` first, then `npx -y betterwright setup`.
13. My own `compare_terra.py` asserted "baseline CTX must be 0", mislabelled 7 samples as
    `BAD`, and via an `elif` chain skipped their remaining checks — turning a PASS into a spurious
    `REVIEW`. The data was fine (§6 correction 2).

## 6. Corrections to earlier documents/claims

1. **Was (mine, earlier this session):** "Dockstore won't create the entry from `.dockstore.yml`;
   there is no entry." **Now:** the entry **was** created and was **unpublished**; the owner set the
   default version and published it manually. Unpublished entries are invisible to the anonymous API
   *and* to Rawls, which is the exact 404 I misread. 001's push-triggered registration was never broken.
2. **Was (mine, this session's check logic):** "the baseline arm must have `CTX == 0`." **Now:**
   production itself contains **7 `CTX_*` records with the flag off** (6 `CTX_PP/QQ` + 1 `CTX_PQ/QP`,
   one each on 7 samples) and the baseline arm reproduces those 7 exactly. So `CTX` is not
   exclusively flag-driven; this change accounts for 27,267 of 27,274.
3. **Was (owner hypothesis):** "the `.dockstore.yml` indentation issue fixed recently is why it
   didn't register." **Refuted with evidence:** `fc014646 "fix indentation in dockstore yml (#962)"`
   is already an ancestor of my head, and my entry uses the corrected shape (`name:`/`filters:` at
   4 spaces, `branches:` at 6, list items at 8). The real cause is correction 1.
4. **Was (001 §4, row 4):** "flag-on docker e2e (**NOT run yet**) … production-stack divergence
   (pysam 0.15.4 path untested locally)". **Now run and PASS** — 3 samples on the real
   `sv-pipeline` image (§1b) and 156 samples on Terra (§1c). 001's row and its matching open-item
   checkbox are annotated as superseded.
5. **Was (001):** "**Cost / time / size:** pipeline runtime impact — **not measured**". **Now
   measured:** +2.4% cpu-VM-minutes on the tloc task, +2.8% cohort complex output. 001 annotated.
6. **Was (001 open item):** "Open PR from `mw_manta_tloc_autoresolve` (base `main` …)". **Superseded:**
   PR #968 comes from a *separate* branch (`mw_manta_tloc_autoresolve_pr`) onto #966's branch, with the
   harness and docs deliberately excluded.
7. **Was (001 open item):** "Decide whether `test_single_tloc.py` joins the PR (recommended: yes)".
   **Owner decision this session: no** — the repo has no python test directory to follow; evidence
   goes in the PR description. The file stays untracked in `wt/manta_tloc_autoresolve`.

## 7. Mutations ledger (this session)

| mutation | reversible? | verify / undo |
|---|---|---|
| pushed `bc89a7c9` (Dockstore branch tag) + `78dd16ab` (harness WDL + `.dockstore.yml` entry) to `mw_manta_tloc_autoresolve` | reversible | force-push back to `42b9877c` (pre-harness); Dockstore entry then stale → unpublish |
| created + pushed branch `mw_manta_tloc_autoresolve_pr` (`ad80380f`) | reversible | `git push origin :mw_manta_tloc_autoresolve_pr` |
| opened PR #968 | reversible | `gh pr close 968` |
| committed `002` + annotated `001` on `mw_manta_tloc_autoresolve` | reversible | revert that commit; **never merge docs/handoff to `main`** |
| Terra: created configs `04t-baseline`/`04t-branch`, entity `tloc_run/run1`, 2 submissions; deleted throwaway configs `04t-probe`, `04t-probe2` and the two stale invalid configs | mostly reversible | submissions are terminal (no undo, ~795 cpu-VM-min); configs/entity/bucket deletable |
| copied 317 objects (8.65 GiB) into `gs://fc-5fb1cc39-…/tloc_frozen/` | reversible | `gcloud storage rm -r` that prefix |
| published `markw/gatk-sv/sv-pipeline:mw-manta-tloc-autoresolve-bc89a7` | reversible | delete tag; **no shared/production tag was touched** |
| created then deleted GCE instance `gsv-mw-manta-tloc-autoresolve-bc89a7c9` (+150 GB pd-ssd) | done | verified `describe` rc=1 "was not found", `disks list` empty |
| owner-side (not me): published the Dockstore `TlocResolveOnly` entry + default version | reversible | unpublish |
| local installs: brew `bun`, BetterChromium | reversible | `brew uninstall bun`; `rm -rf ~/.betterwright/chromium` |
| probes created and deleted during debugging | deleted, reported | `GET .../methodconfigs` shows only `04t-baseline`/`04t-branch` |

## 8. Deliverables (files)

| file | what it is |
|---|---|
| gatk-sv `docs/handoff/002_manta_tloc_autoresolve_pr_and_terra_ab.md` | this doc |
| gatk-sv `docs/handoff/manta_tloc_autoresolve.md` (001) | annotated: 3 rows/items marked superseded |
| `…/manta_tloc_testkit/AGENTS.md` | standing notes: Terra/Rawls + Dockstore + gcloud facts with failure text, measured cost figures |
| `…/results/AB_RESULTS_TERRA.md` | Terra setup, result, scope limits, full Terra-API facts list |
| `…/results/AB_RESULTS.md` | docker A/B on both bases + harness-bug postmortems (+ pointer to Terra doc) |
| `…/results/AB_RESULTS_TERRA.json` | machine-readable census + verdict |
| `…/compare_terra.py`, `compare_newbase.py`, `build_tloc_wdl.py`, `terra_write.py`, `terra_probe.py`, `freeze_copy.py` | tooling; `terra_write.py` refuses writes without `--confirm` and refuses the baseline workspace |
| `…/wdl/TlocResolveOnly.wdl` | the trimmed harness (generated; also on the branch) |
| `…/freeze_verified.json`, `freeze_manifest.json`, `tloc_compare_sets.json` | reproducibility records |
| `…/tloc_vcf/{prod,baseline,branch}/` | 156×3 downloaded complex VCFs (~60 MiB) so `compare_terra.py` runs offline |

## Open items / next steps

- [ ] Retarget #968 to `main` when #966 merges; re-read `changedFiles` (expect 3) and read the CI that
      finally appears (`docs.yml`).
- [ ] Ask the #966 owner to look at #968 (stacked on their branch) — nothing in this session pinged anyone.
- [ ] Confirm with the eph tloc workflow owner: the `EVIDENCE=PE` gate is inert inside
      `mantatloccheck.sh` (it stamps PE on every record, incl. SR-only), so `manta_tloc` evidence
      volume now equals all interchromosomal manta BNDs. Needs their ack (carried from 001, **not**
      re-examined this session).
- [ ] Follow-up ticket before enabling the flag anywhere else:
      `src/sv-pipeline/scripts/single_sample/update_variant_representations.py` CTX→4×BND expansion has
      a wrong REF anchor for M4 (~line 115) and a self-referential `MATEID` (~line 157) — dormant
      behind the flag, reproduced by the blast-radius reviewer in 001, **not** re-verified this session.
- [ ] Decide the leftovers: keep or delete `wdl/TlocResolveOnly.wdl` + the Dockstore entry (deleting
      the file breaks the version the two working configs point at); keep or delete the Terra workspace
      + `tloc_frozen/`; delete or commit `test_single_tloc.py`.
- [ ] **Back up `/Users/markw/Work/manta_tloc_testkit` — it is not a git repo.** All tooling, results,
      freeze manifests and the 3-way VCF cache exist only on this laptop. Also note `ab_ed25519` (a
      private key used for the A/B VM) sits in that directory; move it somewhere access-controlled and
      delete it before any sync-to-cloud/commit of the directory.
- [ ] Optional stronger evidence, not requested: run the **full** step-04 branch-vs-baseline
      (31,206 VM-minutes/arm) to catch interactions with the tasks the trim dropped. Not needed for
      this flag's blast radius; the trim's soundness argument is in §1c.
- [ ] Not re-checked this session: `pi-subagents` background/async lanes wanting `subagent-runner.ts`
      on Node 26 (001 §5 gotcha 1). Foreground lanes worked fine throughout.
