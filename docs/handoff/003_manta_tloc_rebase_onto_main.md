# 003 — Both manta-tloc branches rebased onto `main` after #966 squash-merged

**Handoff written:** `2026-10-08T20:18Z`. Clock agreement verified, not assumed:
local `date -u` = `github` (`api.github.com` HTTP Date) = `Thu, 08 Oct 2026 20:18:06 GMT`, 0 s drift.
No Terra leg in this session, so no third clock to reconcile.

Predecessor: `002` = `docs/handoff/002_manta_tloc_autoresolve_pr_and_terra_ab.md`. Read it for the
validation evidence (docker A/B, Terra head-to-head) — none of that evidence is invalidated here,
because the tloc code is byte-identical across the rebase (§3). **`002` §4 has two now-stale rows;
§4 of this doc supersedes them.**

**This doc exists mainly as a hazard warning:** the branches moved and history was rewritten. If you
find a worktree whose `HEAD` disagrees with `origin`, do **not** "restore" it with
`git reset --hard origin/<branch>` — that discards the rebase. The pre-rebase state is kept as
`tmp/pre-rebase-2026-10-08-*` (§1), not on the remote.

---

## 1. Coordinates

| thing | before | after | how to verify |
|---|---|---|---|
| `main` | `e1909d2f` (old merge base) | `743dd9d4` | `git rev-parse origin/main` |
| The bug fixes | — | `a55c498f` = PR **#966 squash-merged** `2026-10-08T15:17:17Z`, plus `339d1e43` (#959 DRAGEN CNV standardize) + two "Update docker images list" bumps | `gh pr view 966 --json state,mergedAt,mergeCommit` |
| `mw_manta_tloc_autoresolve` | `d47cefe0` | **`eeb913b0`** (+ this doc: 1 more) | `git rev-parse mw_manta_tloc_autoresolve` |
| `mw_manta_tloc_autoresolve_pr` (PR #968 head) | `ad80380f` | **`0ed9e243`** | `git rev-parse mw_manta_tloc_autoresolve_pr` |
| Pre-rebase safety refs (local only) | — | `tmp/pre-rebase-2026-10-08-autoresolve` → `d47cefe0`, `tmp/pre-rebase-2026-10-08-pr` → `ad80380f` | `git rev-parse tmp/pre-rebase-2026-10-08-autoresolve` |
| Commits **dropped** as already-in-main | `144a3fae` ConcatBaf, `8c70779b` pkg_resources, `9f2ebe25` single-sample workspace table, `4419315c` Dockstore publish | — (their content is main's, see §2) | `git range-diff 4419315c..d47cefe0 origin/main..eeb913b0` shows 7 `=`, no gaps |
| PR #968 | base auto-retargeted `mw_fix_single_sample_blocking`→`main` when #966 merged; `mergeable: CONFLICTING`, 12 files `+204/-35` | base `main`, **3 files `+136/-11`** | `gh pr view 968 --json baseRefName,changedFiles,additions,deletions,mergeable` |
| Remote refs | `d47cefe0` / `ad80380f` | pushed with `--force-with-lease` bound to exactly those SHAs | `git ls-remote origin refs/heads/mw_manta_tloc_autoresolve refs/heads/mw_manta_tloc_autoresolve_pr` |
| Worktrees | `wt/manta_tloc_autoresolve`, `wt/tloc-pr` (unchanged paths) | same two, `HEAD` == branch, tracked-clean | `git worktree list` |
| Scratch worktrees used for the rehearsal | `wt/rebase-scratch-autoresolve`, `wt/rebase-scratch-pr` | **removed** | `git worktree list` (absent) |
| Repo-local ignore rule added | `wt/` was untracked-but-not-ignored (`?? wt/` in every agent's `git status`) | `.git/info/exclude` now has `wt/` | `git check-ignore -v wt` |

Nothing outside this repo was created or changed: no images, no Terra objects, no Dockstore action, no
compute. The only remote effect is the two branch ref updates + the #968 diff refresh.

## 2. Why a plain `git rebase origin/main` would have been wrong here

`git cherry -v origin/main mw_manta_tloc_autoresolve` lists **all 11** branch commits as `+`
("not upstream"), including the four that #966 delivered. `git cherry` compares **patch-ids**, and a
squash-merge rewrites the patches, so deduplication silently fails. A plain rebase would have replayed
those four onto main's own copy of them — which is precisely what made GitHub report
`CONFLICTING` for #968 after the base retargeted (`8c70779b`'s `Dockerfile`/`rdtest2vcf.py`/
`rdtest.py`/`vcfcluster.py`/`__init__.py` hunks vs main's merged versions of the same edits).

The command that was used, on both branches, cuts at the last already-merged commit:

```
git rebase --onto origin/main 4419315c
```

Zero conflicts on both branches — main's 30 changed files and this branch's 7 do not intersect at all.
The pre-merge version of the `pkg_resources` change was also *superseded*, not merely duplicated:
main's squash of #966 is 8/0 on `svtk/__init__.py` where the branch had 21/2, so the rebased tree
adopts main's stricter version and the branch's older copy disappears. That is the intended outcome,
not a lost edit — verify with `git show 8c70779b --stat`.

`.dockstore.yml` legitimately loses 5 lines in the rebase: the five `branches: -
mw_fix_single_sample_blocking` filters added by `4419315c`. That branch is merged and gone, so the
filters matched nothing. The `+13` remaining against main is only `TlocResolveOnly` +
`GatherBatchEvidence` on `mw_manta_tloc_autoresolve`.

## 3. Evidence that the rebase preserved the work

| check | result |
|---|---|
| `git range-diff 4419315c..d47cefe0 origin/main..eeb913b0` | 7 pairs, **all `=`** — every kept commit replayed unchanged, 1:1, no fixups |
| `git range-diff 4419315c..ad80380f origin/main..0ed9e243` | 3 pairs, all `=` |
| tloc delta vs old base vs new base, `diff <(git diff e1909d2f <old-tip> -- <3 code files>) <(git diff origin/main <new-tip> -- <3 code files>)` | **identical** (280 diff lines) on both branches |
| `git diff --numstat origin/main mw_manta_tloc_autoresolve_pr` | `3` files, `+136/-11` — exactly `002` §1a's PR shape |
| `git hash-object` of the 3 code files in both working trees vs `git rev-parse <old-tip>:<path>` | all 3 blob SHAs match (`mantatloccheck.sh` `853c9ba2`, `resolve.py` `1ba52125`, `complex_sv.py` `4362cf22`) |
| `git diff --name-only e1909d2f origin/main` ∩ branch's changed files | **empty** — no file touched by both sides, which is why no conflict was possible |
| `.venv-tloc/bin/python test_single_tloc.py` (imports svtk from *this worktree*, verified via `svtk.cxsv.complex_sv.__file__`) | **19 `PASS`, 0 `FAIL`, `ALL PASS`, 2 `NOTE`** — matches `002` §4's expectation, so resolve semantics did not regress against main's svtk changes |
| `python -m py_compile resolve.py complex_sv.py` / `bash -n mantatloccheck.sh` | rc=0 / rc=0 |

## 4. What "good" looks like on the next check (supersedes the two stale `002` §4 rows)

| check | expected | what a mismatch means |
|---|---|---|
| `git rev-list --left-right --count origin/main...<branch>` | `0 7` autoresolve, `0 3` pr | someone merged or rebased again; re-read before assuming conflict markers are yours |
| `git ls-remote` vs local for both branches | equal to `HEAD` | a push was lost, or an agent reset one side — see the hazard note above |
| `gh pr view 968 --json baseRefName,changedFiles,additions,deletions` | `main`, `3`, `136`, `11` | if `changedFiles > 3`, the PR branch picked up non-tloc commits again |
| `gh pr view 968 --json mergeable` | `MERGEABLE` | still `CONFLICTING` ⇒ main moved again; rebase the same way (`--onto origin/main <last-merged-commit>`, never a plain rebase after a squash-merge) |
| `gh pr view 966 --json state` | `MERGED` | if reopened, the dropped-commit analysis in §2 must be redone |
| `gh pr checks 968` | checks now **may** exist — `docs.yml` triggers on `pull_request` to `main`, and the base is now `main` (`002` §1a's "no checks reported" was true only while the base was #966's branch) | read them; they are free signal, not a regression |
| `test_single_tloc.py` | `19 PASS`, `0 FAIL`, `ALL PASS`, `2 NOTE` | resolve semantics regressed, or harness drift |
| `git status` in either worktree | tracked-clean; `?? test_single_tloc.py` expected (untracked, not ignored — keep it, `002` §3 runs it) | an agent left edits here; do not `reset --hard` through them |

## 5. Resume here (paste-able)

```
cd /Users/markw/IdeaProjects/gatk-sv && git fetch origin --prune
git worktree list
for b in mw_manta_tloc_autoresolve mw_manta_tloc_autoresolve_pr; do
  printf '%-32s local=%-9s remote=%s\n' "$b" "$(git rev-parse --short $b)" \
    "$(git ls-remote origin refs/heads/$b | cut -c1-9)"; done
git rev-list --left-right --count origin/main...mw_manta_tloc_autoresolve_pr   # expect: 0 3
gh pr view 968 --repo broadinstitute/gatk-sv \
  --json baseRefName,changedFiles,additions,deletions,mergeable                # expect: main 3 136 11 MERGEABLE
gh pr checks 968 --repo broadinstitute/gatk-sv
git -c core.pager=cat range-diff 4419315c..tmp/pre-rebase-2026-10-08-autoresolve \
    origin/main..mw_manta_tloc_autoresolve                                     # expect: 7 x '='
cd /Users/markw/IdeaProjects/gatk-sv/wt/manta_tloc_autoresolve
.venv-tloc/bin/python test_single_tloc.py                                      # 19 PASS + 2 NOTE
```

The validation legs in `002` §1b/§1c (`compare_terra.py`, `compare_newbase.py`, `twatch.py`,
`freeze_copy.py verify`) were **not rerun** in this session — they are offline/cached checks against
frozen artefacts, and the byte-identity evidence in §3 says re-running them would reproduce the same
numbers. If you want a fresh head-to-head on post-#966 images, that is a new A/B, not a rerun: the
production baseline image in `002` §1b predates main's `a55c498f`, so the CONTROL arm would need to be
re-frozen.

## 6. Gotchas found (only ones actually hit, with the error text)

1. **`git cherry` is blind to squash-merges.** It printed `+` for commits that are in main. Trust
   `git diff --stat origin/main <branch>` (does the tree still differ?) over patch-id presence, and
   check `gh pr view <n> --json state,mergedAt` before deciding a commit is "mine to replay".
2. **A merged PR base silently retargets its children.** `002` §1a recorded #968's base as
   `mw_fix_single_sample_blocking`; after #966 merged and its branch vanished, GitHub retargeted #968
   to `main` and its diff grew from 3 files to 12 while turning `CONFLICTING`. Always re-read
   `baseRefName` after a stack parent merges — do not trust a remembered base.
3. **`git worktree add -b` / `git branch -f` refuse on a branch another worktree owns**, and
   `git update-ref` happily *does* it, which leaves the owning worktree's index+worktree describing
   the old commits (so `git status` there reports hundreds of phantom edits). Rebase *inside* the
   owning worktree, or `git reset --hard <verified-sha>` inside it — never move the ref from outside.
4. **`wt/` was not ignored** in this repo (`.gitignore` has no `wt/` entry, so `git status` showed
   `?? wt/` with every agent's worktrees inside it — one `git add -A` from disaster). Fixed locally in
   `.git/info/exclude`, which is deliberately *not* committed to keep this session off the shared
   `main` checkout other agents are sitting on. The tracked `.gitignore` still needs an owner commit.
5. **Parallelism is live, not hypothetical.** During this session `wt/trio-denovo` moved
   `34dd0010`→`e9d7ba73` and a new `wt/gatk-sv-rebase` appeared, both while I worked. Everything here
   ran in `wt/*` worktrees; the shared checkout at `IdeaProjects/gatk-sv` was never written to.
