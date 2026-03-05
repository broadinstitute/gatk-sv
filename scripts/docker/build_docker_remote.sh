#!/usr/bin/env bash
# build_docker_remote.sh
#
# Runs the full GATK-SV Docker build workflow on a remote Linux VM.
# Usage: ./scripts/docker/build_docker_remote.sh [options]
#
# This script must be run from the root of the gatk-sv repo on your local machine.
# The Broad VPN must be active.

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────────
DEFAULT_SSH_HOST="your_host_name"
DEFAULT_DOCKER_REPO="us.gcr.io/your-repo-name"
DEFAULT_REMOTE_REPO_DIR="~/gatk-sv"
DEFAULT_BASE_GIT_COMMIT="origin/main"

# ── Usage ─────────────────────────────────────────────────────────────────────────
usage() {
  cat <<HELP
Usage: $(basename "$0") [options]

Commits and pushes the current branch, SSHes to a remote build VM, runs the
Docker build, updates inputs/values/dockers.json, and pulls the result locally.

Options:
  -H, --ssh-host HOST         SSH hostname of the build VM  (default: $DEFAULT_SSH_HOST)
  -r, --docker-repo REPO      GCR/Docker repo prefix        (default: $DEFAULT_DOCKER_REPO)
  -d, --remote-dir DIR        Path to gatk-sv repo on VM    (default: $DEFAULT_REMOTE_REPO_DIR)
  -b, --base-commit REF       Base git ref for change detection (default: $DEFAULT_BASE_GIT_COMMIT)
  (git on the remote VM is assumed to be already configured)
  -m, --commit-message MSG    Local commit message if there are uncommitted changes
                              (default: "Trigger docker build")
  -n, --no-local-commit       Skip the local commit/push step (use if already pushed)
  -D, --dry-run               Dry run: verify connectivity and show commands (no actions)
  -h, --help                  Show this help message

Examples:
  # Typical usage – all defaults:
  ./scripts/docker/build_docker_remote.sh

  # Custom VM and repo:
  ./scripts/docker/build_docker_remote.sh \\
      --ssh-host my-vm.example.com \\
      --docker-repo us.gcr.io/my-project/my-repo

HELP
  exit 0
}

# ── Argument parsing ──────────────────────────────────────────────────────────────
SSH_HOST="$DEFAULT_SSH_HOST"
DOCKER_REPO="$DEFAULT_DOCKER_REPO"
REMOTE_DIR="$DEFAULT_REMOTE_REPO_DIR"
BASE_GIT_COMMIT="$DEFAULT_BASE_GIT_COMMIT"
# git identity is assumed configured on the remote VM
COMMIT_MESSAGE="Trigger docker build"
DO_LOCAL_COMMIT=true
DRY_RUN=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    -H|--ssh-host)        SSH_HOST="$2";        shift 2 ;;
    -r|--docker-repo)     DOCKER_REPO="$2";     shift 2 ;;
    -d|--remote-dir)      REMOTE_DIR="$2";      shift 2 ;;
    -b|--base-commit)     BASE_GIT_COMMIT="$2"; shift 2 ;;
    # remote git identity flags removed; assume VM already configured
    -m|--commit-message)  COMMIT_MESSAGE="$2";  shift 2 ;;
    -n|--no-local-commit) DO_LOCAL_COMMIT=false; shift ;;
    -D|--dry-run)         DRY_RUN=true;         shift ;;
    -h|--help)            usage ;;
    *) echo "Unknown option: $1" >&2; usage ;;
  esac
done

# NOTE: we assume git on the remote VM is already configured with a user.name/email

# ── Safety: refuse to build from main ────────────────────────────────────────────
BRANCH="$(git rev-parse --abbrev-ref HEAD)"
if [[ "$BRANCH" == "main" ]]; then
  echo "ERROR: You are on the 'main' branch. Switch to a feature branch before building." >&2
  exit 1
fi

# ── Step 1: Commit and push local changes ─────────────────────────────────────────
  if [[ "$DO_LOCAL_COMMIT" == true ]]; then
  echo "==> [1/5] Committing and pushing branch '$BRANCH'..."
  if [[ "$DRY_RUN" == true ]]; then
    if [[ -n "$(git status --porcelain)" ]]; then
      echo "    DRY RUN: would git add, commit '$COMMIT_MESSAGE', and push to origin/$BRANCH"
    else
      echo "    DRY RUN: no local changes to commit. Would still push origin/$BRANCH if needed."
    fi
    echo "    DRY RUN: skipping actual git push."
  else
    if [[ -n "$(git status --porcelain)" ]]; then
      git add -A
      git commit -m "$COMMIT_MESSAGE"
      echo "    Committed local changes."
    else
      echo "    No local changes to commit."
    fi
    git push --force origin "$BRANCH"
    echo "    Force-pushed to origin/$BRANCH."
  fi
else
  echo "==> [1/5] Skipping local commit/push (--no-local-commit)."
fi

# ── Compute image tag: branch-name (underscores→hyphens) + first-6-sha ───────────
SHA="$(git rev-parse --short=6 HEAD)"
IMAGE_TAG="$(echo "$BRANCH" | tr '_' '-')-${SHA}"
echo "    Image tag: $IMAGE_TAG"

# ── Step 2: Verify SSH connectivity ──────────────────────────────────────────────
echo "==> [2/5] Verifying SSH connectivity to '$SSH_HOST'..."
if ! ssh -o BatchMode=yes -o ConnectTimeout=10 "$SSH_HOST" true 2>/dev/null; then
  echo "ERROR: Cannot reach '$SSH_HOST'. Make sure you are on the Broad VPN." >&2
  exit 1
fi
echo "    SSH OK."

# If dry-run: do lightweight remote checks and print the build command, then exit
if [[ "$DRY_RUN" == true ]]; then
  echo "==> [DRY RUN] Performing remote sanity checks (no build or commits will run)..."
  ssh -o BatchMode=yes -o ConnectTimeout=10 "$SSH_HOST" bash -s "$REMOTE_DIR" "$BRANCH" "$IMAGE_TAG" "$DOCKER_REPO" <<'REMOTE_DRY'
set -euo pipefail
REMOTE_DIR="$1"; BRANCH="$2"; IMAGE_TAG="$3"; DOCKER_REPO="$4"
eval REMOTE_DIR="$REMOTE_DIR"
if [[ ! -d "$REMOTE_DIR" ]]; then
  echo "  Remote dir $REMOTE_DIR not found"
  exit 0
fi
cd "$REMOTE_DIR"
echo "  Remote repo: $(pwd)"
echo "  Remote: git rev-parse --short origin/$BRANCH -> $(git rev-parse --short origin/$BRANCH 2>/dev/null || echo 'not found')"
echo "  Remote: python3 -> $(which python3 2>/dev/null || echo 'not found')"
REMOTE_DRY

  echo "DRY RUN: Would run on remote:"
  echo "  cd $REMOTE_DIR"
  echo "  git checkout origin/$BRANCH"
  echo "  cd scripts/docker"
  echo "  python3 build_docker.py --skip-cleanup --docker-repo \"$DOCKER_REPO\" --image-tag \"$IMAGE_TAG\" --base-git-commit \"$BASE_GIT_COMMIT\" --disable-git-protect"
  echo "DRY RUN complete."; exit 0
fi

# ── Step 3: Run the Docker build on the VM ────────────────────────────────────────
echo "==> [3/5] Running Docker build on '$SSH_HOST' (this may take several minutes)..."
# Pass args as positional parameters so the heredoc is self-contained
if [[ "$DRY_RUN" != true ]]; then
ssh "$SSH_HOST" bash -s "$REMOTE_DIR" "$BRANCH" "$IMAGE_TAG" "$DOCKER_REPO" "$BASE_GIT_COMMIT" << 'REMOTE_BUILD'
set -euo pipefail
REMOTE_DIR="$1"; BRANCH="$2"; IMAGE_TAG="$3"; DOCKER_REPO="$4"; BASE_GIT_COMMIT="$5"

# expand ~ in REMOTE_DIR
eval REMOTE_DIR="$REMOTE_DIR"
cd "$REMOTE_DIR"

echo "  Cleaning untracked inputs/values/ files that would block checkout..."
git clean -fd inputs/values/ 2>/dev/null || true

echo "  Updating main..."
git checkout main && git pull --quiet
git fetch origin

echo "  Checking out origin/$BRANCH (detached HEAD)..."
git checkout "origin/$BRANCH"
echo "  HEAD: $(git rev-parse --short HEAD)"

echo "  Running build_docker.py..."
cd scripts/docker
python3 build_docker.py \
  --skip-cleanup \
  --docker-repo "$DOCKER_REPO" \
  --image-tag "$IMAGE_TAG" \
  --base-git-commit "$BASE_GIT_COMMIT" \
  --disable-git-protect
echo "  Build complete."
REMOTE_BUILD
else
  echo "==> [3/5] Skipped Docker build due to --dry-run."
fi

# ── Step 4: Commit and push updated dockers.json from the VM ──────────────────────
if [[ "$DRY_RUN" == true ]]; then
  echo "==> [4/5] DRY RUN: Would commit updated dockers.json on VM (skipped)."
else
  echo "==> [4/5] Committing updated dockers.json on VM..."
  ssh "$SSH_HOST" bash -s "$REMOTE_DIR" "$BRANCH" << 'REMOTE_COMMIT'
set -euo pipefail
REMOTE_DIR="$1"; BRANCH="$2"

eval REMOTE_DIR="$REMOTE_DIR"
cd "$REMOTE_DIR"

DOCKERS_JSON="inputs/values/dockers.json"
if [[ ! -f "$DOCKERS_JSON" ]]; then
  echo "  WARNING: $DOCKERS_JSON not found; skipping commit." >&2
  exit 0
fi

# Check if the file was actually modified compared to the checked-out commit
if git diff --quiet HEAD -- "$DOCKERS_JSON"; then
  echo "  dockers.json unchanged; nothing to commit."
  exit 0
fi

git add "$DOCKERS_JSON"
git commit -m "Update docker tags after build"
git push origin "HEAD:refs/heads/$BRANCH"
echo "  Pushed updated dockers.json to origin/$BRANCH."
REMOTE_COMMIT
fi

# ── Step 5: Pull updated branch locally ───────────────────────────────────────────
if [[ "$DRY_RUN" == true ]]; then
  echo "==> [5/5] DRY RUN: Would pull updated branch locally (skipped)."
else
  echo "==> [5/5] Pulling updated branch locally..."
  git pull
fi
echo ""
echo "✅ Docker build complete!"
echo "   Branch:     $BRANCH"
echo "   Image tag:  $IMAGE_TAG"
echo "   Repository: $DOCKER_REPO"
