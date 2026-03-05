#!/bin/bash

# Script to run gd_cnv_pyro.py on MARK-RYZEN9 remote machine via SCP and SSH
# Only transfers files that have changed (compares hashes)
# Usage: ./run_gd_cnv_pyro_remote.sh <work_dir>
# Expected structure: <work_dir>/input/ and <work_dir>/model/

set -e

REMOTE_USER="mwalk"
REMOTE_HOST="MARK-RYZEN9"
REMOTE_BASE="/mnt/c/rsync"
LOCAL_SCRIPT_DIR="/Users/markw/IdeaProjects/gatk-sv/src/sv-pipeline/scripts"
SCRIPT_NAME="gd_cnv_pyro.py"

# Parse arguments
WORK_DIR="${1}"

if [[ -z "$WORK_DIR" ]]; then
    echo "Usage: $0 <work_dir>"
    echo "Example: $0 /Users/markw/Work/talkowski/sv-pipe-testing/mw_gd/gd_pyro"
    exit 1
fi

if [[ ! -d "$WORK_DIR" ]]; then
    echo "Error: Work directory does not exist: $WORK_DIR"
    exit 1
fi

if [[ ! -d "$WORK_DIR/input" ]]; then
    echo "Error: $WORK_DIR/input directory does not exist"
    exit 1
fi

WORK_BASENAME=$(basename "$WORK_DIR")
REMOTE_WORK_DIR="$REMOTE_BASE/$WORK_BASENAME"

echo "=========================================="
echo "Running gd_cnv_pyro.py on MARK-RYZEN9"
echo "=========================================="
echo "Work directory: $WORK_DIR"
echo ""

# Function to transfer file only if hash differs
transfer_if_changed() {
    local local_file="$1"
    local remote_file="$2"
    local file_name=$(basename "$local_file")
    local large_file_threshold=$((10 * 1024 * 1024))  # 10MB in bytes

    if [[ ! -f "$local_file" ]]; then
        echo "  ✗ File not found: $local_file"
        return 1
    fi

    local local_size=$(stat -f%z "$local_file")
    local remote_size=$(echo "stat -c%s '$remote_file' 2>/dev/null || echo 0" | ssh "$REMOTE_USER@$REMOTE_HOST" bash || echo "0")

    if [[ "$local_size" -gt "$large_file_threshold" ]]; then
        # Large file: compare by size only
        if [[ "$local_size" == "$remote_size" ]]; then
            echo "  ✓ $file_name (same size ${local_size} bytes, skipping)"
            return 0
        fi
    else
        # Small file: compare by MD5 hash
        local local_hash=$(md5sum "$local_file" | awk '{print $1}')
        local remote_hash=$(echo "md5sum '$remote_file' 2>/dev/null | awk '{print \$1}'" | ssh "$REMOTE_USER@$REMOTE_HOST" bash || echo "")
        if [[ "$local_hash" == "$remote_hash" ]]; then
            echo "  ✓ $file_name (unchanged, skipping)"
            return 0
        fi
    fi

    if [[ "$remote_size" == "0" ]]; then
        echo "  → $file_name (new file, transferring)"
    else
        echo "  ↻ $file_name (changed, updating)"
    fi

    # Pipe command then file content together to ssh bash stdin.
    # base64 avoids binary issues; the first line is the bash command,
    # remainder is stdin for that command (decoded back to binary).
    { printf "base64 -d > '%s'\n" "$remote_file"; base64 < "$local_file"; } \
        | ssh "$REMOTE_USER@$REMOTE_HOST" bash
}

# Fetch a single file from remote to local via SSH+bash (avoids SCP /mnt/c path issues)
fetch_from_remote() {
    local remote_file="$1"
    local local_file="$2"
    echo "base64 '$remote_file'" | ssh "$REMOTE_USER@$REMOTE_HOST" bash | base64 -d > "$local_file"
}

# Step 1: Create remote directory and sync input files
echo "[1/4] Setting up remote directory..."
echo "mkdir -p $REMOTE_BASE/$WORK_BASENAME/input $REMOTE_BASE/$WORK_BASENAME/model" | ssh "$REMOTE_USER@$REMOTE_HOST" bash

echo "[2/4] Transferring input files (checking hashes)..."
transfer_if_changed "$WORK_DIR/input/rebatch_1.10kbp_bins.rd.txt.gz" "$REMOTE_WORK_DIR/input/rebatch_1.10kbp_bins.rd.txt.gz"
transfer_if_changed "$WORK_DIR/input/GenomicDisorderRegions_hg38_2025-12-05.with_bp.tsv" "$REMOTE_WORK_DIR/input/GenomicDisorderRegions_hg38_2025-12-05.with_bp.tsv"
transfer_if_changed "$WORK_DIR/input/hg38_SD.bed.gz" "$REMOTE_WORK_DIR/input/hg38_SD.bed.gz"
transfer_if_changed "$WORK_DIR/input/rebatch_1.RD.txt.gz" "$REMOTE_WORK_DIR/input/rebatch_1.RD.txt.gz"

# Step 2: Transfer Python script
echo "[3/4] Transferring Python script..."
transfer_if_changed "$LOCAL_SCRIPT_DIR/$SCRIPT_NAME" "$REMOTE_WORK_DIR/$SCRIPT_NAME"

# Step 3: Run the script on remote
echo "[4/4] Running script on remote machine..."
ssh "$REMOTE_USER@$REMOTE_HOST" bash << EOF
    set -e
    cd $REMOTE_WORK_DIR
    python3 $REMOTE_WORK_DIR/$SCRIPT_NAME \
        -i input/rebatch_1.10kbp_bins.rd.txt.gz \
        -g input/GenomicDisorderRegions_hg38_2025-12-05.with_bp.tsv \
        -s input/hg38_SD.bed.gz \
        -o model/model \
        --max-iter 3000 \
        --min-delta 1000 \
        --n-discrete-samples 1000 \
        --high-res-counts input/rebatch_1.RD.txt.gz \
        --device cuda \
        --var-sample 0.0001
EOF

# Step 4: Transfer results back (check for changed files)
echo ""
echo "Syncing results back to local machine (checking hashes)..."
mkdir -p "$WORK_DIR/model"

# Get list of output files from remote (top-level only, no subdirs)
REMOTE_FILES=$(echo "find $REMOTE_WORK_DIR/model -maxdepth 1 -type f 2>/dev/null" | ssh "$REMOTE_USER@$REMOTE_HOST" bash || echo "")

# Batch fetch all remote sizes and hashes in a single SSH call
REMOTE_STATS=$(echo "for f in $REMOTE_FILES; do echo \"\$f \$(stat -c%s \"\$f\" 2>/dev/null || echo 0) \$(md5sum \"\$f\" 2>/dev/null | awk '{print \$1}')\"; done" | ssh "$REMOTE_USER@$REMOTE_HOST" bash || echo "")

while IFS=' ' read -r remote_file remote_size remote_hash; do
    [[ -z "$remote_file" ]] && continue
    filename=$(basename "$remote_file")
    local_file="$WORK_DIR/model/$filename"
    large_file_threshold=$((10 * 1024 * 1024))
    local_size=$(stat -f%z "$local_file" 2>/dev/null || echo "0")

    if [[ "$remote_size" -gt "$large_file_threshold" ]]; then
        if [[ "$local_size" == "$remote_size" ]]; then
            echo "  ✓ $filename (same size, skipping)"
            continue
        fi
    else
        local_hash=$(md5sum "$local_file" 2>/dev/null | awk '{print $1}' || echo "")
        if [[ "$local_hash" == "$remote_hash" ]]; then
            echo "  ✓ $filename (unchanged, skipping)"
            continue
        fi
    fi

    if [[ "$local_size" == "0" ]]; then
        echo "  → $filename (new file, transferring)"
    else
        echo "  ↻ $filename (updated, fetching)"
    fi
    fetch_from_remote "$remote_file" "$local_file"
done <<< "$REMOTE_STATS"

echo ""
echo "=========================================="
echo "✓ Complete! Results are in: $WORK_DIR/model"
echo "=========================================="
