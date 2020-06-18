#!/bin/bash

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Intelligently shards a VCF prior to complex resolution (for parallelization)

# Subsetted to first half, to just output lists of which variants should go in each shard

set -Eeu -o pipefail

# ARGS defaults and hard-coded values
# generous 1kb clustering, 10% RO clustering
DEFAULT_DIST=1000
DEFAULT_RECIP=0.1
DEFAULT_MIN_LINES_PER_SHARD=10
DEFAULT_MAX_SHARDS=100
DEFAULT_NONCLUSTER_SHARDS=30
DEFAULT_PREFIX="vcf_shard"
DEFAULT_BREAKPOINT_PADDING=5000
DEFAULT_IGNORE_SV_TYPES=false
DEFAULT_ADD_SINGLE_REC=false
DEFAULT_SHARD_LARGE_CLUSTERS=true
SCRIPT_NAME=$(basename "${BASH_SOURCE[0]}")

BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $BIN/shardVCF_backend_part1.sh
