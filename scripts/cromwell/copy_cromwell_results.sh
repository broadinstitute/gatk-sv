#!/bin/bash
# USAGE: copy_cromwell_results.sh WORKFLOW_ID [DESTINATION_PATH] [SUB_FOLDER]
# Copy all output files from cromwell workflow to flat folder at
# destination path.
# -If cromshell is not available on the path as "cromshell"
#  then define the ENV variable CROMSHELL with the appropriate path.
# -If DESTINATION_PATH is not specified, it is assumed to be
#   "gs://broad-sv-dev-data/caller_benchmark/callers"
# -if SUB_FOLDER is not defined, copy directly to DESTINATION_PATH
set -Eeuo pipefail

WORKFLOW_ID=$1
DESTINATION_PATH=${DESTINATION_PATH:-"gs://broad-sv-dev-data/caller_benchmark/callers"}
DESTINATION_PATH=${2:-$DESTINATION_PATH}
SUB_FOLDER=${3:-""}

if [ -n "$SUB_FOLDER" ]; then
  GCS_OUTPUT_DIR="${DESTINATION_PATH%/}/${SUB_FOLDER%/}"
else
  GCS_OUTPUT_DIR="${DESTINATION_PATH%/}"
fi

# find the outputs
CROMSHELL=${CROMSHELL:-cromshell}
$CROMSHELL -t 200 metadata $WORKFLOW_ID > ./metadata.json
OUTPUTS=$(jq '.outputs | .[]' ./metadata.json \
          | sed -e 's/ //g' -e 's/^"//' -e 's/,$//' -e 's/"$//' \
          | grep -v '\]\|\[' | grep -v "null" \
          | sort)
# copy the outputs to the requested GCS location
printf "$OUTPUTS\n./metadata.json\n" | gsutil -m cp -I "$GCS_OUTPUT_DIR/"
