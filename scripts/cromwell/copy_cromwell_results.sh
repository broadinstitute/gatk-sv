#!/bin/bash
# USAGE: copy_cromwell_results.sh WORKFLOW_ID OUTPUT_FOLDER_NAME [DESTINATION_PATH]
# Copy all output files from cromwell workflow to flat folder at
# destination path.
# -If cromshell is not located at ${HOME}/Documents/cromshell/cromwell
#  then define the ENV variable CROMSHELL with the appropriate path.
# -If DESITNATION_PATH is not specified, it is assumed to be
#   "gs://broad-sv-dev-data/caller_benchmark/callers"
set -Eeuo pipefail

WORKFLOW_ID=$1
CALLER_DIR=$2
DESTINATION_PATH=${DESTINATION_PATH:-"gs://broad-sv-dev-data/caller_benchmark/callers"}
DESTINATION_PATH=${3:-DESTINATION_PATH}

GCS_OUTPUT_DIR="${DESTINATION_PATH%/}/${CALLER_DIR}"

# find the outputs
CROMSHELL=${CROMSHELL:-${HOME}/Documents/cromshell/cromwell}
METADATA=$($CROMSHELL -t 200 metadata $WORKFLOW_ID)
OUTPUTS=$(echo "$METADATA" | jq '.outputs | .[]' | sed -e 's/ //g' -e 's/^"//' -e 's/,$//' -e 's/"$//' | grep -v '\]\|\[' | sort)
# copy the outputs to the requested GCS location
echo "$OUTPUTS" | gsutil -m cp -I "$GCS_OUTPUT_DIR"
