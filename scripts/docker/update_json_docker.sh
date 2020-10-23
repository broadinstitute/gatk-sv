#!/bin/bash


####################################
####################################
# Updates GATKSV json file dockers
####################################
####################################

set -euo pipefail

#################################################
# Utility functions
#################################################

function usage() {
  printf "Usage: \n \
    %s -d <GATKSV_BASE_DIR> -r <ROOT> -n <NAME> -t <TAG> \n \
    <GATKSV_BASE_DIR> \t path to gatk-sv-v1 base directory \n \
    <ROOT> \t docker root/user (e.g. \"gatksv\") \n \
    <NAME> \t docker image name (e.g. \"sv-pipeline\", \"ALL\") \n \
    <TAG> \t docker tag \n" "$1"
}

if [[ "$#" == 0 ]]; then
  usage "$0"; exit 0;
fi

#################################################
# Parsing arguments
#################################################
while getopts "d:r:n:t:" option; do
  case "$option" in
    d) BASE_DIR="$OPTARG" ;;
    r) DOCKER_ROOT="$OPTARG" ;;
    n) DOCKER_NAME="$OPTARG" ;;
    t) DOCKER_TAG="$OPTARG" ;;
    *) usage "$0" && exit 1 ;;
  esac
done

if [[ -z "$BASE_DIR" ]] || [[ -z "$DOCKER_ROOT" ]] || [[ -z "$DOCKER_NAME" ]] || [[ -z "$DOCKER_TAG" ]]; then
    usage "$0"
    exit 1
fi

if [[ ! -d "$BASE_DIR" ]]; then
   echo "Invalid directory: $BASE_DIR"
   exit 1
fi

if [[ ${DOCKER_NAME} == "ALL" ]]; then
  DOCKER_NAME_ARR=("cnmops" "delly" "manta" "samtools-cloud" "sv-base" "sv-base-mini" "sv-pipeline" "sv-pipeline-base" \
    "sv-pipeline-qc" "sv-pipeline-rdtest" "wham")
  echo "Warning: MELT dockers will not be updated."
else
  DOCKER_NAME_ARR=("$DOCKER_NAME")
fi

#################################################
# Update jsons
#################################################

shopt -s nullglob
JSON_ARR=(${BASE_DIR}/*.json ${BASE_DIR}/test/*/*.json ${BASE_DIR}/inputs/*.json)
JSONS=$(printf "%s "  "${JSON_ARR[@]}")

for name in "${DOCKER_NAME_ARR[@]}"; do
  docker_var=$(tr '-' '_' <<< "${name}")"_docker"
  docker_regex='\.'"${docker_var}"'"\s*:\s*".+\/.+:.+"'
  cmd="perl -pi -e 's|${docker_regex}|.${docker_var}\": \"${DOCKER_ROOT}\/${name}:${DOCKER_TAG}\"|g' ${JSONS}"
  echo "${cmd}"
  eval "${cmd}"
done
