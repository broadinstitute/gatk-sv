#!/bin/bash

#################################################
#################################################
# Validates all WDL JSON input files
#################################################
#################################################

set -e

#################################################
# Utility functions
#################################################

function usage() {
  printf "Usage: \n \
    %s -d <REPO_BASE_DIR> -j <WOMTOOL_JAR> \n \
    <REPO_BASE_DIR> \t path to gatk-sv base directory \n \
    <WOMTOOL_JAR> \t path to womtool jar (downloaded from https://github.com/broadinstitute/cromwell/releases) \n" "$1"
}

if [[ "$#" == 0 ]]; then
  usage "$0"; exit 0;
fi

#################################################
# Parsing arguments
#################################################
while getopts "j:d:" option; do
  case "$option" in
    j) WOMTOOL_JAR="$OPTARG" ;;
    d) BASE_DIR="$OPTARG" ;;
    *) usage "$0" && exit 1 ;;
  esac
done

if [ -z "$WOMTOOL_JAR" ] || [ -z "$BASE_DIR" ] ; then
    usage "$0"
    exit 1
fi

if [[ ! -d "$BASE_DIR" ]]; then
   echo "Invalid directory: $BASE_DIR"
   exit 1
fi

if [[ ! -f "$WOMTOOL_JAR" ]]; then
   echo "Invalid file: $BASE_DIR"
   exit 1
fi

#################################################
# For each WDL, test all jsons with matching name
#################################################
shopt -s nullglob
shopt -s extglob

WDLS=(wdl/*.wdl)

COUNTER=0
for wdl in "${WDLS[@]}"
do
  name=$(basename $wdl .wdl)
  JSONS=(${BASE_DIR}/test_inputs/small/*/${name}*(.*).json ${BASE_DIR}/test_inputs/large/*/${name}*(.*).json ${BASE_DIR}/test_inputs/single_sample/*/${name}*(.*).json ${BASE_DIR}/inputs/${name}*(.*).json)
  for json in "${JSONS[@]}"
  do
    cmd="java -jar ${WOMTOOL_JAR} validate ${wdl} -i ${json}"
    echo $cmd
    eval $cmd
    echo "PASS"
    echo ""
    COUNTER=$((COUNTER+1))
  done
done

if [ ${COUNTER} -eq 0 ]; then
  echo "Did not find any WDLs with matching jsons, check -d argument."
  exit 1
fi

echo ""
echo "#############################################################"
echo "${COUNTER} TESTS PASSED SUCCESSFULLY!"
