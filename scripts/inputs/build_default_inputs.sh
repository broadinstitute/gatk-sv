#!/bin/bash

function usage() {
  printf "Usage: \n \
    %s -d <REPO_BASE_DIR> \n \
    <REPO_BASE_DIR> \t path to gatk-sv base directory \n" "$1"
}

if [[ "$#" == 0 ]]; then
  usage "$0"; exit 0;
fi

#################################################
# Parsing arguments
#################################################
while getopts "d:" option; do
  case "$option" in
    d) BASE_DIR="$OPTARG" ;;
    *) usage "$0" && exit 1 ;;
  esac
done

if [ -z "$BASE_DIR" ] ; then
    usage "$0"
    exit 1
fi

if [[ ! -d "$BASE_DIR" ]]; then
   echo "Invalid directory: $BASE_DIR"
   exit 1
fi

scripts/inputs/build_inputs.py ${BASE_DIR}/input_values ${BASE_DIR}/input_templates ${BASE_DIR}/inputs -a '{"ref_panel" : "ref_panel_1kg_v2"}'

scripts/inputs/build_inputs.py ${BASE_DIR}/input_values ${BASE_DIR}/test_input_templates ${BASE_DIR}/test_inputs_small -a '{"test_batch" : "test_batch_small"}'
scripts/inputs/build_inputs.py ${BASE_DIR}/input_values ${BASE_DIR}/test_input_templates ${BASE_DIR}/test_inputs_large -a '{"test_batch" : "test_batch_large"}'

scripts/inputs/build_inputs.py ${BASE_DIR}/input_values ${BASE_DIR}/test_input_templates ${BASE_DIR}/test_inputs_single_sample -a '{ "test_batch" : "test_single_sample_NA19240", "ref_panel" : "ref_panel_v1b" }'

scripts/inputs/build_inputs.py ${BASE_DIR}/input_values ${BASE_DIR}/terra_workspace_templates ${BASE_DIR}/terra_workspace -a '{"ref_panel" : "ref_panel_1kg_v2"}'