#!/bin/bash

function usage() {
  printf "Usage: \n \
    %s -d <REPO_BASE_DIR> -c <CLOUD_ENV> \n \
    <REPO_BASE_DIR> \t path to gatk-sv base directory \n \
    <CLOUD_ENV> \t name of cloud environment json (e.g. 'google_cloud.my' for  inputs/values/google_cloud.my.json)" "$1"
}

if [[ "$#" == 0 ]]; then
  usage "$0"; exit 0;
fi

#################################################
# Parsing arguments
#################################################
while getopts "d:c:" option; do
  case "$option" in
    d) BASE_DIR="$OPTARG" ;;
    c) CLOUD_ENV="$OPTARG" ;;
    *) usage "$0" && exit 1 ;;
  esac
done

if [ -z "$BASE_DIR" ] ; then
    echo "x"
    usage "$0"
    exit 1
fi

if [ -z "$CLOUD_ENV" ] ; then
    echo "xy"
    usage "$0"
    exit 1
fi

if [[ ! -d "$BASE_DIR" ]]; then
   echo "Invalid directory: $BASE_DIR"
   exit 1
fi


# Clean any existing outputs
bash scripts/inputs/clean_default_inputs.sh -d ${BASE_DIR}

echo "########## Building ref_panel_1kg test ##########"
scripts/inputs/build_inputs.py ${BASE_DIR}/inputs/values ${BASE_DIR}/inputs/templates/test ${BASE_DIR}/inputs/build/ref_panel_1kg/test \
  -a '{ "test_batch" : "ref_panel_1kg", "cloud_env" : "'$CLOUD_ENV'" }'

echo "########## Building ref_panel_1kg cohort Terra workspace ##########"
scripts/inputs/build_inputs.py ${BASE_DIR}/inputs/values ${BASE_DIR}/inputs/templates/terra_workspaces/cohort_mode ${BASE_DIR}/inputs/build/ref_panel_1kg/terra \
  -a '{ "test_batch" : "ref_panel_1kg", "cloud_env" : "'$CLOUD_ENV'" }'

# Note CLOUD_ENV is not currently required for the single-sample workflow
echo "########## Building NA19240 single-sample test ##########"
scripts/inputs/build_inputs.py ${BASE_DIR}/inputs/values ${BASE_DIR}/inputs/templates/test/GATKSVPipelineSingleSample ${BASE_DIR}/inputs/build/NA19240/test \
  -a '{ "single_sample" : "test_single_sample_NA19240", "ref_panel" : "ref_panel_1kg" }'

echo "########## Building NA12878 single-sample test ##########"
scripts/inputs/build_inputs.py ${BASE_DIR}/inputs/values ${BASE_DIR}/inputs/templates/test/GATKSVPipelineSingleSample ${BASE_DIR}/inputs/build/NA12878/test \
  -a '{ "single_sample" : "test_single_sample_NA12878", "ref_panel" : "ref_panel_1kg" }'

echo "########## Building NA12878 single-sample Terra workspace ##########"
scripts/inputs/build_inputs.py ${BASE_DIR}/inputs/values ${BASE_DIR}/inputs/templates/terra_workspaces/single_sample ${BASE_DIR}/inputs/build/NA12878/terra \
  -a '{ "single_sample" : "test_single_sample_NA12878", "ref_panel" : "ref_panel_1kg" }'
