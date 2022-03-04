#/bin/bash

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

rm -r ${BASE_DIR}/inputs/build
