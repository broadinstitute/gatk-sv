#!/bin/bash

# Wrapper to manage external SV benchmarking plotting from output tarball

set -e

###USAGE
usage(){
cat <<EOF

usage: plotQC_external_benchmarking.helper.sh [-h] TARBALL PREFIX

Wrapper to manage external SV benchmarking plotting from output tarball

Positional arguments:
  TARBALL         .tar.gz file containing results of external benchmarking from
                  collectQC_external_benchmarking WDL
  PREFIX          String prefix for labeling outputs

Optional arguments:
  -h  HELP        Show this help message and exit

EOF
}


###PARSE ARGS
TARBALL=$1
PREFIX=$2
TARBALL_DIR=$3


###PROCESS ARGS & OPTIONS
#Check for required input
if [ -z ${TARBALL} ]; then
  echo -e "\nERROR: input .tar.gz file not specified\n"
  usage
  exit 0
fi
if ! [ -s ${TARBALL} ]; then
  echo -e "\nERROR: input .tar.gz either empty or not found\n"
  usage
  exit 0
fi
#Check for required input
if [ -z ${PREFIX} ]; then
  echo -e "\nERROR: output prefix not specified\n"
  usage
  exit 0
fi
#Check for required input
if [ -z ${TARBALL_DIR} ]; then
  echo -e "\nERROR: tarball directory name not specified\n"
  usage
  exit 0
fi


###SET BIN
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


###PREP INPUT FILES
#Unarchive tarball
tar -xzvf ${TARBALL} --directory `pwd`
#Prep directory prefix
OUTDIR=`pwd`/${TARBALL_DIR}
#Prep plot directory
if ! [ -e ${OUTDIR}/plots/ ]; then
  mkdir ${OUTDIR}/plots/
fi


###GENERATE PLOTS FOR EACH POPULATION
find ${OUTDIR}/ -name "*.overlaps.bed.gz" | sed '/^$/d' > input_data.tsv
cat input_data.tsv
while read compdat; do
  outprefix=$( basename ${compdat} | sed 's/\.overlaps\.bed\.gz//g' )
  if [ -e ${compdat} ] && [ -s ${compdat} ]; then
    #Print status
    echo -e "$( date ) - VCF QC STATUS: Plotting benchmarking for ${outprefix} subset from ${PREFIX}"
    #Plot benchmarking
    ${BIN}/plot_callset_comparison.R \
      -p ${outprefix} \
      ${compdat} \
      ${OUTDIR}/plots/${outprefix}/
  fi
done < input_data.tsv

