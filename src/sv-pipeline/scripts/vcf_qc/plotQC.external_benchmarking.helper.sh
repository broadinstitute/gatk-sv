#!/bin/bash

# Wrapper to manage external SV benchmarking plotting from output tarball

set -e

###USAGE
usage(){
cat <<EOF

usage: plotQC_external_benchmarking.helper.sh [-h] TARBALL COMPARATOR

Wrapper to manage external SV benchmarking plotting from output tarball

Positional arguments:
  TARBALL         .tar.gz file containing results of external benchmarking from
                  collectQC_external_benchmarking WDL
  COMPARATOR      Comparison dataset used in benchmarking. Specify one of the
                  following: 'ASC_Werling', 'HGSV_Chaisson', or '1000G_Sudmant'

Optional arguments:
  -h  HELP        Show this help message and exit

EOF
}


###PARSE ARGS
TARBALL=$1
COMPARATOR=$2


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
if [ ${COMPARATOR} != "ASC_Werling" ] && [ ${COMPARATOR} != "HGSV_Chaisson" ] && \
   [ ${COMPARATOR} != "1000G_Sudmant" ]; then
  echo -e "\nERROR: COMPARATOR must be one of 'ASC_Werling', 'HGSV_Chaisson', or '1000G_Sudmant'\n"
  usage
  exit 0
fi


###SET BIN
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


###PREP INPUT FILES
#Unarchive tarball
tar -xzvf ${TARBALL} --directory `pwd`
#Prep directory prefix
OUTDIR=`pwd`/collectQC_benchmarking_${COMPARATOR}_output/
#Prep plot directory
if ! [ -e ${OUTDIR}/plots/ ]; then
  mkdir ${OUTDIR}/plots/
fi


###GENERATE PLOTS FOR EACH POPULATION
#Set carrier frequency flag for ASC Werling and HGSV Chaisson
if [ ${COMPARATOR} == "ASC_Werling" ] || [ ${COMPARATOR} == "HGSV_Chaisson" ]; then
  carrierFlag=" -C"
else
  carrierFlag=""
fi
#Iterate over all populations and generate plots
for pop in ALL AFR AMR EAS EUR SAS OTH; do
  compdat=${OUTDIR}/data/${COMPARATOR}.SV.${pop}.overlaps.bed.gz
  if [ -e ${compdat} ] && [ -s ${compdat} ]; then
    #Print status
    echo -e "$( date ) - VCF QC STATUS: Plotting benchmarking for ${pop} samples in ${COMPARATOR}"
    #Plot benchmarking
    ${BIN}/plot_callset_comparison.R \
      ${carrierFlag} \
      -p ${COMPARATOR}_${pop} \
      ${compdat} \
      ${OUTDIR}/plots/${COMPARATOR}_${pop}_samples/
  fi
done

