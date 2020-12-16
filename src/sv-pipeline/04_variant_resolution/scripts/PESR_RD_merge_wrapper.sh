#!/bin/bash

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Wrapper for PESR + depth merging in module 04b

set -ex

###USAGE
usage(){
cat <<EOF

usage: PESR_RD_merge_wrapper.sh [-h] PESR_VCF RD_VCF CHROM OUT_VCF

Wrapper to handle PESR + RD merging in module 04b

Positional arguments:
  PESR_VCF                 Input PESR VCF, bgzipped
  RD_VCF                   Input RD VCF, bgzipped
  CHROM                    Chromosome being processed
  OUT_VCF                  Output merged VCF, bgzipped

Optional arguments:
  -h  HELP                 Show this help message and exit

Notes:
  1. All input files must be compressed with bgzip.

EOF
}


###PARSE ARGS
while getopts ":h" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
  esac
done
shift $(( ${OPTIND} - 1))
PESR=$1
RD=$2
chrom=$3
OUT=$4


###PROCESS ARGS & OPTIONS
#Check for required input
if [ -z ${PESR} ]; then
  echo -e "\nERROR: input PESR VCF not specified\n"
  usage
  exit 0
fi
if ! [ -s ${PESR} ]; then
  echo -e "\nERROR: input PESR VCF either empty or not found\n"
  usage
  exit 0
fi
if [ $( file ${PESR} | fgrep "gzip" | wc -l ) -lt 1 ]; then
  echo -e "\nERROR: input PESR VCF must be bgzipped\n"
  usage
  exit 0
fi
if [ -z ${RD} ]; then
  echo -e "\nERROR: input RD VCF not specified\n"
  usage
  exit 0
fi
if ! [ -s ${RD} ]; then
  echo -e "\nERROR: input RD VCF either empty or not found\n"
  usage
  exit 0
fi
if [ $( file ${RD} | fgrep "gzip" | wc -l ) -lt 1 ]; then
  echo -e "\nERROR: input RD VCF must be bgzipped\n"
  usage
  exit 0
fi
if [ -z ${OUT} ]; then
  echo -e "\nERROR: path to output VCF not specified\n"
  usage
  exit 0
fi
#Set path to execution directory
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#Prepares temporary directory
PROCDIR=`mktemp -d`


#Get list of calls with 100% ref genotypes
#This is necessary to protect against X/Y false negatives before cleanVCF
#svtk vcf2bed ${PESR} ${PROCDIR}/PESR.bed
#awk -v FS="\t" '{ if ($6=="") print $4 }' ${PROCDIR}/PESR.bed \
#  | sort | uniq \
#  > ${PROCDIR}/VIDs_with_no_samples.pesr.list
#
#svtk vcf2bed ${RD} ${PROCDIR}/RD.bed
#awk -v FS="\t" '{ if ($6=="") print $4 }' ${PROCDIR}/RD.bed \
#  | sort | uniq \
#  > ${PROCDIR}/VIDs_with_no_samples.depth.list

# Merge VCFs (doesn't merge variants)
bcftools concat -a --output-type z \
  ${PESR} ${RD} --output ${PROCDIR}/pesr_depth_unmerged.${chrom}.vcf.gz

#Run merging script (merges variants)
${BIN}/merge_pesr_depth.py \
  --prefix pesr_depth_merged_${chrom} \
  ${PROCDIR}/pesr_depth_unmerged.${chrom}.vcf.gz \
  ${PROCDIR}/pesr_depth_merged_unsorted.${chrom}.vcf

rm ${PROCDIR}/pesr_depth_unmerged.${chrom}.vcf.gz

# Sort output
bcftools sort ${PROCDIR}/pesr_depth_merged_unsorted.${chrom}.vcf -Oz -o ${OUT}

##Get list of members in merged VCF
#if [ $( cat ${PROCDIR}/VIDs_with_no_samples.pesr.list | wc -l ) -gt 0 ] || \
#   [ $( cat ${PROCDIR}/VIDs_with_no_samples.depth.list | wc -l ) -gt 0 ]; then
#  svtk vcf2bed -i MEMBERS --no-samples \
#    ${PROCDIR}/pesr_depth_merged.${chrom}.vcf \
#    ${PROCDIR}/pesr_depth_merged.${chrom}.bed
#  awk '{ print $NF }' ${PROCDIR}/pesr_depth_merged.${chrom}.bed \
#    | sed 's/,/\n/g' | sort -Vk1,1 | uniq \
#    > ${PROCDIR}/VIDs_in_merged_vcf.list
#fi
#
#
##Salvage 100% ref PESR calls
#if [ $( cat ${PROCDIR}/VIDs_with_no_samples.pesr.list | wc -l ) -gt 0 ]; then
#  cat <(zgrep '^#' ${PESR}) \
#    <( fgrep -wf <(fgrep -wvf ${PROCDIR}/VIDs_in_merged_vcf.list ${PROCDIR}/VIDs_with_no_samples.pesr.list) <( zgrep -v "^#" ${PESR} ) ) \
#    | ${BIN}/salvage_allRef_calls.py -p salvaged_pesr_allref /dev/stdin /dev/stdout \
#    | grep -v "^#" \
#    >> ${PROCDIR}/pesr_depth_merged.${chrom}.vcf
#fi
#
#
##Salvage 100% ref depth calls
#if [ $( cat ${PROCDIR}/VIDs_with_no_samples.depth.list | wc -l ) -gt 0 ]; then
#  cat <( zgrep "^#" ${RD} ) \
#      <( fgrep -wf <(fgrep -wvf ${PROCDIR}/VIDs_in_merged_vcf.list ${PROCDIR}/VIDs_with_no_samples.depth.list) <( zgrep -v "^#" ${RD} ) ) \
#    | ${BIN}/salvage_allRef_calls.py -p salvaged_depth_allref /dev/stdin /dev/stdout \
#    | grep -v "^#" \
#    >> ${PROCDIR}/pesr_depth_merged.${chrom}.vcf
#fi
#
#
##Sort and compress final output
#vcf-sort ${PROCDIR}/pesr_depth_merged.${chrom}.vcf \
#  | bgzip -c \
#  > ${OUT}


#Clean up
rm -rf ${PROCDIR}

