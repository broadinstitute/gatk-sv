#!/bin/bash

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Resolve redundancies between simple CNVs and unbalanced complex SV in mod04b

# TODO : Missing pipefail
set -e

###USAGE
usage(){
cat <<EOF

usage: resolve_CPX_CNV_redundancies.sh [-h] INVCF OUTVCF

Resolve redundancies between simple CNVs and unbalanced complex SV in mod04b

Positional arguments:
  INVCF                    Original input VCF prior to regenotyping
  OUTVCF                   Full path to output VCF after relabeling

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
INVCF=$1
OUTVCF=$2


###PROCESS ARGS & OPTIONS
#Check for required input
if [ -z ${INVCF} ]; then
  echo -e "\nERROR: input VCF not specified\n"
  usage
  exit 0
fi
if ! [ -s ${INVCF} ]; then
  echo -e "\nERROR: input VCF either empty or not found\n"
  usage
  exit 0
fi
if [ $( file ${INVCF} | fgrep "gzip" | wc -l ) -lt 1 ]; then
  echo -e "\nERROR: input VCF must be bgzipped\n"
  usage
  exit 0
fi
if [ -z ${OUTVCF} ]; then
  echo -e "\nERROR: path to output VCF not specified\n"
  usage
  exit 0
fi
#Set path to execution directory
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#Prepares temporary directory
PROCDIR=`mktemp -d`


###PREP FILES
#Convert full VCF to BED intervals
#Ignore CPX events with UNRESOLVED filter status
svtk vcf2bed --split-cpx --info SVTYPE \
  <(bcftools view -e 'INFO/SVTYPE == "CPX" && FILTER == "UNRESOLVED"' ${INVCF}) - \
  | grep -e '^#\|DEL\|DUP\|CNV\|CPX' \
  | awk -v OFS="\t" '{ if ($5=="CN0") print $1, $2, $3, $4, "DEL", $5"\n"$1, $2, $3, $4, "DUP", $5; \
                       else if ($5=="DEL" || $5=="DUP") print $1, $2, $3, $4, $6, $5 }' \
  | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
  | bgzip -c \
  > ${PROCDIR}/intervals.preclustered.bed.gz


###REMOVE CNVS REDUNDANT WITH COMPLEX EVENTS
#Subset to only variants that share some overlap (at least 10% recip) with at least one CPX variant
bedtools intersect -wa -r -f 0.1 \
  -a ${PROCDIR}/intervals.preclustered.bed.gz \
  -b <( zcat ${PROCDIR}/intervals.preclustered.bed.gz | fgrep "CPX" ) \
  | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
  | uniq \
  | bgzip -c \
  > ${PROCDIR}/intervals.preclustered.subset.bed.gz
#Melt subsetted variants
while read chr start end VID samples CNV; do
  echo -e "${samples}" \
    | sed 's/,/\n/g' \
    | awk -v OFS="\t" -v chr=${chr} -v start=${start} -v end=${end} -v VID=${VID} -v CNV=${CNV} \
      '{ print chr, start, end, VID, $1, CNV }'
done < <( zcat ${PROCDIR}/intervals.preclustered.subset.bed.gz ) \
  | bgzip -c \
  > ${PROCDIR}/intervals.preclustered.subset.melted.bed.gz
#Cluster BED intervals (50% RO)
svtk bedcluster -f 0.5 \
  ${PROCDIR}/intervals.preclustered.subset.melted.bed.gz - \
  | bgzip -c > \
  ${PROCDIR}/intervals.clustered.bed.gz
#Get list of all variants that cluster with a complex variant, 
# evaluate sample overlap from original intervals file, 
# and, if overlap >50%, write that ID to be stripped from the output VCF
while read VIDs; do
  #Get nonredundant list of sample IDs involved in any clustered variant
  echo -e "${VIDs}" | sed 's/,/\n/g' \
    | fgrep -wf - <( zcat ${PROCDIR}/intervals.preclustered.bed.gz ) \
    | cut -f5 | sort | uniq \
    > ${PROCDIR}/nonredundant_samples.list
  #Iterate over VIDs and print non-CPX VID if sample overlap >50%
  while read VID samples; do
    #Get list of samples in variant
    echo -e "${samples}" | sed 's/,/\n/g' \
      | sort | uniq > ${PROCDIR}/query_samples.list
    nsamp=$( cat ${PROCDIR}/query_samples.list | wc -l )

    #Compare
    frac=$( fgrep -wf ${PROCDIR}/query_samples.list \
              ${PROCDIR}/nonredundant_samples.list | wc -l \
              | awk -v nsamp=${nsamp} '{ print 100*($1/nsamp) }' \
              | cut -f1 -d\. )
    if [ ${frac} -ge 50 ]; then
      echo "${VID}"
    fi

    #Clean up
    rm ${PROCDIR}/query_samples.list

  done < <( echo -e "${VIDs}" | sed 's/,/\n/g' \
              | fgrep -wf - <( zcat ${PROCDIR}/intervals.preclustered.bed.gz ) \
              | cut -f4,5 | sort | uniq | fgrep -v "CPX" )

  #Clean up
  rm ${PROCDIR}/nonredundant_samples.list

done < <( zcat ${PROCDIR}/intervals.clustered.bed.gz \
            | cut -f7 | fgrep "CPX" | grep -e "DEL\|DUP" ) \
  | sort -V | uniq \
  > ${PROCDIR}/VIDs_to_remove.list


###FIND REMAINING REDUNDANT CNVS WITH STRONG (80%) OVERLAP IN SAMPLES AND SIZE
#Find CNV intervals that have 80% reciprocal overlap
bedtools intersect -wa -wb -r -f 0.8 \
  -a ${PROCDIR}/intervals.preclustered.bed.gz \
  -b ${PROCDIR}/intervals.preclustered.bed.gz \
  | awk -v FS="\t" '{ if ($4!=$10 && $6==$12) print $0 }' \
  | awk -v OFS="\t" '$4 ~ /DEL|DUP/ { print $0 }' \
  | awk -v OFS="\t" '$10 ~ /DEL|DUP/ { print $0 }' \
  | bgzip -c \
  > ${PROCDIR}/step2.intervals.preclustered.subset.bed.gz
#Determine which events share 80% sample overlap
while read VIDa sa VIDb sb; do
  na=$( echo -e "${sa}" | sed 's/,/\n/g' | sort | uniq | wc -l )
  nb=$( echo -e "${sb}" | sed 's/,/\n/g' | sort | uniq | wc -l )
  denom=$( echo -e "${sa},${sb}" | sed 's/,/\n/g' | sort | uniq | wc -l )
  numer=$( echo -e "${sa}" | sed 's/,/\n/g' | fgrep -wf - \
            <( echo -e "${sb}" | sed 's/,/\n/g' ) \
            | sort | uniq | wc -l )
  if [ ${denom} -gt 0 ]; then
    ovr=$(( 100 * ${numer} / ${denom} ))
  fi
  if [ -z ${ovr} ]; then
    ovr=0
  fi
  if [ ${ovr} -ge 80 ]; then
    echo -e "${VIDa}\n${VIDb}" \
      | sort | uniq | paste -s -d,
  fi
done < <( zcat ${PROCDIR}/step2.intervals.preclustered.subset.bed.gz \
            | cut -f4,5,10,11 ) \
  | sort | uniq \
  > ${PROCDIR}/step2.variants_to_resolve.list
#Iterate over variants, pick info & coords from variant with largest N, 
# and consolidate genotypes
sed 's/,/\n/g' ${PROCDIR}/step2.variants_to_resolve.list \
  | sort | uniq \
  > ${PROCDIR}/step2.variants_to_resolve.melted.list
if [ -e ${PROCDIR}/records_to_add.vcf ]; then
  rm ${PROCDIR}/records_to_add.vcf
fi
until [ $( cat ${PROCDIR}/step2.variants_to_resolve.melted.list | wc -l ) -eq 0 ]; do
  #get next variant
  VID=$( head -n1 ${PROCDIR}/step2.variants_to_resolve.melted.list )
  #get all other variants from clusters containing this variant
  fgrep -w ${VID} ${PROCDIR}/step2.variants_to_resolve.list \
    | sed 's/,/\n/g' | sort | uniq \
    > ${PROCDIR}/step2.partners.tmp
  #Print all genotypes to tmp file
  zcat ${INVCF} | fgrep -v "#" \
    | fgrep -wf ${PROCDIR}/step2.partners.tmp | cut -f10- \
    > ${PROCDIR}/gts.tmp
  #Select best genotypes to keep
  ${BIN}/selectBestGT.R ${PROCDIR}/gts.tmp ${PROCDIR}/gts.best.tmp
  #Select record with greatest total number of samples
  bVID=$( zcat ${PROCDIR}/intervals.preclustered.bed.gz \
            | fgrep -wf ${PROCDIR}/step2.partners.tmp \
            | cut -f4-5 | sed 's/,/\t/g' \
            | awk -v OFS="\t" '{ print $1, NF }' \
            | sort -nrk2,2 \
            | cut -f1 \
            | head -n1 )
  #Add new record to final append tmp file
  paste <( zcat ${INVCF} | fgrep -w ${bVID} | cut -f1-9 ) \
        ${PROCDIR}/gts.best.tmp \
    >> ${PROCDIR}/records_to_add.vcf
  #Write list of variants to exclude from original VCF
  cat ${PROCDIR}/step2.partners.tmp >> ${PROCDIR}/VIDs_to_remove.list
  #Exclude variants from list of VIDs to resolve
  fgrep -wvf ${PROCDIR}/step2.partners.tmp \
    ${PROCDIR}/step2.variants_to_resolve.melted.list \
    > ${PROCDIR}/step2.variants_to_resolve.melted.list2 \
    || true
  mv ${PROCDIR}/step2.variants_to_resolve.melted.list2 \
    ${PROCDIR}/step2.variants_to_resolve.melted.list
done


###CLEAN UP FINAL OUTPUT
zcat ${INVCF} \
  | fgrep -wvf ${PROCDIR}/VIDs_to_remove.list \
  | cat - ${PROCDIR}/records_to_add.vcf \
  | vcf-sort \
  | bgzip -c \
  > ${OUTVCF}


###CLEAN UP
rm -rf ${PROCDIR}


