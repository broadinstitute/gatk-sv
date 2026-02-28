#!/bin/bash

# Stitch fragmented CNVs present in a VCF

###USAGE
usage(){
cat <<EOF

usage: stitch_fragmented_CNVs.sh [-h] [-d] [-e <float>] [-m <int>] [-X <float>] INVCF OUTVCF

Stitch fragmented CNVs present in the same samples in a VCF

Positional arguments:
  INVCF                    Original input VCF
  OUTVCF                   Full path to output VCF

Optional arguments:
  -h  HELP                 Show this help message and exit
  -d  DONT CHECK EVIDENCE  Don't check the EVIDENCE annotations (default: only process calls which have 'RD' or 'BAF'
                           evidence and not 'PE' or 'SR')
  -e  EXTEND               Fraction of total variant size to pad at 
                           ends of variants (default: 20%)
  -m  MAX PADDING          Maximum overlap (in bp) to be padded at the
                           ends of variants (default: 200kb)
  -x  MAX OVERLAP          Allow no more than X% overlap between two 
                           variants (default: 20%)


Notes:
  1. All input files must be compressed with bgzip.

EOF
}


#Strict error reporting
set -eu


###PARSE ARGS
EXTEND=0.2
MAXPAD=200000
MAXOVR=0.2
CHECK_EVIDENCE=true
while getopts ":e:m:x:h" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    d)
      CHECK_EVIDENCE=false
      ;;
    e)
      EXTEND=${OPTARG}
      ;;
    m)
      MAXPAD=${OPTARG}
      ;;
    x)
      MAXOVR=${OPTARG}
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
if [ -z ${OUTVCF} ]; then
  echo -e "\nERROR: path to output VCF not specified\n"
  usage
  exit 0
fi
#Set path to execution directory
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#Prepares temporary directory
PROCDIR=`mktemp -d`


#Subset VCF to biallelic CNVs without PE or SR support
zcat ${INVCF} \
  | fgrep -v "#" \
  | grep -e 'SVTYPE=DEL\|SVTYPE=DUP' \
  | fgrep -v "MULTIALLELIC" \
  | awk -v OFS="\t" '{ print $3, $8 }' \
  | sed 's/\;/\t/g' \
  | cut -f1 \
  | fgrep -wf - <( zcat ${INVCF} ) \
  | cat <( zcat ${INVCF} | fgrep "#" ) - \
  | vcf-sort \
  | bgzip -c \
  > ${PROCDIR}/biallelic_depth_CNVs.vcf.gz


#Convert to BED
svtk vcf2bed \
  ${PROCDIR}/biallelic_depth_CNVs.vcf.gz \
  ${PROCDIR}/biallelic_depth_CNVs.bed


#Find CNVs appearing in the same samples that meet merging criteria
for CNV in DEL DUP; do
  awk -v OFS="\t" -v CNV=${CNV} -v MAXOVR=${MAXOVR} \
    '{ if ($5==CNV) print $1"_"$6, $2, $3, $4, $2, $3, $3-$2, MAXOVR*($3-$2) }' \
    ${PROCDIR}/biallelic_depth_CNVs.bed \
    | sed 's/\.[0-9]+\t/\t/g' \
    | awk -v OFS="\t" -v MAXPAD=${MAXPAD} '{ if ($NF>MAXPAD) $NF=MAXPAD; print $1, $2, $3, $4, $5, $6, $7, $8 }' \
    | awk -v OFS="\t" '{ print $1, $2-$NF, $3+$NF, $4, $5, $6, $7 }' \
    | awk -v OFS="\t" '{ if ($2<1) $2=1; printf "%s\t%.0f\t%.0f\t%s\t%.0f\t%.0f\t%.0f\n", $1, $2, $3, $4, $5, $6, $7 }' \
    | sort -Vk1,1 -k2,2n -k3,3n \
    > ${PROCDIR}/ovr.tmp.bed
  if [ $( cat ${PROCDIR}/ovr.tmp.bed | wc -l ) -gt 0 ]; then
    bedtools intersect -wa -wb \
      -a ${PROCDIR}/ovr.tmp.bed \
      -b ${PROCDIR}/ovr.tmp.bed \
      | awk -v OFS="\t" -v MAXOVR=${MAXOVR} \
        '{ if ($4!=$11 && $6-$12<=(MAXOVR*$7) && $6-$12<=(MAXOVR*$14) ) print $1, $2, $3, $4"\n"$8, $9, $10, $11 }' \
      | sort -Vk1,1 -k2,2n -k3,3n \
      > ${PROCDIR}/int_to_merge.bed
    if [ $( cat ${PROCDIR}/int_to_merge.bed | wc -l ) -gt 0 ]; then
      bedtools merge -c 4 -o distinct \
        -i ${PROCDIR}/int_to_merge.bed \
        | awk '{ if ($4~",") print $4 }'
    fi
  fi
done > ${PROCDIR}/variants_to_merge.list


#Merge CNVs that meet criteria and write modified variants to new set of records
if [ -e ${PROCDIR}/records_to_add.vcf ]; then
  rm ${PROCDIR}/records_to_add.vcf
fi
if [ -e ${PROCDIR}/records_to_remove.VIDs.list ]; then
  rm ${PROCDIR}/records_to_remove.VIDs.list
fi
if [ $( cat ${PROCDIR}/variants_to_merge.list | wc -l ) -gt 0 ]; then
  while read cluster; do
    #Get new start and end coordinates and new length
    nstart=$( echo "${cluster}" \
                | sed 's/,/\n/g' \
                | fgrep -wf - ${PROCDIR}/biallelic_depth_CNVs.bed \
                | cut -f2 \
                | sort -nk1,1 \
                | head -n1 )
    nend=$( echo "${cluster}" \
              | sed 's/,/\n/g' \
              | fgrep -wf - ${PROCDIR}/biallelic_depth_CNVs.bed \
              | cut -f3 \
              | sort -nrk1,1 \
              | head -n1 )
    nlen=$(( ${nend} - ${nstart} ))
    #Take first record from VCF, modify start & end coordinates, and save to be added later
    echo "${cluster}" \
      | sed 's/,/\n/g' \
      | fgrep -wf - <( zcat ${PROCDIR}/biallelic_depth_CNVs.vcf.gz ) \
      > ${PROCDIR}/tmp_vars_to_cluster.vcf
    head -n 1 ${PROCDIR}/tmp_vars_to_cluster.vcf \
      | awk -v FS="\t" -v OFS="\t" -v nstart=${nstart} \
        '{ $2=nstart; print }' \
      | sed -r -e "s/END=[0-9]*;/END=${nend};/g" \
        -e "s/SVLEN=[0-9]*;/SVLEN=${nlen};/g" \
      >> ${PROCDIR}/records_to_add.vcf
    #Add all record IDs to be stripped from the original VCF
    echo "${cluster}" \
      | sed 's/,/\n/g' \
      >> ${PROCDIR}/records_to_remove.VIDs.list
  done < ${PROCDIR}/variants_to_merge.list
fi


###FORMAT FINAL VCF
if [ $( cat ${PROCDIR}/variants_to_merge.list | wc -l ) -gt 0 ]; then
  zcat ${INVCF} \
    | fgrep -wvf ${PROCDIR}/records_to_remove.VIDs.list \
    | cat - ${PROCDIR}/records_to_add.vcf \
    | vcf-sort \
    | bgzip -c \
    > ${OUTVCF}
else
  cp ${INVCF} ${OUTVCF}
fi


###CLEAN UP
rm -rf ${PROCDIR}

