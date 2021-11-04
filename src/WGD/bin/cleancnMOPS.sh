#!/bin/bash

# Simple script to split cn.MOPS calls by sample and copy state and merge
# across multiple runs (e.g. at different resolutions)

set -e

#Usage statement
usage(){
cat <<EOF

usage: cleancnMOPS.sh [-h] [-z] [-o OUTDIR] [-S SUBTRACT] SAMPLES GFFS

Helper tool to split cn.MOPS calls by sample and copy state and merge
across multiple runs (e.g. at different resolutions)

Positional arguments:
  SAMPLES   list of sample IDs
  GFFS      list of full paths to all cn.MOPS outputs to be considered

Optional arguments:
  -h  HELP        Show this help message and exit
  -z  GZIP        Gzip output file
  -o  OUTDIR      Output directory
  -S  SUBTRACT    Intervals to subtract from any overlapping calls 
                  (e.g. N-masked reference gaps) 

EOF
}

#Parse arguments
OUTDIR=`pwd`
GZ=0
SUBTRACT=0
while getopts ":o:zS:h" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    o)
      OUTDIR=${OPTARG}
      ;;
    z)
      GZ=1
      ;;
    S)
      SUBTRACT=${OPTARG}
      ;;
  esac
done
shift $(( ${OPTIND} - 1))
SAMPLES=$1
GFFS=$2

#Check for required input
if [ -z ${SAMPLES} ] || [ -z ${GFFS} ]; then
  usage
  exit 0
fi

#Attempts to create OUTDIR
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi

#Splits raw GFFs by copy state and formats to BED
#Deletions
DEL_MASTER=`mktemp`
while read gff; do
  sed 's/\;/\t/g' ${gff} | fgrep -v "#" | sed 's/CN=CN//g' | \
  awk -v OFS="\t" '{ if ($NF<2) print $1, $4, $5, $9 }' | 
  sed 's/sampleName=//g'
done < ${GFFS} | sort -Vk1,1 -k2,2n -k3,3n -k4,4 > ${DEL_MASTER}
#Duplications
DUP_MASTER=`mktemp`
while read gff; do
  sed 's/\;/\t/g' ${gff} | fgrep -v "#" | sed 's/CN=CN//g' | \
  awk -v OFS="\t" '{ if ($NF>2) print $1, $4, $5, $9 }' | 
  sed 's/sampleName=//g'
done < ${GFFS} | sort -Vk1,1 -k2,2n -k3,3n -k4,4 > ${DUP_MASTER}

#Iterates over samples & merges across GFFs
while read ID; do
  #Make output directory per sample
  if ! [ -e ${OUTDIR}/${ID} ]; then
    mkdir ${OUTDIR}/${ID}
  fi
  #Merge deletions
  awk -v ID="${ID}" -v OFS="\t" '{ if ($4==ID) print $1, $2, $3 }' ${DEL_MASTER} | \
  bedtools merge -i - | \
  awk -v ID="${ID}" -v OFS="\t" '{ print $1, $2, $3, ID }' > ${OUTDIR}/${ID}/${ID}.cnMOPS.DEL.bed
  #Merge duplications
  awk -v ID="${ID}" -v OFS="\t" '{ if ($4==ID) print $1, $2, $3 }' ${DUP_MASTER} | \
  bedtools merge -i - | \
  awk -v ID="${ID}" -v OFS="\t" '{ print $1, $2, $3, ID }' > ${OUTDIR}/${ID}/${ID}.cnMOPS.DUP.bed
  #Subtracts intervals (if optioned)
  if [ ${SUBTRACT} != "0" ]; then
    for CNV in DEL DUP; do
      bedtools subtract \
      -a ${OUTDIR}/${ID}/${ID}.cnMOPS.${CNV}.bed \
      -b ${SUBTRACT} | sort -Vk1,1 -k2,2n -k3,3n > \
      ${OUTDIR}/${ID}/${ID}.cnMOPS.${CNV}.bed2 
      mv ${OUTDIR}/${ID}/${ID}.cnMOPS.${CNV}.bed2 \
      ${OUTDIR}/${ID}/${ID}.cnMOPS.${CNV}.bed
    done
  fi
  #Gzip (if optioned)
  if [ ${GZ} -eq 1 ]; then
    gzip -f ${OUTDIR}/${ID}/${ID}.cnMOPS.DEL.bed
    gzip -f ${OUTDIR}/${ID}/${ID}.cnMOPS.DUP.bed
  fi
done < ${SAMPLES}

#Clean up
#rm ${DEL_MASTER} ${DUP_MASTER}

