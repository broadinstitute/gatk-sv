#!/bin/bash

# Cleans fam file to only include samples present in an input list

set -e

###USAGE
usage(){
cat <<EOF

usage: cleanFamFile.sh [-h] SAMPLES FAMFILE OUTFILE

Clean a fam file to only include samples present in an input list

Positional arguments:
  SAMPLES         List of samples used during filtering of fam file
  FAMFILE         Fam file to be filtered
  OUTFILE         Output path for filtered fam file

Optional arguments:
  -h  HELP        Show this help message and exit

EOF
}


###PARSE ARGS
SAMPLES=$1
FAMFILE=$2
OUTFILE=$3


###PROCESS ARGS & OPTIONS
#Check for required input
if [ -z ${SAMPLES} ]; then
  echo -e "\nERROR: input samples list not specified\n"
  usage
  exit 0
fi
if ! [ -s ${SAMPLES} ]; then
  echo -e "\nERROR: input samples list either empty or not found\n"
  usage
  exit 0
fi
if [ -z ${FAMFILE} ]; then
  echo -e "\nERROR: input fam file not specified\n"
  usage
  exit 0
fi
if ! [ -s ${FAMFILE} ]; then
  echo -e "\nERROR: input fam file either empty or not found\n"
  usage
  exit 0
fi


###FILTER FAM FILE
#Write header
echo -e "#FAM_ID\tPROBAND\tFATHER\tMOTHER\tPROBAND_SEX\tPHENOTYPE" > ${OUTFILE}
#Write revised list of duos and trios after masking on analysis_samples.list
while read fam pro fa mo sex pheno; do
  #Only consider families with child in analysis set
  if [ $( fgrep -w ${pro} ${SAMPLES} | wc -l ) -gt 0 ]; then
    #Check for father
    if [ $( fgrep -w ${fa} ${SAMPLES} | wc -l ) -eq 0 ]; then
      fa="."
    fi
    
    #Check for mother
    if [ $( fgrep -w ${mo} ${SAMPLES} | wc -l ) -eq 0 ]; then
      mo="."
    fi

    #Only consider families with at least one parent
    if [ ${fa} != "." ] || [ ${mo} != "." ]; then
      #Write as duo if one parent is missing
      if [ ${fa} == "." ] || [ ${mo} == "." ]; then
        echo -e "DUO_\t${pro}\t${fa}\t${mo}\t${sex}\t${pheno}"
      #Otherwise, write as trio
      else
        echo -e "TRIO_\t${pro}\t${fa}\t${mo}\t${sex}\t${pheno}"
      fi
    fi
  fi
done < <( fgrep -v "#" ${FAMFILE} ) | sort -Vk1,1 -k2,2V > ${OUTFILE}.tmp
#Clean output fam file
for f in DUO TRIO; do
  grep -e "^${f}" ${OUTFILE}.tmp | \
  sort -Vk2,2 | awk -v OFS="\t" '{ print $1NR, $2, $3, $4, $5, $6 }'
done >> ${OUTFILE}
rm ${OUTFILE}.tmp

