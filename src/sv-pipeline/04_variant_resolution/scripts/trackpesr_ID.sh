#!/bin/bash
#
# trackpesr_ID.sh
#

##Keep track of PE-SR IDs throughout pipeline##

set -euo pipefail

##gzipped vcf##
vcf=$1
originalidlist=$2
OUTFILE=$3

##append new ids to original list##
svtk vcf2bed $vcf int.bed -i MEMBERS

##remove header and match id one per line##
awk '{if (NR>1) print $4 "\t" $NF}'  int.bed \
  | awk -F'[,\t]' '{for(i=2; i<=NF; ++i) print $i "\t" $1 }' \
  | sort -k1,1\
  > newidlist.txt

join -j 1 -t $'\t' <(awk '{print $NF "\t" $0}' $originalidlist | sort -k1,1) newidlist.txt \
    | cut -f2-  \
    > ${OUTFILE}

