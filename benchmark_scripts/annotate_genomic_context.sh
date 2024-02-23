#!/bin/bash

# Copyright (c) 2024 Talkowski Laboratory
# Contact: Ryan Collins <xzhao12@mgh.harvard.edu>
# Distributed under terms of the MIT license.

# Benchmarks the sensitivity of one SV callset against another

set -e

###USAGE
usage(){
cat <<EOF

usage: compare_callsets.sh [options] INPUT OUTPUT RepMask SegDup SimpRep

Helper tool to annotate genomic context of SVs

Positional arguments:
  INPUT            SV callset to be benchmarked
  OUTPUT            Ground truth SV callset to use for benchmarking
  RepMask     repeat masked regions in bed format
  SegDup      segmental duplicated regions in bed format
  SimpRep     simple repeated regions in bed format


Notes:
  1) INPUT are expected to be BED3+ formatted files
  2) SET1 must have at least the following five columns, in order:
     chr, start, end, SV ID, SV type, SV length
 
EOF
}


###PARSE ARGS

INPUT=$1
OUTPUT=$2
RepMask=$3
SegDup=$4
SimpRep=$5

###SET BIN
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


mkdir US_RM_SD_SR/
zcat ${INPUT} | cut -f1-6  > US_RM_SD_SR/tmp.sites
awk '{print $1,$2,$2,$4,$5,$6}'  US_RM_SD_SR/tmp.sites  | grep -v "#" | sed -e 's/ /\t/g' > US_RM_SD_SR/il_inte.le_bp
awk '{print $1,$3,$3,$4,$5,$6}'  US_RM_SD_SR/tmp.sites  | grep -v "#" | sed -e 's/ /\t/g' > US_RM_SD_SR/il_inte.ri_bp
zcat ${INPUT} | awk '{if ($5=="DEL" || $5=="DUP" || $5=="CNV" ) print}' | awk '{if ($3-$2>5000) print}' | cut -f1-5 >  US_RM_SD_SR/il_inte.lg_cnv

bedtools coverage -a US_RM_SD_SR/il_inte.le_bp -b ${SimpRep}  | awk '{if ($NF>0) print}'> US_RM_SD_SR/il_inte.le_bp.vs.SR
bedtools coverage -a US_RM_SD_SR/il_inte.le_bp -b ${SegDup}   | awk '{if ($NF>0) print}'> US_RM_SD_SR/il_inte.le_bp.vs.SD
bedtools coverage -a US_RM_SD_SR/il_inte.le_bp -b ${RepMask}  | awk '{if ($NF>0) print}'> US_RM_SD_SR/il_inte.le_bp.vs.RM

bedtools coverage -a US_RM_SD_SR/il_inte.ri_bp -b ${SimpRep}  | awk '{if ($NF>0) print}'> US_RM_SD_SR/il_inte.ri_bp.vs.SR
bedtools coverage -a US_RM_SD_SR/il_inte.ri_bp -b ${SegDup}   | awk '{if ($NF>0) print}'> US_RM_SD_SR/il_inte.ri_bp.vs.SD
bedtools coverage -a US_RM_SD_SR/il_inte.ri_bp -b ${RepMask}  | awk '{if ($NF>0) print}'> US_RM_SD_SR/il_inte.ri_bp.vs.RM

bedtools coverage -a US_RM_SD_SR/il_inte.lg_cnv -b ${SimpRep}  > US_RM_SD_SR/il_inte.lg_cnv.vs.SR
bedtools coverage -a US_RM_SD_SR/il_inte.lg_cnv -b ${SegDup}   > US_RM_SD_SR/il_inte.lg_cnv.vs.SD
bedtools coverage -a US_RM_SD_SR/il_inte.lg_cnv -b ${RepMask}  > US_RM_SD_SR/il_inte.lg_cnv.vs.RM

Rscript ${BIN}/annotate_genomic_context_helper.R -i ${INPUT} -o ${OUTPUT} -p "US_RM_SD_SR"

