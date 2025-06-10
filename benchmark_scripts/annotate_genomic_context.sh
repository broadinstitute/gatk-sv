#!/bin/bash

# Copyright (c) 2024 Talkowski Laboratory
# Contact: Ryan Collins <xzhao12@mgh.harvard.edu>
# Distributed under terms of the MIT license.

# Benchmarks the sensitivity of one SV callset against another

set -e

########################################
# USAGE
########################################
usage() {
cat <<EOF

Usage: compare_callsets.sh [options] INPUT OUTPUT RepMask SegDup SimpRep

Helper tool to annotate genomic context of SVs.

Positional arguments:
  INPUT       SV callset to be benchmarked (BED format)
  OUTPUT      Path for output annotated callset
  RepMask     RepeatMasker regions (BED format)
  SegDup      Segmental duplications (BED format)
  SimpRep     Simple repeats (BED format)

Notes:
  - INPUT must be a BED3+ file with at least the following columns:
    [chr, start, end, SV ID, SV type, SV length]
  - OUTPUT is a filename for the annotated result.
  - Requires bedtools and Rscript with annotate_genomic_context_helper.R.

Options:
  -h, --help  Show this help message and exit

Example:
  ./compare_callsets.sh input.bed.gz annotated_output.bed.gz repmask.bed segdup.bed simple_repeats.bed

EOF
}

########################################
# PARSE ARGS
########################################
if [[ "$1" == "-h" || "$1" == "--help" || "$#" -ne 5 ]]; then
  usage
  exit 1
fi

INPUT=$1
OUTPUT=$2
RepMask=$3
SegDup=$4
SimpRep=$5

########################################
# SETUP
########################################
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
OVRTMP=$(mktemp -d)

########################################
# PROCESSING
########################################
zcat "${INPUT}" | cut -f1-6  > "${OVRTMP}/tmp.sites"

awk '{print $1,$2,$2,$4,$5,$6}' "${OVRTMP}/tmp.sites" | grep -v "#" | sed -e 's/ /\t/g' > "${OVRTMP}/il_inte.le_bp"
awk '{print $1,$3,$3,$4,$5,$6}' "${OVRTMP}/tmp.sites" | grep -v "#" | sed -e 's/ /\t/g' > "${OVRTMP}/il_inte.ri_bp"

zcat "${INPUT}" | awk '{if ($5=="DEL" || $5=="DUP" || $5=="CNV" ) print}' | awk '{if ($3-$2>5000) print}' | cut -f1-5 > "${OVRTMP}/il_inte.lg_cnv"

bedtools coverage -a "${OVRTMP}/il_inte.le_bp" -b "${SimpRep}" | awk '{if ($NF>0) print}' > "${OVRTMP}/il_inte.le_bp.vs.SR"
bedtools coverage -a "${OVRTMP}/il_inte.le_bp" -b "${SegDup}"  | awk '{if ($NF>0) print}' > "${OVRTMP}/il_inte.le_bp.vs.SD"
bedtools coverage -a "${OVRTMP}/il_inte.le_bp" -b "${RepMask}" | awk '{if ($NF>0) print}' > "${OVRTMP}/il_inte.le_bp.vs.RM"

bedtools coverage -a "${OVRTMP}/il_inte.ri_bp" -b "${SimpRep}" | awk '{if ($NF>0) print}' > "${OVRTMP}/il_inte.ri_bp.vs.SR"
bedtools coverage -a "${OVRTMP}/il_inte.ri_bp" -b "${SegDup}"  | awk '{if ($NF>0) print}' > "${OVRTMP}/il_inte.ri_bp.vs.SD"
bedtools coverage -a "${OVRTMP}/il_inte.ri_bp" -b "${RepMask}" | awk '{if ($NF>0) print}' > "${OVRTMP}/il_inte.ri_bp.vs.RM"

bedtools coverage -a "${OVRTMP}/il_inte.lg_cnv" -b "${SimpRep}" > "${OVRTMP}/il_inte.lg_cnv.vs.SR"
bedtools coverage -a "${OVRTMP}/il_inte.lg_cnv" -b "${SegDup}"  > "${OVRTMP}/il_inte.lg_cnv.vs.SD"
bedtools coverage -a "${OVRTMP}/il_inte.lg_cnv" -b "${RepMask}" > "${OVRTMP}/il_inte.lg_cnv.vs.RM"

Rscript "${BIN}/annotate_genomic_context_helper.R" -i "${INPUT}" -o "${OUTPUT}" -p "${OVRTMP}"
