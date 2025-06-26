#!/bin/bash

# Copyright (c) 2024 Talkowski Laboratory
# Contact: Ryan Collins <xzhao12@mgh.harvard.edu>
# Distributed under terms of the MIT license.

# Benchmarks the sensitivity of one SV callset against another

set -e

### USAGE
usage(){
cat <<EOF

usage: compare_callsets.sh INPUT OUTPUT [--rm RepMask] [--sd SegDup] [--sr SimpRep]

Helper tool to annotate genomic context of SVs

Positional arguments:
  INPUT            SV callset to be benchmarked
  OUTPUT           Ground truth SV callset to use for benchmarking

Optional arguments:
  --rm FILE        RepeatMasked regions in BED format
  --sd FILE        Segmental duplicated regions in BED format
  --sr FILE        Simple repeated regions in BED format
  -h               Show this help message and exit

Notes:
  1) INPUT and OUTPUT are expected to be BED3+ formatted files
  2) INPUT must have at least the following five columns, in order:
     chr, start, end, SV ID, SV type, SV length

EOF
}

# Parse positional args
if [[ "$1" == "-h" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 2 ]]; then
  echo "ERROR: Missing required INPUT and OUTPUT arguments."
  usage
  exit 1
fi

INPUT=$1
OUTPUT=$2
shift 2

# Set defaults
RepMask=""
SegDup=""
SimpRep=""

# Parse optional arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --rm)
      RepMask=$2
      shift 2
      ;;
    --sd)
      SegDup=$2
      shift 2
      ;;
    --sr)
      SimpRep=$2
      shift 2
      ;;
    -h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

# You can now use $INPUT, $OUTPUT, $RepMask, $SegDup, and $SimpRep in the rest of your script

echo "INPUT: $INPUT"
echo "OUTPUT: $OUTPUT"
echo "RepMask: $RepMask"
echo "SegDup: $SegDup"
echo "SimpRep: $SimpRep"

###SET BIN
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

is_bgzf() {
    local file="$1"
    if file "$file" | grep -q "BGZF"; then
        echo "true"
    else
        echo "false"
    fi
}


TMPPATH=`mktemp -d`

if [[ $(is_bgzf "$1") == "true" ]]; then
    echo "The input file is bgzipped."
    zcat ${INPUT} | cut -f1-6  > ${TMPPATH}/tmp.sites
    zcat ${INPUT} | awk '{if ($5=="DEL" || $5=="DUP" || $5=="CNV" ) print}' | awk '{if ($3-$2>5000 || $3-$2< -5000) print}' | cut -f1-5 >  ${TMPPATH}/il_inte.lg_cnv
else
    echo "The input file is not bgzipped."
    cut -f1-6 ${INPUT} > ${TMPPATH}/tmp.sites
    awk '{if ($5=="DEL" || $5=="DUP" || $5=="CNV" ) print}' ${INPUT} | awk '{if ($3-$2>5000 || $3-$2< -5000) print}' | cut -f1-5 >  ${TMPPATH}/il_inte.lg_cnv
fi

#zcat ${INPUT} | cut -f1-6  > ${TMPPATH}/tmp.sites
awk '{print $1,$2,$2,$4,$5,$6}'  ${TMPPATH}/tmp.sites  | grep -v "#" | sed -e 's/ /\t/g' > ${TMPPATH}/il_inte.le_bp
awk '{print $1,$3,$3,$4,$5,$6}'  ${TMPPATH}/tmp.sites  | grep -v "#" | sed -e 's/ /\t/g' > ${TMPPATH}/il_inte.ri_bp



bedtools coverage -a ${TMPPATH}/il_inte.le_bp -b ${SimpRep}  | awk '{if ($NF>0) print}'> ${TMPPATH}/il_inte.le_bp.vs.SR
bedtools coverage -a ${TMPPATH}/il_inte.le_bp -b ${SegDup}   | awk '{if ($NF>0) print}'> ${TMPPATH}/il_inte.le_bp.vs.SD
bedtools coverage -a ${TMPPATH}/il_inte.le_bp -b ${RepMask}  | awk '{if ($NF>0) print}'> ${TMPPATH}/il_inte.le_bp.vs.RM

bedtools coverage -a ${TMPPATH}/il_inte.ri_bp -b ${SimpRep}  | awk '{if ($NF>0) print}'> ${TMPPATH}/il_inte.ri_bp.vs.SR
bedtools coverage -a ${TMPPATH}/il_inte.ri_bp -b ${SegDup}   | awk '{if ($NF>0) print}'> ${TMPPATH}/il_inte.ri_bp.vs.SD
bedtools coverage -a ${TMPPATH}/il_inte.ri_bp -b ${RepMask}  | awk '{if ($NF>0) print}'> ${TMPPATH}/il_inte.ri_bp.vs.RM

if [ "$(wc -l < ${TMPPATH}/il_inte.lg_cnv)" -gt 1 ]; then
	bedtools coverage -a ${TMPPATH}/il_inte.lg_cnv -b ${SimpRep}  > ${TMPPATH}/il_inte.lg_cnv.vs.SR
	bedtools coverage -a ${TMPPATH}/il_inte.lg_cnv -b ${SegDup}   > ${TMPPATH}/il_inte.lg_cnv.vs.SD
	bedtools coverage -a ${TMPPATH}/il_inte.lg_cnv -b ${RepMask}  > ${TMPPATH}/il_inte.lg_cnv.vs.RM
fi


Rscript ${BIN}/annotate_genomic_context_helper.R -i ${INPUT} -o ${OUTPUT} -p "${TMPPATH}"

