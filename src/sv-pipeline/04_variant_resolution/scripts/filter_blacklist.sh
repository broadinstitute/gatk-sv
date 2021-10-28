#!/bin/bash
#
# filter_blacklist.sh
#
# Print IDs which are removable based on blacklist coverage
#

set -e

vcf=$1
blacklist=$2

bed=${vcf}.tmp.bed
svtk vcf2bed $vcf $bed

# 30% coverage
bedtools coverage -a $bed -b $blacklist \
  | awk '(($5=="DEL" || $5=="DUP" || $5=="CPX") && ($10 >= 0.3)) {print $4}'

# or direct overlap with start
awk -v OFS="\t" '($5=="INV" || $5~"INS") {$3=$2+1; print;}' $bed \
  | bedtools intersect -a stdin -b $blacklist \
  | cut -f4

# or direct overlap with end
awk -v OFS="\t" '($5=="INV" || $5~"INS") {end=$3; $2=end; $3=end+1; print;}' $bed \
  | bedtools intersect -a stdin -b $blacklist \
  | cut -f4

rm $bed
