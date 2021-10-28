#!/bin/bash
#
# merge.sh
#

set -e

quad=$1

ALU=${quad}.ALU.final_comp.vcf
LINE1=${quad}.LINE1.final_comp.vcf
SVA=${quad}.SVA.final_comp.vcf

vcf-concat -c $ALU $LINE1 $SVA |& fgrep -e "do not match"

if [[ $? -eq 0 ]]; then
  echo "Columns do not match for quad ${quad}"
  exit 1
fi

cat merged_header.vcf \
  <(egrep -v -e "^##" $ALU | sed -e 's/.final//g' | sed -f ~/sfari/lists/ID_swap.sed) \
  <(egrep -v -e "^#" $LINE1) \
  <(egrep -v -e "^#" $SVA) \
  | vcf-sort -c \
  | bgzip -c \
  > melt.${quad}.vcf.gz

tabix -p vcf melt.${quad}.vcf.gz
