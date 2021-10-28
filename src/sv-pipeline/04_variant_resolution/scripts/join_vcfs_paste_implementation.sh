#!/bin/bash
#
# join_vcfs.sh
#
# Join VCFs which contain the same set of variants in different samples 
#

set -e

vcflist=$1
prefix=$2

fout=${prefix}.vcf

base=$(head -n1 $vcflist)

if [[ ${base} == *.gz ]]; then
  zcat $base | sed -n -e '/^##/p' > $fout
else
  sed -n -e '/^##/p' $base > $fout
fi

# Paste final line of header (samples) and VCF records together
if [[ ${base} == *.gz ]]; then
  cmd="paste <(zcat $base | sed -e '/^##/d')"
else
  cmd="paste <(sed -e '/^##/d' $base)"
fi

while read vcf; do
  if [[ ${vcf} == *.gz ]]; then
    cmd="$cmd <(zcat $vcf | sed -e '/^##/d' | cut -f 10-)"
  else
    cmd="$cmd <(sed -e '/^##/d' $vcf | cut -f 10-)"
  fi
done < <(sed '1d' $vcflist)

eval "$cmd" >> $fout

bgzip -f $fout