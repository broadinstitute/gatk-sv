#!/bin/bash
#
# localize_bincov.sh
#

set -e

bed=$1
bincov_gs_path=$2

sort -k1,1V -k2,2n $bed | bedtools merge -i stdin -d 1000000 > merged.bed
export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

while read chrom start end; do 
  echo "Fetching ${chrom}:${start}-${end}"
  tabix -h ${bincov_gs_path} "$chrom":"$start"-"$end" \
    | sed -e 's/Chr/chr/g' -e 's/Start/start/g' -e 's/End/end/' \
    | bgzip -c \
    > local_coverage.bed.gz
done < <(head -n1 merged.bed)

while read chrom start end; do 
  echo "Fetching ${chrom}:${start}-${end}"
  tabix -h ${bincov_gs_path} "$chrom":"$start"-"$end" \
    | sed '1d' \
    | bgzip -c \
    >> local_coverage.bed.gz
done < <(sed '1d' merged.bed)

tabix -p bed local_coverage.bed.gz

