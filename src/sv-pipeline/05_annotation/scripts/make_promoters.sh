#!/bin/bash
#
# make_promoters.sh
#
# Derive promoters from gencode GTF
#

set -e

gtf=$1
window=$2

zcat $gtf \
  | awk -v FS="\t" -v OFS="\t" -v window=$window '{
      if ($3=="transcript") {
        split($9, fields, ";");
        split(fields[1], gene_id, " ");
        split(fields[5], gene_name, " ");
        if ($7=="+") {
          print $1, $4-window, $4, gene_name[2], ".", $7;
        } else {
          print $1, $5, $5+window, gene_name[2], ".", $7;
        }
      }}' \
  | sed -e 's/"//g'
