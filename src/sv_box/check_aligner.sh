#!/bin/bash

set -Exeuo pipefail

reads_path=$1
reads_index=$2
sample_id=$3
reference_fasta=${4:-""}
reference_index=${5:-""}
reference_dict=${6:-""}

if [[ -n "$reference_fasta" ]]; then
  ref_opt="-R $reference_fasta"
else
  ref_opt=""
fi

gatk PrintReadsHeader \
  -I "${reads_path}" \
  --read-index "${reads_index}" \
  -O "${sample_id}".header.sam \
  $ref_opt

awk '$0~"@PG" && $0~"ID: DRAGEN SW build" && $0~"VN: 05.021.604.3.7.8"' "${sample_id}".header.sam \
  | wc -l \
  > is_dragen_3_7_8.txt
