#!/bin/bash

# This script is a bash implementation of the following workflow/task in WDL:
# Filename: wdl/GatherSampleEvidence.wdl
# Task: RealignSoftClippedReads

set -Exeuo pipefail

sample_id=${1}
reads_path=${2}
reads_index=${3}
scramble_table=${4}
is_bam=${5}
reference_fasta=${6}
reference_index=${7}
reference_bwa_alt=${8}
reference_bwa_amb=${9}
reference_bwa_ann=${10}
reference_bwa_bwt=${11}
reference_bwa_pac=${12}
reference_bwa_sa=${13}
outputs_json_filename=${14}

working_dir=$(mktemp -d wd_realign_soft_clipped_reads_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
output_dir=$(mktemp -d output_realign_soft_clipped_reads_XXXXXXXX)
output_dir="$(realpath ${output_dir})"
cd "${working_dir}"

zcat "${scramble_table}" \
  | sed 1d \
  | cut -f1 \
  | tr ':' '\t' \
  | awk -F'\t' -v OFS='\t' '{print $1,$2,$2+1}' \
  | sort -k1,1V -k2,2n \
  | bedtools slop -i - -g "${reference_index}" -b 150 \
  | bedtools merge \
  > intervals.bed

samtools view --header-only "${reads_path}" > header.sam
N_CORES=$(nproc)
time samtools view --no-header \
  -T "${reference_fasta}" \
  -ML intervals.bed \
  "${reads_path}" \
  | awk -F'\t' -v OFS='\t' '$6~"S"' \
  | sort -u \
  | cat header.sam - \
  | samtools fastq \
  > reads.fastq
bwa mem -H header.sam -K 100000000 -v 3 -t ${N_CORES} -Y "${reference_fasta}" reads.fastq \
  | samtools sort -T "${working_dir}" \
  | samtools view -1 -h -O BAM -o "${sample_id}.realign_soft_clipped_reads.bam"
samtools index "-@${N_CORES}" "${sample_id}.realign_soft_clipped_reads.bam"


out_filename="${output_dir}/${sample_id}.realign_soft_clipped_reads.bam"
out_index_filename="${output_dir}/${sample_id}.realign_soft_clipped_reads.bam.bai"

mv "${sample_id}.realign_soft_clipped_reads.bam" "${out_filename}"
mv "${sample_id}.realign_soft_clipped_reads.bam.bai" "${out_index_filename}"

outputs_filename="${output_dir}/outputs.json"
outputs_json=$(jq -n \
  --arg ofname "${out_filename}" \
  --arg oifname "${out_index_filename}" \
  '{out: $ofname, out_index: $oifname}' )
echo "${outputs_json}" > "${outputs_filename}"
cp "${outputs_filename}" "${outputs_json_filename}"
