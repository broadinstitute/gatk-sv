#!/bin/bash

set -Eeuo pipefail

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

zcat "${scramble_table}" \
  | sed 1d \
  | cut -f1 \
  | tr ':' '\t' \
  | awk -F'\t' -v OFS='\t' '{print $1,$2,$2+1}' \
  | sort -k1,1V -k2,2n \
  | bedtools slop -i - -g "${reference_index}" -b 150 \
  | bedtools merge \
  > intervals.bed

TMPDIR=`mktemp -d -p .` || exit 1

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
  | samtools sort -T "${TMPDIR}" \
  | samtools view -1 -h -O BAM -o "${sample_id}.realign_soft_clipped_reads.bam"
samtools index -@${N_CORES} "${sample_id}.realign_soft_clipped_reads.bam"
