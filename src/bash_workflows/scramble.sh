#!/bin/bash

set -Eeuo pipefail

bam_or_cram_file=${1}
bam_or_cram_index=${2}
original_bam_or_cram_file=${3}
original_bam_or_cram_index=${4}
counts_file=${5}
manta_vcf=${6}
sample_name=${7}
reference_fasta=${8}
reference_index=${9}
regions_list=${10}

# Critical parameter for sensitivity/specificity
# Recommended values for aligners:
#   BWA-MEM: 90
#   DRAGEN-3.7.8: 60
alignment_score_cutoff=${11}

mei_bed=${12}

min_clipped_reads_fraction=${13:-0.22}
percent_align_cutoff=${14:-70}

#File? scramble_vcf_script
#String? make_scramble_vcf_args

part2_threads=${15:-7}


# ScramblePart1
# -------------

# Calibrate clipped reads cutoff based on median coverage
zcat "${counts_file}" \
  | awk '$0!~"@"' \
  | sed 1d \
  | awk 'NR % 100 == 0' \
  | cut -f4 \
  | Rscript -e "cat(round(${min_clipped_reads_fraction}*median(data.matrix(read.csv(file(\"stdin\"))))))" \
  > cutoff.txt
MIN_CLIPPED_READS=$(cat cutoff.txt)
echo "MIN_CLIPPED_READS: ${MIN_CLIPPED_READS}"

# Identify clusters of split reads
while read region; do
  time /app/scramble-gatk-sv/cluster_identifier/src/build/cluster_identifier -l -s ${MIN_CLIPPED_READS} -r "${region}" -t "${reference_fasta}" "${bam_or_cram_file}" \
    | gzip >> "${sample_name}".scramble_clusters.tsv.gz
done < "${regions_list}"