#!/bin/bash

# ./gather_sample_evidence.sh test NA12878.final.cram NA12878.final.cram.crai Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.fai Homo_sapiens_assembly38.dict primary_contigs.list preprocessed_intervals.interval_list

set -Eeuo pipefail

sample_id=${1}
bam_or_cram_file=${2}
bam_or_cram_index=${3}
reference_fasta=${4}
reference_index=${5}
reference_dict=${6}
primary_contigs_list=${7}
preprocessed_intervals=${8}
disabled_read_filters=${9:-"MappingQualityReadFilter"}
collect_coverage=${10:-true}
run_scramble=${11:-true}


if [[ "${collect_coverage}" == true || "${run_scramble}" == true ]]; then
  # Collects read counts at specified intervals.
  # The count for each interval is calculated by counting the number of
  # read starts that lie in the interval.
  ./collect_counts.sh \
    "${preprocessed_intervals}" \
    "${bam_or_cram_file}" \
    "${bam_or_cram_index}" \
    "${sample_id}" \
    "${reference_fasta}" \
    "${reference_index}" \
    "${reference_dict}" \
    "/root/gatk.jar" \
    "${disabled_read_filters}"
fi

