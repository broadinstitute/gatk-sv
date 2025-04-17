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
manta_regions_bed=${9}
manta_regions_bed_index=${10}
disabled_read_filters=${11:-"MappingQualityReadFilter"}
collect_coverage=${12:-true}
run_scramble=${13:-true}
run_manta=${14:-true}


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

if [[ "${run_manta}" == true ]]; then
  ./run_manta.sh \
    "${sample_id}" \
    "${bam_or_cram_file}" \
    "${bam_or_cram_index}" \
    "${reference_fasta}" \
    "${manta_regions_bed}"
fi

