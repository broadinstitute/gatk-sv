#!/bin/bash

#./gather_sample_evidence.sh test NA12878.final.cram NA12878.final.cram.crai Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.fai Homo_sapiens_assembly38.dict primary_contigs.list preprocessed_intervals.interval_list primary_contigs_plus_mito.bed.gz primary_contigs_plus_mito.bed.gz Homo_sapiens_assembly38.dbsnp138.vcf hg38.repeatmasker.mei.with_SVA.pad_50_merged.bed.gz
#./gather_sample_evidence.sh test downsampled_HG00096.final.cram downsampled_HG00096.final.cram.crai Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.fai Homo_sapiens_assembly38.dict downsampled_primary_contigs.list downsampled_preprocessed_intervals.interval_list primary_contigs_plus_mito.bed.gz primary_contigs_plus_mito.bed.gz.tbi

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
sd_locs_vcf=${11}
mei_bed=${12}
disabled_read_filters=${13:-"MappingQualityReadFilter"}
collect_coverage=${14:-true}
run_scramble=${15:-true}
run_manta=${16:-true}
collect_pesr=${17:-true}
scramble_alignment_score_cutoff=${18:-90}


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
    "${manta_regions_bed}" \
    "${manta_regions_bed_index}"
fi

if [[ "${collect_pesr}" == true ]]; then
  ./collect_sv_evidence.sh \
    "${sample_id}" \
    "${bam_or_cram_file}" \
    "${bam_or_cram_index}" \
    "${reference_fasta}" \
    "${reference_index}" \
    "${reference_dict}" \
    "${sd_locs_vcf}" \
    "${preprocessed_intervals}"
fi


# TODO: change all the scripts to take the name of the outputs they are expected produce,
#       in the input, then change all the following to make sure such output names are used
#       consistently, for instance `"${sample_id}.counts.tsv.gz"` or manta vcf in the following.

if [[ "${run_scramble}" == true ]]; then
  ./scramble.sh \
    "${sample_id}" \
    "${bam_or_cram_file}" \
    "${bam_or_cram_index}" \
    "${bam_or_cram_file}" \
    "${bam_or_cram_index}" \
    "${sample_id}.counts.tsv.gz" \
    "/manta/${sample_id}.manta.vcf.gz" \
    "${reference_fasta}" \
    "${reference_index}" \
    "${primary_contigs_list}" \
    "${scramble_alignment_score_cutoff}" \
    "${mei_bed}"
fi
