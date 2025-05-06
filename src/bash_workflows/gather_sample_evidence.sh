#!/bin/bash

# chmod +x run_manta.sh collect_counts.sh collect_sv_evidence.sh gather_sample_evidence.sh scramble.sh run_whamg.sh standardize_vcf.sh
# ./gather_sample_evidence.sh test NA12878.final.cram NA12878.final.cram.crai Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.fai Homo_sapiens_assembly38.dict primary_contigs.list contig.fai preprocessed_intervals.interval_list primary_contigs_plus_mito.bed.gz primary_contigs_plus_mito.bed.gz Homo_sapiens_assembly38.dbsnp138.vcf hg38.repeatmasker.mei.with_SVA.pad_50_merged.bed.gz wham_whitelist.bed

set -Eeuo pipefail

sample_id=${1}
bam_or_cram_file=${2}
bam_or_cram_index=${3}
reference_fasta=${4}
reference_index=${5}
reference_dict=${6}
primary_contigs_list=${7}

# the wdl version sets this optional and requires it only if run_module_metrics is set,
# however conditional inputs like that are confusing and need additional check and docs.
# So, making it required here to keep the interface simpler.
primary_contigs_fai=${8}
preprocessed_intervals=${8}
manta_regions_bed=${9}
manta_regions_bed_index=${10}
sd_locs_vcf=${11}
mei_bed=${12}
include_bed_file=${13}
disabled_read_filters=${14:-"MappingQualityReadFilter"}
collect_coverage=${15:-true}
run_scramble=${16:-true}
run_manta=${17:-true}
run_wham=${18:-true}
collect_pesr=${19:-true}
scramble_alignment_score_cutoff=${20:-90}
run_module_metrics=${21:-true}
min_size=${22:-50}


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

if [[ "${run_wham}" == true ]]; then
  ./run_whamg.sh \
    "${sample_id}" \
    "${bam_or_cram_file}" \
    "${bam_or_cram_index}" \
    "${reference_fasta}" \
    "${reference_index}" \
    "${include_bed_file}" \
    "${primary_contigs_list}"
fi

if [[ "${run_module_metrics}" == true ]]; then
  if [[ "${run_manta}" == true ]]; then
    ./standardize_vcf.sh \
      "${sample_id}" \
      "/manta/${sample_id}.manta.vcf.gz" \
      "manta" \
      "${primary_contigs_fai}" \
      "${min_size}"
  fi
fi
