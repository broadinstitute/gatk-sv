#!/bin/bash

# This script is a bash implementation of the following workflow/task in WDL:
# Filename: wdl/GatherSampleEvidence.wdl
# Workflow: GatherSampleEvidence

# Examples on running this script:
# bash gather_sample_evidence.sh test NA12878.final.cram NA12878.final.cram.crai Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.fai Homo_sapiens_assembly38.dict primary_contigs.list contig.fai preprocessed_intervals.interval_list primary_contigs_plus_mito.bed.gz primary_contigs_plus_mito.bed.gz Homo_sapiens_assembly38.dbsnp138.vcf hg38.repeatmasker.mei.with_SVA.pad_50_merged.bed.gz wham_whitelist.bed Homo_sapiens_assembly38.fasta.64.alt Homo_sapiens_assembly38.fasta.64.amb Homo_sapiens_assembly38.fasta.64.ann Homo_sapiens_assembly38.fasta.64.bwt Homo_sapiens_assembly38.fasta.64.pac Homo_sapiens_assembly38.fasta.64.sa
# example running using downsampled data.
# bash gather_sample_evidence.sh test downsampled_HG00096.final.cram downsampled_HG00096.final.cram.crai Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.fai Homo_sapiens_assembly38.dict downsampled_primary_contigs.list downsampled_contig.fai downsampled_preprocessed_intervals.interval_list downsampled_primary_contigs_plus_mito.bed.gz downsampled_primary_contigs_plus_mito.bed.gz downsampled_Homo_sapiens_assembly38.dbsnp138.vcf hg38.repeatmasker.mei.with_SVA.pad_50_merged.bed.gz downsampled_wham_whitelist.bed Homo_sapiens_assembly38.fasta.64.alt Homo_sapiens_assembly38.fasta.64.amb Homo_sapiens_assembly38.fasta.64.ann Homo_sapiens_assembly38.fasta.64.bwt Homo_sapiens_assembly38.fasta.64.pac Homo_sapiens_assembly38.fasta.64.sa

# Implementation notes:
# This script closely reproduces the GatherSampleEvidence workflow.
# However, there are few adjustments in the implementation and the pipeline
# to better suit the use case; for instance, we decided to skip running
# GatherSampleEvidenceMetrics module.

# For details: https://serverfault.com/a/103569
# saves original stdout to &3 and original stderr to &4
# without this, the logs of subprocess can get mixed with the parent's logs.
exec 3>&1 4>&2

set -Exeuo pipefail

RED='\033[0;31m'
BOLD_RED="\033[1;31m"
GREEN='\033[0;32m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color


CURRENT_STDERR_FILE="N/A"
# See https://stackoverflow.com/a/35800451 for details.
log_err() {
  local exit_code=$?
  local line_no=$1
  echo -e "${RED}Script on line ${line_no} exited with status ${exit_code}.${NC}" >&4
  echo -e "${BOLD_RED}Stdout: ${CURRENT_STDERR_FILE}${NC}" >&4
}
trap 'log_err $LINENO' ERR


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
preprocessed_intervals=${9}
manta_regions_bed=${10}
manta_regions_bed_index=${11}
sd_locs_vcf=${12}
mei_bed=${13}
include_bed_file=${14}
reference_bwa_alt=${15}
reference_bwa_amb=${16}
reference_bwa_ann=${17}
reference_bwa_bwt=${18}
reference_bwa_pac=${19}
reference_bwa_sa=${20}
disabled_read_filters=${21:-"MappingQualityReadFilter"}
collect_coverage=${22:-true}
run_scramble=${23:-true}
run_manta=${24:-true}
run_wham=${25:-true}
collect_pesr=${26:-true}
scramble_alignment_score_cutoff=${27:-90}
run_module_metrics=${28:-true}
min_size=${29:-50}
output_dir=${30:-""}


if [[ "${output_dir}" == "" ]]; then
  output_dir=$(mktemp -d output_gather_sample_evidence_XXXXXXXX)
  output_dir="$(realpath ${output_dir})"
fi

gather_sample_evidence_stdout="${output_dir}/gather_sample_evidence_stdout.txt"
gather_sample_evidence_stderr="${output_dir}/gather_sample_evidence_stderr.txt"
touch "${gather_sample_evidence_stdout}"
touch "${gather_sample_evidence_stderr}"

# The following directs all the trace output to a file descriptor 100 & also shows them on terminal.
# See https://stackoverflow.com/a/26611009
# If you want to skip showing trace on terminal, you may remove tee (copy-paste the solution from the above link).
exec 100> >(tee -a "${gather_sample_evidence_stdout}" >&2)
export BASH_XTRACEFD=100

echo -e "${MAGENTA}gather_sample_evidence.sh logs at: stdout:${gather_sample_evidence_stdout} and stderr:${gather_sample_evidence_stderr}${NC}"
gather_sample_evidence_start_time=`date +%s`


bam_or_cram_file="$(realpath ${bam_or_cram_file})"
bam_or_cram_index="$(realpath ${bam_or_cram_index})"
reference_fasta="$(realpath ${reference_fasta})"
reference_index="$(realpath ${reference_index})"
reference_dict="$(realpath ${reference_dict})"
primary_contigs_list="$(realpath ${primary_contigs_list})"
primary_contigs_fai="$(realpath ${primary_contigs_fai})"
preprocessed_intervals="$(realpath ${preprocessed_intervals})"
manta_regions_bed="$(realpath ${manta_regions_bed})"
manta_regions_bed_index="$(realpath ${manta_regions_bed_index})"
sd_locs_vcf="$(realpath ${sd_locs_vcf})"
mei_bed="$(realpath ${mei_bed})"
include_bed_file="$(realpath ${include_bed_file})"
reference_bwa_alt="$(realpath ${reference_bwa_alt})"
reference_bwa_amb="$(realpath ${reference_bwa_amb})"
reference_bwa_ann="$(realpath ${reference_bwa_ann})"
reference_bwa_bwt="$(realpath ${reference_bwa_bwt})"
reference_bwa_pac="$(realpath ${reference_bwa_pac})"
reference_bwa_sa="$(realpath ${reference_bwa_sa})"


if [[ "${collect_coverage}" == true || "${run_scramble}" == true ]]; then
  # Collects read counts at specified intervals.
  # The count for each interval is calculated by counting the number of
  # read starts that lie in the interval.

  collect_counts_outputs_json_filename=$(mktemp --suffix=.json "${output_dir}/collect_counts_XXXXXX")
  collect_counts_stdout=$(mktemp --suffix=.txt "${output_dir}/collect_counts_stdout_XXXXXX")
  collect_counts_stderr=$(mktemp --suffix=.txt "${output_dir}/collect_counts_stderr_XXXXXX")
  echo -e "${CYAN}Running collect_counts.sh ... stdout:${collect_counts_stdout} and stderr:${collect_counts_stderr}${NC}" | tee -a "${gather_sample_evidence_stdout}"
  collect_counts_start_time=`date +%s`

  CURRENT_STDERR_FILE="${collect_counts_stderr}"
  bash collect_counts.sh \
    "${preprocessed_intervals}" \
    "${bam_or_cram_file}" \
    "${bam_or_cram_index}" \
    "${sample_id}" \
    "${reference_fasta}" \
    "${reference_index}" \
    "${reference_dict}" \
    "${collect_counts_outputs_json_filename}" \
    "/root/gatk.jar" \
    "${disabled_read_filters}" > "${collect_counts_stdout}" 2> "${collect_counts_stderr}"

  collect_counts_end_time=`date +%s`
  collect_counts_et=$((collect_counts_end_time-collect_counts_start_time))
  echo -e "${GREEN}Successfully finished running collect_counts.sh in ${collect_counts_et} seconds.${NC}" | tee -a "${gather_sample_evidence_stdout}"
fi

if [[ "${run_manta}" == true ]]; then

  manta_outputs_json_filename=$(mktemp --suffix=.json "${output_dir}/manta_XXXXXX")
  manta_stdout=$(mktemp --suffix=.txt "${output_dir}/manta_stdout_XXXXXX")
  manta_stderr=$(mktemp --suffix=.txt "${output_dir}/manta_stderr_XXXXXX")
  echo -e "${CYAN}Running run_manta.sh ... stdout:${manta_stdout} and stderr:${manta_stderr}${NC}" | tee -a "${gather_sample_evidence_stdout}"
  manta_start_time=`date +%s`

  CURRENT_STDERR_FILE="${manta_stderr}"
  bash run_manta.sh \
    "${sample_id}" \
    "${bam_or_cram_file}" \
    "${bam_or_cram_index}" \
    "${reference_fasta}" \
    "${manta_regions_bed}" \
    "${manta_regions_bed_index}" \
    "${manta_outputs_json_filename}" > "${manta_stdout}" 2> "${manta_stderr}"

  manta_end_time=`date +%s`
  manta_et=$((manta_end_time-manta_start_time))
  echo -e "${GREEN}Successfully finished running run_manta.sh in ${manta_et} seconds.${NC}" | tee -a "${gather_sample_evidence_stdout}"
fi

if [[ "${collect_pesr}" == true ]]; then

  collect_pesr_outputs_json_filename=$(mktemp --suffix=.json "${output_dir}/collect_pesr_XXXXXX")
  collect_pesr_stdout=$(mktemp --suffix=.txt "${output_dir}/collect_pesr_stdout_XXXXXX")
  collect_pesr_stderr=$(mktemp --suffix=.txt "${output_dir}/collect_pesr_stderr_XXXXXX")
  echo -e "${CYAN}Running collect_sv_evidence.sh ... stdout:${collect_pesr_stdout} and stderr:${collect_pesr_stderr}${NC}" | tee -a "${gather_sample_evidence_stdout}"
  collect_pesr_start_time=`date +%s`

  CURRENT_STDERR_FILE="${collect_pesr_stderr}"
  bash collect_sv_evidence.sh \
    "${sample_id}" \
    "${bam_or_cram_file}" \
    "${bam_or_cram_index}" \
    "${reference_fasta}" \
    "${reference_index}" \
    "${reference_dict}" \
    "${sd_locs_vcf}" \
    "${preprocessed_intervals}" \
    "${collect_pesr_outputs_json_filename}" > "${collect_pesr_stdout}" 2> "${collect_pesr_stderr}"

  collect_pesr_end_time=`date +%s`
  collect_pesr_et=$((collect_pesr_end_time-collect_pesr_start_time))
  echo -e "${GREEN}Successfully finished running collect_sv_evidence.sh in ${collect_pesr_et} seconds.${NC}" | tee -a "${gather_sample_evidence_stdout}"
fi


if [[ "${run_scramble}" == true && "${run_manta}" == true ]]; then

  scramble_p1_outputs_json_filename=$(mktemp --suffix=.json "${output_dir}/scramble_p1_XXXXXX")
  scramble_stdout=$(mktemp --suffix=.txt "${output_dir}/scramble_stdout_XXXXXX")
  scramble_stderr=$(mktemp --suffix=.txt "${output_dir}/scramble_stderr_XXXXXX")
  echo -e "${CYAN}Running scramble.sh (part 1 & 2)... stdout:${scramble_stdout} and stderr:${scramble_stderr}${NC}" | tee -a "${gather_sample_evidence_stdout}"
  scramble_start_time=`date +%s`

  CURRENT_STDERR_FILE="${scramble_stderr}"

  {
    bash scramble.sh \
      "${sample_id}" \
      "${bam_or_cram_file}" \
      "${bam_or_cram_index}" \
      "${bam_or_cram_file}" \
      "${bam_or_cram_index}" \
      "$([ "${collect_pesr}" = "false" ] && echo "" || jq -r ".counts" "${collect_pesr_outputs_json_filename}")" \
      "$([ "${run_manta}" = "false" ] && echo "" || jq -r ".vcf" "${manta_outputs_json_filename}")" \
      "${reference_fasta}" \
      "${reference_index}" \
      "${primary_contigs_list}" \
      "${scramble_alignment_score_cutoff}" \
      "${mei_bed}" \
      "${scramble_p1_outputs_json_filename}"

    realign_soft_clipped_reads_json_filename=$(mktemp --suffix=.json "${output_dir}/realign_soft_clipped_reads_XXXXXX")
    # addresses bug in Dragen v3.7.8 where some reads are incorrectly soft-clipped

    bash realign_soft_clipped_reads.sh \
      "${sample_id}" \
      "${bam_or_cram_file}" \
      "${bam_or_cram_index}" \
      $(jq -r ".table" ${scramble_p1_outputs_json_filename}) \
      false \
      "${reference_fasta}" \
      "${reference_index}" \
      "${reference_bwa_alt}" \
      "${reference_bwa_amb}" \
      "${reference_bwa_ann}" \
      "${reference_bwa_bwt}" \
      "${reference_bwa_pac}" \
      "${reference_bwa_sa}" \
      "${realign_soft_clipped_reads_json_filename}"

    # ScrambleRealigned
    scramble_p2_outputs_json_filename=$(mktemp --suffix=.json "${output_dir}/scramble_p2_XXXXXX")
    bash scramble.sh \
      "${sample_id}" \
      $(jq -r ".out" ${realign_soft_clipped_reads_json_filename}) \
      $(jq -r ".out_index" ${realign_soft_clipped_reads_json_filename}) \
      "${bam_or_cram_file}" \
      "${bam_or_cram_index}" \
      "$([ "${collect_pesr}" = "false" ] && echo "" || jq -r ".counts" "${collect_pesr_outputs_json_filename}")" \
      "$([ "${run_manta}" = "false" ] && echo "" || jq -r ".vcf" "${manta_outputs_json_filename}")" \
      "${reference_fasta}" \
      "${reference_index}" \
      "${primary_contigs_list}" \
      "${scramble_alignment_score_cutoff}" \
      "${mei_bed}" \
      "${scramble_p1_outputs_json_filename}"
  } > "${scramble_stdout}" 2> "${scramble_stderr}"

  scramble_end_time=`date +%s`
  scramble_et=$((scramble_end_time-scramble_start_time))
  echo -e "${GREEN}Successfully finished running scramble.sh (part 1 & 2) in ${scramble_et} seconds.${NC}" | tee -a "${gather_sample_evidence_stdout}"
fi

if [[ "${run_wham}" == true ]]; then

  wham_outputs_json_filename=$(mktemp --suffix=.json "${output_dir}/wham_XXXXXX")
  wham_stdout=$(mktemp --suffix=.txt "${output_dir}/wham_stdout_XXXXXX")
  wham_stderr=$(mktemp --suffix=.txt "${output_dir}/wham_stderr_XXXXXX")
  echo -e "${CYAN}Running run_whamg.sh ... stdout:${wham_stdout} and stderr:${wham_stderr}${NC}" | tee -a "${gather_sample_evidence_stdout}"
  wham_start_time=`date +%s`

  CURRENT_STDERR_FILE="${wham_stderr}"
  bash run_whamg.sh \
    "${sample_id}" \
    "${bam_or_cram_file}" \
    "${bam_or_cram_index}" \
    "${reference_fasta}" \
    "${reference_index}" \
    "${include_bed_file}" \
    "${primary_contigs_list}" \
    "${wham_outputs_json_filename}" > "${wham_stdout}" 2> "${wham_stderr}"

  wham_end_time=`date +%s`
  wham_et=$((wham_end_time-wham_start_time))
  echo -e "${GREEN}Successfully finished running run_whamg.sh (part 1 & 2) in ${wham_et} seconds.${NC}" | tee -a "${gather_sample_evidence_stdout}"
fi

outputs_filename="${output_dir}/gather_sample_evidence_outputs.json"
outputs_json=$(jq -n \
  --arg coverage_counts "$([ "${collect_coverage}" = "false" ] && echo "" || jq -r ".counts" "${collect_counts_outputs_json_filename}")" \
  --arg manta_vcf "$([ "${run_manta}" = "false" ] && echo "" || jq -r ".vcf" "${manta_outputs_json_filename}")" \
  --arg manta_index "$([ "${run_manta}" = "false" ] && echo "" || jq -r ".index" "${manta_outputs_json_filename}")" \
  --arg scramble_vcf "$([ "${run_scramble}" = "false" ] && echo "" || jq -r ".vcf" "${scramble_p2_outputs_json_filename}")" \
  --arg scramble_index "$([ "${run_scramble}" = "false" ] && echo "" || jq -r ".index" "${scramble_p2_outputs_json_filename}")" \
  --arg scramble_clusters "$([ "${run_scramble}" = "false" ] && echo "" || jq -r ".clusters" "${scramble_p2_outputs_json_filename}")" \
  --arg scramble_table "$([ "${run_scramble}" = "false" ] && echo "" || jq -r ".table" "${scramble_p2_outputs_json_filename}")" \
  --arg pesr_disc "$([ "${collect_pesr}" = "false" ] && echo "" || jq -r ".disc_out" "${collect_pesr_outputs_json_filename}")" \
  --arg pesr_disc_index "$([ "${collect_pesr}" = "false" ] && echo "" || jq -r ".disc_out_index" "${collect_pesr_outputs_json_filename}")" \
  --arg pesr_split "$([ "${collect_pesr}" = "false" ] && echo "" || jq -r ".split_out" "${collect_pesr_outputs_json_filename}")" \
  --arg pesr_split_index "$([ "${collect_pesr}" = "false" ] && echo "" || jq -r ".split_out_index" "${collect_pesr_outputs_json_filename}")" \
  --arg pesr_sd "$([ "${collect_pesr}" = "false" ] && echo "" || jq -r ".sd_out" "${collect_pesr_outputs_json_filename}")" \
  --arg pesr_sd_index "$([ "${collect_pesr}" = "false" ] && echo "" || jq -r ".sd_out_index" "${collect_pesr_outputs_json_filename}")" \
  --arg wham_vcf "$([ "${run_wham}" = "false" ] && echo "" || jq -r ".vcf" "${wham_outputs_json_filename}")" \
  --arg wham_index "$([ "${run_wham}" = "false" ] && echo "" || jq -r ".index" "${wham_outputs_json_filename}")" \
  '{
     "coverage_counts": $coverage_counts,
     "manta_vcf": $manta_vcf,
     "manta_index": $manta_index,
     "scramble_vcf": $scramble_vcf,
     "scramble_index": $scramble_index,
     "scramble_clusters": $scramble_clusters,
     "scramble_table": $scramble_table,
     "pesr_disc": $pesr_disc,
     "pesr_disc_index": $pesr_disc_index,
     "pesr_split": $pesr_split,
     "pesr_split_index": $pesr_split_index,
     "pesr_sd": $pesr_sd,
     "pesr_sd_index": $pesr_sd_index,
     "wham_vcf": $wham_vcf,
     "wham_index": $wham_index
   }' \
)
echo "${outputs_json}" > "${outputs_filename}"

gather_sample_evidence_end_time=`date +%s`
gather_sample_evidence_et=$((gather_sample_evidence_end_time-gather_sample_evidence_start_time))
echo -e "${GREEN}Successfully finished running gather_sample_evidence in ${gather_sample_evidence_et} seconds. Outputs are serialized to: ${outputs_filename} ${NC}" | tee -a "${gather_sample_evidence_stdout}"
