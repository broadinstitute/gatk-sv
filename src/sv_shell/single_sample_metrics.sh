#!/bin/bash

set -Exeuo pipefail

# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------

input_json=${1}
output_json_filename=${2:-""}
output_dir=${3:-""}

input_json="$(realpath ${input_json})"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d /output_single_sample_metrics_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir="$(mktemp -d /wd_single_sample_metrics_XXXXXXXX)"
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"


name=$(jq -r '.name' "$input_json")
case_sample=$(jq -r '.case_sample' "$input_json")
ref_samples=$(jq -r '.ref_samples[]' "$input_json")
sample_sr=$(jq -r '.sample_sr // ""' "$input_json")
sample_pe=$(jq -r '.sample_pe // ""' "$input_json")
sample_counts=$(jq -r '.sample_counts // ""' "$input_json")
cleaned_vcf=$(jq -r '.cleaned_vcf // ""' "$input_json")
final_vcf=$(jq -r '.final_vcf // ""' "$input_json")
contig_list=$(jq -r '.contig_list' "$input_json")
genotyped_pesr_vcf=$(jq -r '.genotyped_pesr_vcf // ""' "$input_json")
genotyped_depth_vcf=$(jq -r '.genotyped_depth_vcf // ""' "$input_json")
non_genotyped_unique_depth_calls_vcf=$(jq -r '.non_genotyped_unique_depth_calls_vcf // ""' "$input_json")
wgd_scores=$(jq -r '.wgd_scores // ""' "$input_json")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

samples_list="$(realpath "case_sample.txt")"
echo "${case_sample}" > "${samples_list}"

samples=("${case_sample}" ${ref_samples})
all_samples_list="$(realpath "all_samples.txt")"
printf "%s\n" "${samples[@]}" > "${all_samples_list}"

prefix="${case_sample}"

# SRMetrics
# ----------------------------------------------------------------------------------------------------------------------
SRMetrics_out=""
if [[ -n "${sample_sr}" ]]; then
  SRMetrics_out="$(realpath "${prefix}.sr-file.tsv")"
  svtest sr-file "${sample_sr}" "${samples_list}" > "${SRMetrics_out}"
fi


# PEMetrics
# ----------------------------------------------------------------------------------------------------------------------
PEMetrics_out=""
if [[ -n "${sample_pe}" ]]; then
  PEMetrics_out="$(realpath "${prefix}.pe-file.tsv")"
  svtest pe-file "${sample_pe}" "${samples_list}" > "${PEMetrics_out}"
fi


# CountsMetrics
# ----------------------------------------------------------------------------------------------------------------------
CountsMetrics_out=""
if [[ -n "${sample_counts}" ]]; then
  CountsMetrics_out="$(realpath "${case_sample}.raw-counts.tsv")"
  svtest raw-counts "${sample_counts}" "${case_sample}" > "${CountsMetrics_out}"
fi


# Cleaned_VCF_Metrics
# ----------------------------------------------------------------------------------------------------------------------
Cleaned_VCF_Metrics_out=""
if [[ -n "${cleaned_vcf}" ]]; then
  _prefix="cleaned"
  Cleaned_VCF_Metrics_out="$(realpath "${_prefix}.vcf.tsv")"

  svtest vcf \
    "${cleaned_vcf}" \
    "${contig_list}" \
    "${all_samples_list}" \
    "DEL,DUP,INS,INV,CTX,CNV,CPX,BND" \
    "${_prefix}" \
    > "${Cleaned_VCF_Metrics_out}"
fi


# Final_VCF_Metrics
# ----------------------------------------------------------------------------------------------------------------------
Final_VCF_Metrics_out=""
if [[ -n "${final_vcf}" ]]; then
  _prefix="final"
  Final_VCF_Metrics_out="$(realpath "${_prefix}.vcf.tsv")"

  svtest vcf \
    "${final_vcf}" \
    "${contig_list}" \
    "${all_samples_list}" \
    "DEL,DUP,INS,INV,CTX,CNV,CPX,BND" \
    "${_prefix}" \
    > "${Final_VCF_Metrics_out}"
fi


# Genotyped_PESR_VCF_Metrics
# ----------------------------------------------------------------------------------------------------------------------
Genotyped_PESR_VCF_Metrics_out=""
if [[ -n "${genotyped_pesr_vcf}" ]]; then
  _prefix="genotyped_pesr"
  Genotyped_PESR_VCF_Metrics_out="$(realpath "${_prefix}.vcf.tsv")"

  svtest vcf \
    "${genotyped_pesr_vcf}" \
    "${contig_list}" \
    "${all_samples_list}" \
    "DEL,DUP,INS,INV,BND" \
    "${_prefix}" \
    > "${Genotyped_PESR_VCF_Metrics_out}"
fi


# Genotyped_Depth_VCF_Metrics
# ----------------------------------------------------------------------------------------------------------------------
Genotyped_Depth_VCF_Metrics_out=""
if [[ -n "${genotyped_depth_vcf}" ]]; then
  _prefix="genotyped_depth"
  Genotyped_Depth_VCF_Metrics_out="$(realpath "${_prefix}.vcf.tsv")"

  svtest vcf \
    "${genotyped_depth_vcf}" \
    "${contig_list}" \
    "${all_samples_list}" \
    "DEL,DUP,CNV" \
    "${_prefix}" \
    > "${Genotyped_Depth_VCF_Metrics_out}"
fi


# NonGenotypedUniqueDepthCallsVCFMetrics
# ----------------------------------------------------------------------------------------------------------------------
NonGenotypedUniqueDepthCallsVCFMetrics_out=""
if [[ -n "${non_genotyped_unique_depth_calls_vcf}" ]]; then
  _prefix="non_genotyped_uniq_depth"
  NonGenotypedUniqueDepthCallsVCFMetrics_out="$(realpath "${_prefix}.vcf.tsv")"

  svtest vcf \
    "${non_genotyped_unique_depth_calls_vcf}" \
    "${contig_list}" \
    "${all_samples_list}" \
    "DEL,DUP" \
    "${_prefix}" \
    > "${NonGenotypedUniqueDepthCallsVCFMetrics_out}"
fi


# SingleSampleWGDMetrics
# ----------------------------------------------------------------------------------------------------------------------
SingleSampleWGDMetrics_out=""
if [[ -n "${wgd_scores}" ]]; then
  SingleSampleWGDMetrics_out="$(realpath "wgd.tsv")"
  SCORE=$(zcat "${wgd_scores}" | tail -n1 | cut -f2)
  echo "wgd_score_sample	${SCORE}" > "${SingleSampleWGDMetrics_out}"
fi


# CatMetrics
# ----------------------------------------------------------------------------------------------------------------------
CatMetrics_out="single_sample.${name}.metrics.tsv"

candidate_metric_files=(
  "${SRMetrics_out}"
  "${PEMetrics_out}"
  "${CountsMetrics_out}"
  "${Cleaned_VCF_Metrics_out}"
  "${Final_VCF_Metrics_out}"
  "${Genotyped_PESR_VCF_Metrics_out}"
  "${Genotyped_Depth_VCF_Metrics_out}"
  "${NonGenotypedUniqueDepthCallsVCFMetrics_out}"
  "${SingleSampleWGDMetrics_out}"
)

metric_files=()
for f in "${candidate_metric_files[@]}"; do
  if [[ -n "$f" ]]; then
    metric_files+=("$f")
  fi
done

cat "${metric_files[@]}" > "${CatMetrics_out}"
sed -i "s/${case_sample}/sample/g" "${CatMetrics_out}"



# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

CatMetrics_out_output="${output_dir}/$(basename "${CatMetrics_out}")"
mv "${CatMetrics_out}" "${CatMetrics_out_output}"

jq -n \
  --arg metrics_file "${CatMetrics_out_output}" \
  '{
     "metrics_file": $metrics_file
   }' > "${output_json_filename}"

echo "Finished single sample metrics successfully, output json filename: ${output_json_filename}"
