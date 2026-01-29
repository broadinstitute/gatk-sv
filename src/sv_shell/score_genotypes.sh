#!/bin/bash

set -Exeuo pipefail


# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------


input_json=${1}
output_json_filename=${2-""}
output_dir=${3:-""}

input_json="$(realpath ${input_json})"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/output_score_genotypes_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_score_genotypes_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
echo "Filter Genotypes working directory: ${working_dir}"

vcf=$(jq -r ".vcf" "$input_json")
genome_tracks=($(jq -r ".genome_tracks[]" "$input_json"))
gq_recalibrator_model_file=$(jq -r ".gq_recalibrator_model_file" "$input_json")
recalibrate_gq_args=($(jq -r ".recalibrate_gq_args[]" "$input_json"))


function getJavaMem() {
  # get JVM memory in MiB by getting total memory from /proc/meminfo
  # and multiplying by java_mem_fraction

  local mem_fraction=${java_mem_fraction:=0.85}
  cat /proc/meminfo | \
    awk -v MEM_FIELD="$1" -v frac="${mem_fraction}" '{
      f[substr($1, 1, length($1)-1)] = $2
    } END {
      printf "%dM", f[MEM_FIELD] * frac / 1024
    }'
}
JVM_MAX_MEM=$(getJavaMem MemTotal)
echo "JVM memory: $JVM_MAX_MEM"


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

tabix -f "${vcf}"

base_name=${vcf%.vcf.gz}
base_name=${base_name%.vcf}
filtered_vcf_name="${base_name}_gq_recalibrated.vcf.gz"
filtered_vcf_name="$(realpath "${filtered_vcf_name}")"

genome_tracks_arg=""
if [ ${#genome_tracks[@]} -gt 0 ]; then
  for track in "${genome_tracks[@]}"; do
    genome_tracks_arg+=" --genome-track ${track}"
  done
fi

recalibrate_gq_args_combined=""
if [ ${#recalibrate_gq_args[@]} -gt 0 ]; then
  for gq_arg in "${recalibrate_gq_args[@]}"; do
    recalibrate_gq_args_combined+="${gq_arg} "
  done
fi

java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar XGBoostMinGqVariantFilter \
  --mode "Filter" \
  --variant "${vcf}" \
  ${genome_tracks_arg} \
  --model-file "${gq_recalibrator_model_file}" \
  --output "${filtered_vcf_name}" \
  ${recalibrate_gq_args_combined}

# gatk indices still have problems, overwrite with tabix
tabix -f "${filtered_vcf_name}"


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

filtered_vcf_name_output_dir="${output_dir}/$(basename "${filtered_vcf_name}")"
mv "${filtered_vcf_name}" "${filtered_vcf_name_output_dir}"
mv "${filtered_vcf_name}.tbi" "${filtered_vcf_name_output_dir}.tbi"


jq -n \
  --arg unfiltered_recalibrated_vcf "${filtered_vcf_name_output_dir}" \
  --arg unfiltered_recalibrated_vcf_index "${filtered_vcf_name_output_dir}.tbi" \
  '{
      unfiltered_recalibrated_vcf: $unfiltered_recalibrated_vcf,
      unfiltered_recalibrated_vcf_index: $unfiltered_recalibrated_vcf_index
  }' > "${output_json_filename}"

echo "Finished Score Genotypes successfully, output json filename: ${output_json_filename}"
