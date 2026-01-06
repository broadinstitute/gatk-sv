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
  output_dir=$(mktemp -d /output_sv_concordance_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir="$(mktemp -d /wd_sv_concordance_XXXXXXXX)"
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

reference_dict=$(jq -r ".reference_dict" "$input_json")
eval_vcf=$(jq -r ".eval_vcf" "$input_json")
truth_vcf=$(jq -r ".truth_vcf" "$input_json")
output_prefix=$(jq -r ".output_prefix" "$input_json")
java_mem_fraction=$(jq -r '.java_mem_fraction // ""' "$input_json")


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

SVConcordanceTask_out="${output_prefix}.concordance.vcf.gz"
SVConcordanceTask_out_index="${output_prefix}.concordance.vcf.gz.tbi"

java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar SVConcordance \
  --sequence-dictionary "${reference_dict}" \
  --eval "${eval_vcf}" \
  --truth "${truth_vcf}" \
  -O "${SVConcordanceTask_out}"


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------


SVConcordanceTask_out_output_dir="${output_dir}/$(basename "${SVConcordanceTask_out}")"
mv "${SVConcordanceTask_out}" "${SVConcordanceTask_out_output_dir}"
mv "${SVConcordanceTask_out}.tbi" "${SVConcordanceTask_out_output_dir}.tbi"


outputs_json=$(jq -n \
  --arg concordance_vcf "${SVConcordanceTask_out_output_dir}" \
  --arg concordance_vcf_index "${SVConcordanceTask_out_output_dir}.tbi" \
  '{
     "concordance_vcf": $concordance_vcf,
     "concordance_vcf_index": $concordance_vcf_index
   }' \
)
echo "${outputs_json}" > "${output_json_filename}"

echo "Successfully finished SVConcordance, output json filename: ${output_json_filename}"