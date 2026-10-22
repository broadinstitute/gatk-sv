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
  output_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/output_genotype_svs_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(realpath $(mktemp -d "${SV_SHELL_BASE_DIR}/wd_genotype_svs_XXXXXXXX"))
cd "${working_dir}"
echo "Genotype SVs Working directory: ${working_dir}"

reference_dict=$(jq -r ".reference_dict" "${input_json}")
output_prefix=$(jq -r ".output_prefix" "${input_json}")


function getJavaMem() {
  local mem_fraction=${java_mem_fraction:=0.6}
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


java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar PrintSVEvidence \
  --sequence-dictionary "${reference_dict}" \
  --evidence-file $(jq -r ".rd_file" "${input_json}") \
  -O local.rd.txt.gz

java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar PrintSVEvidence \
  --sequence-dictionary "${reference_dict}" \
  --evidence-file $(jq -r ".pe_file" "${input_json}") \
  -O local.pe.txt.gz

java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar PrintSVEvidence \
  --sequence-dictionary "${reference_dict}" \
  --evidence-file $(jq -r ".sr_file" "${input_json}") \
  -O local.sr.txt.gz

output_filename="${output_prefix}.vcf.gz"
java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar GenotypeSVs \
  -V $(jq -r ".vcf" "${input_json}") \
  -O "${output_filename}" \
  --median-coverage $(jq -r ".median_coverage" "${input_json}") \
  --rd-file local.rd.txt.gz \
  --discordant-pairs-file local.pe.txt.gz \
  --split-reads-file local.sr.txt.gz \
  --sequence-dictionary "${reference_dict}" \
  --ploidy-table $(jq -r ".ploidy_table" "${input_json}") \
  --pesr-exclusion-intervals $(jq -r ".pesr_exclusion_intervals" "${input_json}") \
  --depth-exclusion-intervals $(jq -r ".depth_exclusion_intervals" "${input_json}") \
  --rd-table $(jq -r ".rd_table" "${input_json}") \
  --pe-table $(jq -r ".pe_table" "${input_json}") \
  --sr-table $(jq -r ".sr_table" "${input_json}")


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

output_filename_output_dir="$(realpath "${output_dir}/$(basename "${output_filename}")")"
mv "${output_filename}" "${output_filename_output_dir}"
mv "${output_filename}.tbi" "${output_filename_output_dir}.tbi"

jq -n \
  --arg out "${output_filename_output_dir}" \
  --arg out_index "${output_filename_output_dir}.tbi" \
  '{
      out: $out,
      out_index: $out_index
  }' > "${output_json_filename}"

echo "Finished genotype SVs: ${output_json_filename}"