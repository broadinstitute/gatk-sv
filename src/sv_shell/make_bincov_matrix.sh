#!/bin/bash

# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------

set -Exeuo pipefail

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

input_json=${1}
output_json_filename=${2-""}
output_dir=${3:-""}

input_json="$(realpath ${input_json})"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/output_make_bincov_matrix_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_make_bincov_matrix_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

batch_name=$(jq -r ".batch" "${input_json}")

readarray -t count_files < <(jq -r '.count_files[]' "${input_json}")

bin_size=$(jq -r ".bin_size // 100" "${input_json}")
skip_bin_size_filter=$(jq -r ".skip_bin_size_filter // false" "${input_json}")

# These files need to have the '.list' extension (gatk requirement)
evidence_files_list="$(realpath "evidence_files.list")"
samples_filename="$(realpath "samples.list")"

reference_dict=$(jq -r ".reference_dict" "${input_json}")

jq -r ".samples[]" "${input_json}" >> "${samples_filename}"

# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

shopt -s nocasematch # for case-insensitive extension matching
for filename in "${count_files[@]}"; do
  case "${filename}" in
    *.counts.tsv.gz)
      sample_id="${filename%.counts.tsv.gz}"
      output_filename="$(realpath ${sample_id}.rd.txt.gz)"

      java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar \
        ConvertCountsToDepthFile \
        --counts-filename "${filename}" \
        --output "${output_filename}"

      echo "${output_filename}" >> "${evidence_files_list}"
      ;;

    *.rd.tsv.gz)
      echo "${filename}" >> "${evidence_files_list}"
      ;;

    *.rd.txt.gz)
      echo "${filename}" >> "${evidence_files_list}"
      ;;

    *)
      echo "${filename} extension does not match *.counts.tsv.gz or *.rd.tsv.gz or *.rd.txt.gz."
      ;;
  esac
done

merged_bincov_pre_filter="$(realpath ${batch_name}_pre_filter.RD.txt.gz)"
java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar \
  PrintSVEvidence \
    -F "${evidence_files_list}" \
    --sample-names "${samples_filename}" \
    --sequence-dictionary "${reference_dict}" \
    --output "${merged_bincov_pre_filter}"

# Note that the following removes bins that do not have
# the same bin size as the given bin_size variable.
# This is a quick patch, and ultimately we want to update
# the PrintSVEvidence tool to take bin_size and bin_locus
# input arguments.
merged_bincov="$(realpath ${batch_name}.RD.txt.gz)"
if [[ "${skip_bin_size_filter}" == "false" ]]; then
  zcat "${merged_bincov_pre_filter}" | \
    awk -v bin_size="${bin_size}" 'NR==1 || ($3 - $2) == bin_size' | \
    bgzip > "${merged_bincov}"

  tabix -p bed "${merged_bincov}"
else
  mv "${merged_bincov_pre_filter}" "${merged_bincov}"
  mv "${merged_bincov_pre_filter}.tbi" "${merged_bincov}.tbi"
fi

# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

merged_bincov_output="${output_dir}/$(basename "${merged_bincov}")"
mv "${merged_bincov}" "${merged_bincov_output}"

merged_bincov_idx_output="${output_dir}/$(basename "${merged_bincov}.tbi")"
mv "${merged_bincov}.tbi" "${merged_bincov_idx_output}"

jq -n \
  --arg merged_bincov "${merged_bincov_output}" \
  --arg merged_bincov_idx "${merged_bincov_idx_output}" \
  '{
      merged_bincov: $merged_bincov,
      merged_bincov_idx: $merged_bincov_idx
  }' > "${output_json_filename}"

echo "Finished make bincov matrix, output json filename: ${output_json_filename}"
