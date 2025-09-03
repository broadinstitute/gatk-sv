#!/bin/bash

# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------

set -Eeuo pipefail

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
  output_dir=$(mktemp -d output_make_bincov_matrix_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d wd_make_bincov_matrix_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

batch_name=$(jq -r ".batch" "${input_json}")

readarray -t count_files < <(jq -r '.count_files[]' "${input_json}")

# These files need to have the `.list` extension (gatk requirement)
evidence_files_list="evidence_files.list"
samples_filename="samples.list"

reference_dict=$(jq -r ".reference_dict" "${input_json}")

bincov_matrix=$(jq -r ".bincov_matrix" "${input_json}")

jq -r ".samples[]" "${input_json}" >> "${samples_filename}"
# TODO: adding the folloiwng results in getting values for all the samples in the bincov matrix,
#  and the output is generates does not match the current output (it mostly contains -1
#  for all the columns except the sample in samples list).
#jq -r ".bincov_matrix_samples[]" "${input_json}" >> "${samples_filename}"

# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

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
      echo filename >> "${evidence_files_list}"
      ;;

    *)
      echo "File extension does not match *.counts.tsv.gz or *.rd.tsv.gz."
      ;;
  esac
done

merged_bincov="$(realpath ${batch_name}.RD.txt.gz)"
java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar \
  PrintSVEvidence \
    -F "${evidence_files_list}" \
    --sample-names "${samples_filename}" \
    --sequence-dictionary "${reference_dict}" \
    --output "${merged_bincov}"


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

merged_bincov_output="${output_dir}/$(basename "${merged_bincov}")"
mv "${merged_bincov}" "${merged_bincov_output}"

merged_bincov_idx_output="${output_dir}/$(basename "${merged_bincov}.tbi")"
mv "${merged_bincov}.tbi" "${merged_bincov_idx_output}"

outputs_json=$(jq -n \
  --arg merged_bincov "${merged_bincov_output}" \
  --arg merged_bincov_idx "${merged_bincov_idx_output}" \
  '{merged_bincov: $merged_bincov, merged_bincov_idx: $merged_bincov_idx}' )
echo "${outputs_json}" > "${output_json_filename}"
