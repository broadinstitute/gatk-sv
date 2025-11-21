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
  output_dir=$(mktemp -d /output_cnv_germline_case_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_cnv_germline_case_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

contig_ploidy_model_tar=$(jq -r ".contig_ploidy_model_tar" "${input_json}")
counts=($(jq -r ".counts[]" "${input_json}"))
ploidy_mapping_error_rate=$(jq -r '.ploidy_mapping_error_rate // "0.01"' "${input_json}")
ploidy_sample_psi_scale=$(jq -r '.ploidy_sample_psi_scale // "0.0001"' "${input_json}")


cpu_core_count="$(nproc)"

function getJavaMem() {
  # get JVM memory in MiB by getting total memory from /proc/meminfo
  # and multiplying by java_mem_fraction

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


source /opt/gatk_miniconda3/etc/profile.d/conda.sh
set +u
conda activate gatk
set -u

# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

# DetermineGermlineContigPloidyCaseMode
# ---------------------------------------------------------------------------------------------------------------------

DetermineGermlineContigPloidyCaseMode_wd=$(mktemp -d /wd_DetermineGermlineContigPloidyCaseMode_XXXXXXXX)
DetermineGermlineContigPloidyCaseMode_wd="$(realpath ${DetermineGermlineContigPloidyCaseMode_wd})"
cd "${DetermineGermlineContigPloidyCaseMode_wd}"


export MKL_NUM_THREADS="${cpu_core_count}"
export OMP_NUM_THREADS="${cpu_core_count}"

mkdir input-contig-ploidy-model
tar xzf "${contig_ploidy_model_tar}" -C input-contig-ploidy-model

read_count_files_list="counts_list.tsv"
printf "%s\n" "${counts[@]}" > "${read_count_files_list}"


grep 'gz$' "${read_count_files_list}" | while IFS= read -r filename; do
    output_filename=$(basename "$filename" .gz)
    zcat "$filename" > "$output_filename"
done
awk -F'/' '{ sub(/\.gz$/, "", $NF); print "--input " $NF }' "${read_count_files_list}" > read_count_files.args


java "-Xmx${JVM_MAX_MEM}" -jar /opt/gatk.jar DetermineGermlineContigPloidy \
  --arguments_file read_count_files.args \
  --model input-contig-ploidy-model \
  --output out \
  --output-prefix case \
  --verbosity DEBUG \
  --mapping-error-rate "${ploidy_mapping_error_rate}" \
  --sample-psi-scale "${ploidy_sample_psi_scale}"

tar c -C out/case-calls . | gzip -1 > case-contig-ploidy-calls.tar.gz


set +u
conda activate
set -u

echo "---------- finished"
