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
  output_dir=$(mktemp -d /output_join_raw_calls_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_join_raw_calls_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

prefix=$(jq -r ".prefix" "${input_json}")

clustered_depth_vcf=$(jq -r '.clustered_depth_vcf // ""' "${input_json}")
clustered_depth_vcf_index=$(jq -r '.clustered_depth_vcf_index // ""' "${input_json}")

clustered_dragen_vcf=$(jq -r '.clustered_dragen_vcf // ""' "${input_json}")
clustered_dragen_vcf_index=$(jq -r '.clustered_dragen_vcf_index // ""' "${input_json}")

clustered_manta_vcf=$(jq -r '.clustered_manta_vcf // ""' "${input_json}")
clustered_manta_vcf_index=$(jq -r '.clustered_manta_vcf_index // ""' "${input_json}")

clustered_scramble_vcf=$(jq -r '.clustered_scramble_vcf // ""' "${input_json}")
clustered_scramble_vcf_index=$(jq -r '.clustered_scramble_vcf_index // ""' "${input_json}")

clustered_wham_vcf=$(jq -r '.clustered_wham_vcf // ""' "${input_json}")
clustered_wham_vcf_index=$(jq -r '.clustered_wham_vcf_index // ""' "${input_json}")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


# ConcatInputVcfs
# ---------------------------------------------------------------------------------------------------------------------

vcf_files=(
  "${clustered_depth_vcf}"
  "${clustered_dragen_vcf}"
  "${clustered_manta_vcf}"
  "${clustered_scramble_vcf}"
  "${clustered_wham_vcf}"
)

vcf_files_list="vcfs_list.txt"
touch "${vcf_files_list}"
for vcf_file_path in "${vcf_files[@]}"; do
  if [[ -n "${vcf_file_path}" ]]; then
    echo "${vcf_file_path}" >> "${vcf_files_list}"
  fi
done

concat_vcf="${prefix}.join_raw_calls.concat.vcf.gz"
concat_vcf_idx="${prefix}.join_raw_calls.concat.vcf.gz.tbi"
VCFS="~{write_lines(vcfs)}"

bcftools concat --no-version --allow-overlaps -Oz --file-list "${vcf_files_list}" > "${concat_vcf}"

tabix "${concat_vcf}"
