#!/bin/bash

set -Eeuo pipefail

# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------

input_json=${1}
output_json_filename=${2:-""}
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


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

make_bincov_matrix_input_json="$(realpath "${output_dir}/make_bincov_matrix_input.json")"
jq '
  {
    batch,
    samples,
    count_files,
    bincov_matrix_samples,
    bincov_matrix,
    reference_dict
  }
' "${input_json}" |
  jq '.' > "${make_bincov_matrix_input_json}"
make_bincov_matrix_output_json="$(realpath "${output_dir}/make_bincov_matrix_output.json")"
bash /make_bincov_matrix.sh "${make_bincov_matrix_input_json}" "${make_bincov_matrix_output_json}"
