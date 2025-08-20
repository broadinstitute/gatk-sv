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
  output_dir=$(mktemp -d output_evidence_qc_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d wd_evidence_qc_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

batch=$(jq -r ".batch" "${input_json}")
bincov_matrix=$(jq -r ".bincov_matrix" "${input_json}")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

# --- MakeBincovMatrix
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


# --- MedianCov
merged_bincov=$(jq -r ".merged_bincov" "${make_bincov_matrix_output_json}")
zcat "${merged_bincov}" > "${batch}_fixed.bed"
Rscript /opt/WGD/bin/medianCoverage.R "${batch}_fixed.bed" -H "${batch}_medianCov.bed"
Rscript -e "x <- read.table(\"${batch}_medianCov.bed\",check.names=FALSE); xtransposed <- t(x[,c(1,2)]); write.table(xtransposed,file=\"${batch}_medianCov.transposed.bed\",sep=\"\\t\",row.names=F,col.names=F,quote=F)"

medianCov="$(realpath "${batch}_medianCov.transposed.bed")"