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
  output_dir=$(mktemp -d /output_median_cov_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_median_cov_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
echo "Median coverage Working directory: ${working_dir}"

bincov_matrix=($(jq -r '.bincov_matrix' "$input_json"))
cohort_id=($(jq -r '.cohort_id' "$input_json"))


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


zcat "${bincov_matrix}" > "${cohort_id}_fixed.bed"
Rscript /opt/WGD/bin/medianCoverage.R "${cohort_id}_fixed.bed" -H "${cohort_id}_medianCov.bed"
Rscript -e "x <- read.table(\"${cohort_id}_medianCov.bed\",check.names=FALSE); xtransposed <- t(x[,c(1,2)]); write.table(xtransposed,file=\"${cohort_id}_medianCov.transposed.bed\",sep=\"\\t\",row.names=F,col.names=F,quote=F)"


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

median_cov_out="${output_dir}/$(basename "${cohort_id}_medianCov.transposed.bed")"
mv "${cohort_id}_medianCov.transposed.bed" "${median_cov_out}"

outputs_json=$(jq -n \
  --arg medcov "${median_cov_out}" \
  '{medianCov: $medcov}' )
echo "${outputs_json}" > "${output_json_filename}"

echo "Finished MedianCoverage successfully, output json filename: ${output_json_filename}"
