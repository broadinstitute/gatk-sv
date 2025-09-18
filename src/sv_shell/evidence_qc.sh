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
run_ploidy=$(jq -r ".run_ploidy" "${input_json}")


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
merged_bincov=$(jq -r ".merged_bincov" "${make_bincov_matrix_output_json}")
merged_bincov_idx=$(jq -r ".merged_bincov_idx" "${make_bincov_matrix_output_json}")


# --- MedianCov
zcat "${merged_bincov}" > "${batch}_fixed.bed"
Rscript /opt/WGD/bin/medianCoverage.R "${batch}_fixed.bed" -H "${batch}_medianCov.bed"
Rscript -e "x <- read.table(\"${batch}_medianCov.bed\",check.names=FALSE); xtransposed <- t(x[,c(1,2)]); write.table(xtransposed,file=\"${batch}_medianCov.transposed.bed\",sep=\"\\t\",row.names=F,col.names=F,quote=F)"

medianCov="$(realpath "${batch}_medianCov.transposed.bed")"

# ---- WGD
wgd_input_json="$(realpath "${output_dir}/wgd_input.json")"
jq '
  {
    batch,
    bincov_matrix,
    wgd_scoring_mask
  }
' "${input_json}" |
  jq '.' > "${wgd_input_json}"
wgd_output_json="$(realpath "${output_dir}/wgd_output.json")"
bash /wgd.sh "${wgd_input_json}" "${wgd_output_json}"

WGD_scores=$(jq -r ".WGD_scores" "${wgd_output_json}")
WGD_dist=$(jq -r ".WGD_dist" "${wgd_output_json}")
WGD_matrix=$(jq -r ".WGD_matrix" "${wgd_output_json}")



# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

outputs_json=$(jq -n \
  --arg WGD_dist "${WGD_dist}" \
  --arg WGD_matrix "${WGD_matrix}" \
  --arg WGD_scores "${WGD_scores}" \
  --arg bincov_matrix "${merged_bincov}" \
  --arg bincov_matrix_index "${merged_bincov_idx}" \
  --arg bincov_median "${medianCov}" \
  '{
     "WGD_dist": $WGD_dist,
     "WGD_matrix": $WGD_matrix,
     "WGD_scores": $WGD_scores,
     "bincov_matrix": $bincov_matrix,
     "bincov_matrix_index": $bincov_matrix_index,
     "bincov_median": $bincov_median
   }' \
)
echo "${outputs_json}" > "${output_json_filename}"
