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
  output_dir=$(mktemp -d output_wgd_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d wd_wgd_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

batch=$(jq -r ".batch" "${input_json}")
wgd_scoring_mask=$(jq -r ".wgd_scoring_mask" "${input_json}")
bincov_matrix=$(jq -r ".bincov_matrix" "${input_json}")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

# --- build matrix

# Note that in the WDL version, 'set -eu' so zcat exiting
# after head does not result in pipe fail. Since we want to
# use 'set -Eeuo pipefail' consistently throughout all the
# scripts, we're using '|| true' in the following to avoid
# pipe fail. There are a few other alternatives, but this
# is the easiest.
zcat "${bincov_matrix}" | head -n 1 > header.txt || true
sed -i 's/#//g' header.txt

zcat "${bincov_matrix}" \
| bedtools intersect -f 0.49 -wa -u \
  -a - \
  -b "${wgd_scoring_mask}" \
| sort -Vk1,1 -k2,2n -k3,3n \
> "${batch}_WGD_scoring_matrix.bed"

wgd_matrix="${batch}_WGD_scoring_matrix_output.bed.gz"
cat header.txt "${batch}_WGD_scoring_matrix.bed" \
| bgzip -c \
> "${wgd_matrix}"


# --- WGD Score
Rscript /opt/WGD/bin/scoreDosageBiases.R -z -O . "${wgd_matrix}" "${wgd_scoring_mask}"

wgd_scores="${batch}_WGD_scores.txt.gz"
wgd_dist="${batch}_WGD_score_distributions.pdf"
mv WGD_scores.txt.gz "${wgd_scores}"
mv WGD_score_distributions.pdf "${wgd_dist}"


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

wgd_matrix_output="${output_dir}/$(basename "${wgd_matrix}")"
mv "${wgd_matrix}" "${wgd_matrix_output}"

wgd_scores_output="${output_dir}/$(basename "${wgd_scores}")"
mv "${wgd_scores}" "${wgd_scores_output}"

wgd_dist_output="${output_dir}/$(basename "${wgd_dist}")"
mv "${wgd_dist}" "${wgd_dist_output}"

outputs_json=$(jq -n \
  --arg wgd_matrix "${wgd_matrix_output}" \
  --arg wgd_scores "${wgd_scores_output}" \
  --arg wgd_dist "${wgd_dist_output}" \
  '{WGD_dist: $wgd_dist, WGD_matrix: $wgd_matrix, WGD_scores:$wgd_scores}' )
echo "${outputs_json}" > "${output_json_filename}"

echo "finished successfully"
