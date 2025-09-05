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
  output_dir=$(mktemp -d output_gather_batch_evidence_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d wd_gather_batch_evidence_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

batch=$(jq -r ".batch" "${input_json}")
samples=$(jq -r ".samples[]" "${input_json}")
ped_file=$(jq -r ".ped_file" "${input_json}")
ref_panel_samples=($(jq -r '.ref_panel_samples[]' "$input_json"))
append_first_sample_to_ped=($(jq -r '.append_first_sample_to_ped' "$input_json"))

sample_bincov_matrix=$(jq -r ".sample_bincov_matrix" "${input_json}")
ref_panel_bincov_matrix=$(jq -r ".ref_panel_bincov_matrix" "${input_json}")
reference_dict=$(jq -r ".reference_dict" "${input_json}")

# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


# ---- make binned coverage matrix
counts_files=("${ref_panel_bincov_matrix}" "${sample_bincov_matrix}")
make_bin_cov_matrix_json="make_bincov_matrix.json"
jq -n \
  --argfile s <(jq '.samples' "${input_json}") \
  --arg c "${counts_files[*]}" \
  --argfile r <(jq '.ref_panel_samples' "${input_json}") \
  --arg b "${ref_panel_bincov_matrix}" \
  --arg p "${reference_dict}" \
  --arg t "${batch}" \
  '{samples: $s, count_files: ($c | split(" ")), bincov_matrix_samples: $r, bincov_matrix: $b, reference_dict: $p, batch: $t, skip_bin_size_filter: true}' > "${make_bin_cov_matrix_json}"

bash /make_bincov_matrix.sh "${make_bin_cov_matrix_json}"



# TODO: you will get this file from running ploidy and then from the following direcotry in its output:
# ploidy_est/sample_sex_assignments.txt.gz
sample_sex_assignments="/sample_sex_assignments.txt.gz"


if [[ -n "${ref_panel_samples}" ]]; then
  samples_batch=("${ref_panel_samples[@]}")
else
  samples_batch=("${samples[@]}")
fi
samples_batch_file="samples_batch_file.txt"
printf "%s\n" "${samples_batch[@]}" > "${samples_batch_file}"

# validate PED file
python /opt/sv-pipeline/scripts/validate_ped.py -p "${ped_file}" -s "${samples_batch_file}"

# subset PED file
ped_subset_filename="$(basename "${ped_file}" .ped).${batch}.ped"
awk 'FNR==NR {a[$1]; next}; $2 in a' "${samples_batch_file}" "${ped_file}" > "${ped_subset_filename}"

# Add case sample to PED
combined_ped_file="combined_ped_file.ped"
if [[ "${append_first_sample_to_ped}" == "true" ]]; then

  sample_id=${samples[0]}

  RECORD=$(gunzip -c "${sample_sex_assignments}" | { grep -w "^${sample_id}" || true; })
  if [ -z "$RECORD" ]; then
    >&2 echo "Error: Sample ${sample_id} not found in ploidy calls"
    exit 1
  fi
  SEX=$(echo "$RECORD" | cut -f2)
  # TODO: would not it be better to just accept the fact that ped already has the sample and move on? or we should really error-out?
  awk -v sample="${sample_id}" '$2 == sample { print "ERROR: A sample with the name "sample" is already present in the ped file." > "/dev/stderr"; exit 1; }' < "${ped_subset_filename}"
  awk -v sample="${sample_id}" -v sex="${SEX}" '{print} END {OFS="\t"; print "case_sample",sample,"0","0",sex,"1" }' < "${ped_subset_filename}" > "${combined_ped_file}"
fi







































median_cov_filename="/all_samples_medianCov.transposed.bed"


paste -d'\t' all_samples_medianCov.transposed.bed wd_evidence_qc_XGqRWD02/NA12878_medianCov.transposed.bed > merged_mediancov.bed


zcat "${merged_bincov}" > "${batch}_fixed.bed"
Rscript /opt/WGD/bin/medianCoverage.R "${batch}_fixed.bed" -H "${batch}_medianCov.bed"
Rscript -e "x <- read.table(\"${batch}_medianCov.bed\",check.names=FALSE); xtransposed <- t(x[,c(1,2)]); write.table(xtransposed,file=\"${batch}_medianCov.transposed.bed\",sep=\"\\t\",row.names=F,col.names=F,quote=F)"
