#!/bin/bash

set -Exeuo pipefail

RED='\033[0;31m'
BOLD_RED="\033[1;31m"
GREEN='\033[0;32m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# -------------------------------------------------------
# ==================== Input & Setup ====================
# -------------------------------------------------------

input_json=${1}
output_json_filename=${2:-""}
output_dir=${3:-""}

input_json="$(realpath ${input_json})"

if [ -z "${output_dir}" ]; then
  output_dir=$(mktemp -d /output_gather_batch_evidence_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_gather_batch_evidence_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

batch=$(jq -r ".batch" "${input_json}")
samples=$(jq -r ".samples[]" "${input_json}")
ref_panel_samples=($(jq -r '.ref_panel_samples[]' "$input_json"))
all_samples=($(jq -r '(.samples + .ref_panel_samples)[]' "$input_json"))
ped_file=$(jq -r ".ped_file" "${input_json}")
append_first_sample_to_ped=($(jq -r '.append_first_sample_to_ped' "$input_json"))
sample_bincov_matrix=$(jq -r ".sample_bincov_matrix" "${input_json}")
ref_panel_bincov_matrix=$(jq -r ".ref_panel_bincov_matrix" "${input_json}")
reference_dict=$(jq -r ".reference_dict" "${input_json}")
sd_locs_vcf=$(jq -r ".sd_locs_vcf" "${input_json}")
primary_contigs_fai=$(jq -r ".primary_contigs_fai" "${input_json}")
subset_primary_contigs=$(jq -r ".subset_primary_contigs" "${input_json}")
rename_samples=$(jq -r ".rename_samples" "${input_json}")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------


# make binned coverage matrix
# ---------------------------------------------------------------------------------------------------------------------
echo -e "${MAGENTA}Running make binned coverage matrix.${NC}"
counts_files=("${sample_bincov_matrix}")
make_bin_cov_matrix_inputs_json="$(realpath "${output_dir}/make_bincov_matrix_inputs.json")"
make_bin_cov_matrix_outputs_json="$(realpath "${output_dir}/make_bincov_matrix_outputs.json")"
jq -n \
  --argfile s <(jq '.ref_panel_samples + .samples' "${input_json}") \
  --arg c "${counts_files[*]}" \
  --argfile r <(jq '.ref_panel_samples' "${input_json}") \
  --arg b "${ref_panel_bincov_matrix}" \
  --arg p "${reference_dict}" \
  --arg t "${batch}" \
  '{
      "samples": $s,
      "count_files": ($c | split(" ")),
      "bincov_matrix_samples": $r,
      "bincov_matrix": $b,
      "reference_dict": $p,
      "batch": $t,
      "skip_bin_size_filter": true
  }' > "${make_bin_cov_matrix_inputs_json}"

bash /opt/sv_shell/make_bincov_matrix.sh "${make_bin_cov_matrix_inputs_json}" "${make_bin_cov_matrix_outputs_json}"

echo -e "${GREEN}Successfully finished make binned coverage matrix.${NC}"


# ploidy estimation
# ---------------------------------------------------------------------------------------------------------------------
echo -e "${MAGENTA}Running ploidy estimation.${NC}"
ploidy_estimation_inputs_json="$(realpath "${output_dir}/ploidy_estimation_inputs.json")"
ploidy_estimation_outputs_json="$(realpath "${output_dir}/ploidy_estimation_outputs.json")"
jq -n \
  --arg b "${batch}" \
  --argfile m <(jq '.merged_bincov' "${make_bin_cov_matrix_outputs_json}") \
  '{batch: $b, bincov_matrix: $m}' > "${ploidy_estimation_inputs_json}"

bash /opt/sv_shell/ploidy_estimation.sh "${ploidy_estimation_inputs_json}" "${ploidy_estimation_outputs_json}"

sample_sex_assignments=$(jq -r ".sample_sex_assignments" "${ploidy_estimation_outputs_json}")

if [[ -n "${ref_panel_samples}" ]]; then
  samples_batch=("${ref_panel_samples[@]}")
else
  samples_batch=("${samples[@]}")
fi
samples_batch_file="samples_batch_file.txt"
printf "%s\n" "${samples_batch[@]}" > "${samples_batch_file}"


# validate PED file
# ---------------------------------------------------------------------------------------------------------------------
echo -e "${MAGENTA}Validating, subsetting, and adding sample to PED file.${NC}"
python /opt/sv-pipeline/scripts/validate_ped.py -p "${ped_file}" -s "${samples_batch_file}"


# subset PED file
ped_subset_filename="$(basename "${ped_file}" .ped).${batch}.ped"
awk 'FNR==NR {a[$1]; next}; $2 in a' "${samples_batch_file}" "${ped_file}" > "${ped_subset_filename}"


# Add case sample to PED
combined_ped_file="$(realpath "combined_ped_file.ped")"
if [[ "${append_first_sample_to_ped}" == "true" ]]; then

  sample_id=${samples[0]}

  RECORD=$(gunzip -c "${sample_sex_assignments}" | { grep -w "^${sample_id}" || true; })
  if [ -z "$RECORD" ]; then
    >&2 echo "Error: Sample ${sample_id} not found in ploidy calls"
    exit 1
  fi
  SEX=$(echo "$RECORD" | cut -f2)

  awk -v sample="${sample_id}" '$2 == sample { print "ERROR: A sample with the name "sample" is already present in the ped file." > "/dev/stderr"; exit 1; }' < "${ped_subset_filename}"
  awk -v sample="${sample_id}" -v sex="${SEX}" '{print} END {OFS="\t"; print "case_sample",sample,"0","0",sex,"1" }' < "${ped_subset_filename}" > "${combined_ped_file}"
fi

echo -e "${GREEN}Successfully finished Validating, subsetting, and adding sample to PED file.${NC}"


# Batch evidence merging
# ---------------------------------------------------------------------------------------------------------------------
echo -e "${MAGENTA}Starting batch evidence merging.${NC}"
batch_evidence_merging_inputs_json="$(realpath "${output_dir}/batch_evidence_merging_inputs.json")"
batch_evidence_merging_outputs_json="$(realpath "${output_dir}/batch_evidence_merging_outputs.json")"
jq -n \
  --arg batch "${batch}" \
  --argfile samples <(jq '.samples + .ref_panel_samples' "${input_json}") \
  --argfile pe_files <(jq '.PE_files + (.ref_panel_PE_files // [])' "${input_json}") \
  --argfile sr_files <(jq '.SR_files + (.ref_panel_SR_files // [])' "${input_json}") \
  --argfile sd_files <(jq '.SD_files + (.ref_panel_SD_files // [])' "${input_json}") \
  --arg sd_locs_vcf "${sd_locs_vcf}" \
  --arg reference_dict "${reference_dict}" \
  --arg primary_contigs_fai "${primary_contigs_fai}" \
  --arg subset_primary_contigs "${subset_primary_contigs}" \
  --arg rename_samples "${rename_samples}" \
  '{
      "batch": $batch,
      "samples": $samples,
      "PE_files": $pe_files,
      "SR_files": $sr_files,
      "SD_files": $sd_files,
      "sd_locs_vcf": $sd_locs_vcf,
      "reference_dict": $reference_dict,
      "primary_contigs_fai": $primary_contigs_fai,
      "min_het_probability": 0.5,
      "subset_primary_contigs": $subset_primary_contigs,
      "rename_samples": $rename_samples
  }' > "${batch_evidence_merging_inputs_json}"

bash /opt/sv_shell/batch_evidence_merging.sh \
  "${batch_evidence_merging_inputs_json}" \
  "${batch_evidence_merging_outputs_json}"



# CNMOPS
# ---------------------------------------------------------------------------------------------------------------------
echo -e "${MAGENTA}Starting cnMOPS.${NC}"
cnmops_inputs_json="$(realpath "${output_dir}/cnmops_inputs.json")"
cnmops_outputs_json="$(realpath "${output_dir}/cnmops_outputs.json")"

jq -n \
  --arg r1 "3" \
  --arg r2 "10" \
  --arg batch "${batch}" \
  --argfile samples <(jq '.samples + .ref_panel_samples' "${input_json}") \
  --arg bincov_matrix "${merged_bincov_}" \
  --arg bincov_matrix_index "${merged_bincov_idx_}" \
  --argfile chrom_file <(jq '.cnmops_chrom_file' "${input_json}") \
  --argfile exclude_list <(jq '.cnmops_exclude_list' "${input_json}") \
  --argfile allo_file <(jq '.cnmops_allo_file' "${input_json}") \
  --argfile min_size <(jq '.cnmops_min_size // 1000000' "${input_json}") \
  --arg ped_file "${combined_ped_file}" \
  --arg ref_dict "${reference_dict}" \
  --arg prefix "header" \
  --arg stitch_and_clean_large_events false \
  '{
      "r1": $r1,
      "r2": $r2,
      "batch": $batch,
      "samples": $samples,
      "bincov_matrix": $bincov_matrix,
      "bincov_matrix_index": $bincov_matrix_index,
      "chrom_file": $chrom_file,
      "ped_file": $ped_file,
      "exclude_list": $exclude_list,
      "allo_file": $allo_file,
      "ref_dict": $ref_dict,
      "prefix": $prefix,
      "stitch_and_clean_large_events": $stitch_and_clean_large_events
  }' > "${cnmops_inputs_json}"

bash /opt/sv_shell/cnmops.sh "${cnmops_inputs_json}" "${cnmops_outputs_json}"

echo -e "${GREEN}Successfully finished running cnMOPS.${NC}"


# CNMOPS Large
# ---------------------------------------------------------------------------------------------------------------------
echo -e "${MAGENTA}Starting cnMOPS Large.${NC}"
cnmops_large_inputs_json="$(realpath "${output_dir}/cnmops_large_inputs.json")"
cnmops_large_outputs_json="$(realpath "${output_dir}/cnmops_large_outputs.json")"

jq -n \
  --arg r1 "1000" \
  --arg r2 "100" \
  --arg batch "${batch}" \
  --argfile samples <(jq '.samples + .ref_panel_samples' "${input_json}") \
  --arg bincov_matrix "${merged_bincov_}" \
  --arg bincov_matrix_index "${merged_bincov_idx_}" \
  --argfile chrom_file <(jq '.cnmops_chrom_file' "${input_json}") \
  --argfile exclude_list <(jq '.cnmops_exclude_list' "${input_json}") \
  --argfile allo_file <(jq '.cnmops_allo_file' "${input_json}") \
  --argfile min_size <(jq '.cnmops_large_min_size // 1000000' "${input_json}") \
  --arg ped_file "${combined_ped_file}" \
  --arg ref_dict "${reference_dict}" \
  --arg prefix "large" \
  --arg stitch_and_clean_large_events false \
  '{
      "r1": $r1,
      "r2": $r2,
      "batch": $batch,
      "samples": $samples,
      "bincov_matrix": $bincov_matrix,
      "bincov_matrix_index": $bincov_matrix_index,
      "chrom_file": $chrom_file,
      "ped_file": $ped_file,
      "exclude_list": $exclude_list,
      "allo_file": $allo_file,
      "ref_dict": $ref_dict,
      "prefix": $prefix,
      "stitch_and_clean_large_events": $stitch_and_clean_large_events
  }' > "${cnmops_large_inputs_json}"

bash /opt/sv_shell/cnmops.sh "${cnmops_large_inputs_json}" "${cnmops_large_outputs_json}"

echo -e "${GREEN}Successfully finished running cnMOPS Large.${NC}"
































median_cov_filename="/all_samples_medianCov.transposed.bed"


paste -d'\t' all_samples_medianCov.transposed.bed wd_evidence_qc_XGqRWD02/NA12878_medianCov.transposed.bed > merged_mediancov.bed


zcat "${merged_bincov}" > "${batch}_fixed.bed"
Rscript /opt/WGD/bin/medianCoverage.R "${batch}_fixed.bed" -H "${batch}_medianCov.bed"
Rscript -e "x <- read.table(\"${batch}_medianCov.bed\",check.names=FALSE); xtransposed <- t(x[,c(1,2)]); write.table(xtransposed,file=\"${batch}_medianCov.transposed.bed\",sep=\"\\t\",row.names=F,col.names=F,quote=F)"
