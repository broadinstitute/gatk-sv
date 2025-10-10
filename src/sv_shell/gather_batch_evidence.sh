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
counts=$(jq -r ".counts[]" "${input_json}")
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



# Condense Read Counts
# ---------------------------------------------------------------------------------------------------------------------

# Note that the WDL version implements this as a for loop over all the 'samples' in the input;
# however, since in the single-sample mode the 'samples' list contains one sample only,
# for simplicity, the in the sv-shell version we're not implementing the loop.

echo -e "${MAGENTA}Starting Condense Read Counts.${NC}"
condense_read_counts_inputs_json="$(realpath "${output_dir}/condense_read_counts_inputs.json")"
condense_read_counts_outputs_json="$(realpath "${output_dir}/condense_read_counts_outputs.json")"

jq -n \
  --arg counts "${counts[0]}" \
  --arg sample "${samples[0]}" \
  --argfile min_interval_size <(jq '.min_interval_size' "${input_json}") \
  --argfile max_interval_size <(jq '.max_interval_size' "${input_json}") \
  '{
      "sample": $sample,
      "counts": $counts,
      "min_interval_size": $min_interval_size,
      "max_interval_size": $max_interval_size
  }' > "${condense_read_counts_inputs_json}"

bash /opt/sv_shell/condense_read_counts.sh "${condense_read_counts_inputs_json}" "${condense_read_counts_outputs_json}"

echo -e "${GREEN}Successfully finished running Condense Read Counts.${NC}"



# Merge Depth
# ---------------------------------------------------------------------------------------------------------------------

echo -e "${MAGENTA}Starting Merge Depth.${NC}"
merge_depth_inputs_json="$(realpath "${output_dir}/merge_depth_inputs.json")"
merge_depth_outputs_json="$(realpath "${output_dir}/mergge_depth_outputs.json")"

# The following args are all temporarily skipped as they are gcnv related.
#  --arg genotyped_segments_vcfs
#  --arg contig_ploidy_calls
#  --arg gcnv_qs_cutoff

jq -n \
  --argfile samples <(jq '.samples' "${input_json}") \
  --argfile defragment_max_dist <(jq '.defragment_max_dist // ""' "${input_json}") \
  --argfile std_cnmops_del <(jq '.Del // ""' "${cnmops_outputs_json}") \
  --argfile std_cnmops_dup <(jq '.Dup // ""' "${cnmops_outputs_json}") \
  --argfile large_cnmops_del <(jq '.Del // ""' "${cnmops_large_outputs_json}") \
  --argfile large_cnmops_dup <(jq '.Dup // ""' "${cnmops_large_outputs_json}") \
  --arg batch "${batch}" \
  '{
      "samples": $samples,
      "defragment_max_dist": $defragment_max_dist,
      "std_cnmops_del": $std_cnmops_del,
      "std_cnmops_dup": $std_cnmops_dup,
      "large_cnmops_del": $large_cnmops_del,
      "large_cnmops_dup": $large_cnmops_dup,
      "batch": $batch
  }' > "${merge_depth_inputs_json}"





# Median Cov
# ---------------------------------------------------------------------------------------------------------------------
echo -e "${MAGENTA}Starting Median Cov.${NC}"
merged_median_cov="$(realpath "${working_dir}/merged_median_cov.bed")"

paste -d'\t' \
  $(jq -r ".ref_panel_median_cov" "${input_json}") \
  $(jq -r ".sample_median_cov" "${input_json}") \
  > "${merged_median_cov}"

## zcat "${merged_median_cov}" > "${batch}_fixed.bed"
# since the file is already unzipped, so just cp it as following to keep original untouched,
# use the above if the file was gzipped.
cp "${merged_median_cov}" "${batch}_fixed.bed"

Rscript /opt/WGD/bin/medianCoverage.R "${batch}_fixed.bed" -H "${batch}_medianCov.bed"
Rscript -e "x <- read.table(\"${batch}_medianCov.bed\",check.names=FALSE); xtransposed <- t(x[,c(1,2)]); write.table(xtransposed,file=\"${batch}_medianCov.transposed.bed\",sep=\"\\t\",row.names=F,col.names=F,quote=F)"

output_median_cov="$(realpath "${batch}_medianCov.transposed.bed")"

echo -e "${GREEN}Successfully finished running Median Cov.${NC}"



# Preprocess PE/SR
# ---------------------------------------------------------------------------------------------------------------------

echo -e "${MAGENTA}Starting Preprocess PE/SR.${NC}"
preprocess_pesr_inputs_json="$(realpath "${output_dir}/preprocess_pesr_inputs.json")"
preprocess_pesr_outputs_json="$(realpath "${output_dir}/preprocess_pesr_outputs.json")"

# TODO: the following includes scramble VCFs, but they are not included in the WDL equivalent, bug in WDL or by design?
jq -n \
  --argfile samples <(jq '.samples' "${input_json}") \
  --argfile manta_vcfs <(jq '.manta_vcfs // ""' "${input_json}") \
  --argfile scramble_vcfs <(jq '.scramble_vcfs // ""' "${input_json}") \
  --argfile wham_vcfs <(jq '.wham_vcfs // ""' "${input_json}") \
  --argfile dragen_vcfs <(jq '.dragen_vcfs // []' "${input_json}") \
  --argfile contigs <(jq '.primary_contigs_fai' "${input_json}") \
  --argfile min_svsize <(jq '.min_svsize' "${input_json}") \
  --arg batch "${batch}" \
  '{
      "samples": $samples,
      "manta_vcfs": $manta_vcfs,
      "scramble_vcfs": $scramble_vcfs,
      "wham_vcfs": $wham_vcfs,
      "dragen_vcfs": $dragen_vcfs,
      "contigs": $contigs,
      "min_svsize": $min_svsize,
      "batch": $batch
  }' > "${preprocess_pesr_inputs_json}"

bash /opt/sv_shell/preprocess_pesr.sh "${preprocess_pesr_inputs_json}" "${preprocess_pesr_outputs_json}"

echo -e "${GREEN}Successfully finished running Preprocess PE/SR.${NC}"



# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

merged_BAF_task_out=$(jq -r ".merged_BAF" "${batch_evidence_merging_outputs_json}")
merged_BAF="${output_dir}/$(basename "${merged_BAF_task_out}")"
mv "${merged_BAF_task_out}" "${merged_BAF}"

merged_BAF_index_task_out=$(jq -r ".merged_BAF_index" "${batch_evidence_merging_outputs_json}")
merged_BAF_index="${output_dir}/$(basename "${merged_BAF_index_task_out}")"
mv "${merged_BAF_index_task_out}" "${merged_BAF_index}"

merged_SR_task_out=$(jq -r ".merged_SR" "${batch_evidence_merging_outputs_json}")
merged_SR="${output_dir}/$(basename "${merged_SR_task_out}")"
mv "${merged_SR_task_out}" "${merged_SR}"

merged_SR_index_task_out=$(jq -r ".merged_SR_index" "${batch_evidence_merging_outputs_json}")
merged_SR_index="${output_dir}/$(basename "${merged_SR_index_task_out}")"
mv "${merged_SR_index_task_out}" "${merged_SR_index}"

merged_PE_task_out=$(jq -r ".merged_PE" "${batch_evidence_merging_outputs_json}")
merged_PE="${output_dir}/$(basename "${merged_PE_task_out}")"
mv "${merged_PE_task_out}" "${merged_PE}"

merged_PE_index_task_out=$(jq -r ".merged_PE_index" "${batch_evidence_merging_outputs_json}")
merged_PE_index="${output_dir}/$(basename "${merged_PE_index_task_out}")"
mv "${merged_PE_index_task_out}" "${merged_PE_index}"

merged_bincov_out="${output_dir}/$(basename "${merged_bincov_}")"
merged_bincov_idx_out="${output_dir}/$(basename "${merged_bincov_idx_}")"
mv "${merged_bincov_}" "${merged_bincov_out}"
mv "${merged_bincov_idx_}" "${merged_bincov_idx_out}"

ploidy_matrix_task_out=$(jq -r ".ploidy_matrix" "${ploidy_estimation_outputs_json}")
ploidy_matrix="${output_dir}/$(basename "${ploidy_matrix_task_out}")"
mv "${ploidy_matrix_task_out}" "${ploidy_matrix}"

ploidy_plots_task_out=$(jq -r ".ploidy_plots" "${ploidy_estimation_outputs_json}")
ploidy_plots="${output_dir}/$(basename "${ploidy_plots_task_out}")"
mv "${ploidy_plots_task_out}" "${ploidy_plots}"

combined_ped_file_output="${output_dir}/$(basename "${combined_ped_file}")"
mv "${combined_ped_file}" "${combined_ped_file_output}"

merge_depth_del_task_out=$(jq -r ".Del" "${merge_depth_outputs_json}")
merge_depth_del="${output_dir}/$(basename "${merge_depth_del_task_out}")"
mv "${merge_depth_del_task_out}" "${merge_depth_del}"

merge_depth_dup_task_out=$(jq -r ".Dup" "${merge_depth_outputs_json}")
merge_depth_dup="${output_dir}/$(basename "${merge_depth_dup_task_out}")"
mv "${merge_depth_dup_task_out}" "${merge_depth_dup}"






outputs_json=$(jq -n \
  --arg del "${del_out}" \
  --arg del_index "${del_index_out}" \
  --arg dup "${dup_out}" \
  --arg dup_index "${dup_index_out}" \
  '{Del: $del, Del_idx: $del_index, Dup: $dup, Dup_idx: $dup_index}' )
echo "${outputs_json}" > "${output_json_filename}"

echo "Successfully finished Gather Batch Evidence, output json filename: ${output_json_filename}"