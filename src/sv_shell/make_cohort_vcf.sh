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
  output_dir=$(mktemp -d /output_make_cohort_vcf_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_make_cohort_vcf_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
echo "Make cohort vcf Working directory: ${working_dir}"

bincov_matrix=($(jq -r '.bincov_matrix' "$input_json"))
cohort_id=($(jq -r '.cohort_id' "$input_json"))


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

# CombineBatches
# ---------------------------------------------------------------------------------------------------------------------
CombineBatches_output_dir=$(realpath $(mktemp -d "/output_CombineBatches_XXXXXXXX"))
CombineBatches_inputs_json_filename="${CombineBatches_output_dir}/inputs.json"
CombineBatches_outputs_json_filename="${CombineBatches_output_dir}/outputs.json"

jq -n \
  --slurpfile inputs "${input_json}" \
  '{
      "cohort_name": $inputs[0].cohort_name,
      "batches": $inputs[0].batches,
      "ped_file": $inputs[0].ped_file,
      "pesr_vcfs": $inputs[0].pesr_vcfs,
      "depth_vcfs": $inputs[0].depth_vcfs,
      "raw_sr_bothside_pass_files": $inputs[0].raw_sr_bothside_pass_files,
      "raw_sr_background_fail_files": $inputs[0].raw_sr_background_fail_files,
      "contig_list": $inputs[0].contig_list,
      "min_sr_background_fail_batches": $inputs[0].min_sr_background_fail_batches,
      "clustering_config_part1": $inputs[0].clustering_config_part1,
      "stratification_config_part1": $inputs[0].stratification_config_part1,
      "clustering_config_part2": $inputs[0].clustering_config_part2,
      "stratification_config_part2": $inputs[0].stratification_config_part2,
      "track_names": $inputs[0].track_names,
      "track_bed_files": $inputs[0].track_bed_files,
      "reference_fasta": $inputs[0].reference_fasta,
      "reference_fasta_fai": $inputs[0].reference_fasta_fai,
      "reference_dict": $inputs[0].reference_dict,
      "chr_x": $inputs[0].chr_x,
      "chr_y": $inputs[0].chr_y,
  }' > "${CombineBatches_inputs_json_filename}"

bash /opt/sv_shell/combine_batches.sh \
  "${CombineBatches_inputs_json_filename}" \
  "${CombineBatches_outputs_json_filename}" \
  "${CombineBatches_output_dir}"