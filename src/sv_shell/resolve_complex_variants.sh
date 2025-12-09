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
  output_dir=$(mktemp -d /output_resolve_complex_variants_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir="$(mktemp -d /wd_resolve_complex_variants_XXXXXXXX)"
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"

cohort_name=$(jq -r ".cohort_name" "$input_json")
cluster_vcfs=$(jq -r ".cluster_vcfs" "$input_json")
cluster_bothside_pass_lists=$(jq -r ".cluster_bothside_pass_lists" "$input_json")
cluster_background_fail_lists=$(jq -r ".cluster_background_fail_lists" "$input_json")


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

# SubsetInversions
# Subset inversions from PESR+RD VCF
# ---------------------------------------------------------------------------------------------------------------------

SubsetInversions_filtered_vcf="${cohort_name}.inversions_only.vcf.gz"
bcftools view --no-version --no-update -i 'INFO/SVTYPE="INV"' -O z -o "${SubsetInversions_filtered_vcf}" "${cluster_vcfs}"
tabix "${SubsetInversions_filtered_vcf}"

SubsetInversions_filtered_vcf="$(realpath "${SubsetInversions_filtered_vcf}")"
SubsetInversions_filtered_vcf_idx="${SubsetInversions_filtered_vcf}.tbi"


# ResolveCpxInv
# ---------------------------------------------------------------------------------------------------------------------
ResolveCpxInv_wd=$(mktemp -d "/wd_ResolveCpxInv_XXXXXXXX")
ResolveCpxInv_wd="$(realpath ${ResolveCpxInv_wd})"
ResolveCpxInv_inputs_json="${ResolveCpxInv_wd}/inputs.json"
ResolveCpxInv_outputs_json="${ResolveCpxInv_wd}/outputs.json"

jq -n \
  --slurpfile inputs "${input_json}" \
  --arg vcf "${SubsetInversions_filtered_vcf}" \
  --arg prefix "${cohort_name}.inv_only" \
  '{
    "vcf": $vcf,
    "prefix": $prefix,
    "max_shard_size": $inputs[0].max_shard_size,
    "cytobands": $inputs[0].cytobands,
    "disc_files": $inputs[0].disc_files,
    "mei_bed": $inputs[0].mei_bed,
    "pe_exclude_list": $inputs[0].pe_exclude_list,
    "rf_cutoff_files": $inputs[0].rf_cutoff_files,
    "ref_dict": $inputs[0].ref_dict,
    "precluster_distance": 2000,
    "precluster_overlap_frac": 0.000000001
  }' > "${ResolveCpxInv_inputs_json}"

bash /opt/sv_shell/resolve_complex_sv.sh "${ResolveCpxInv_inputs_json}" "${ResolveCpxInv_outputs_json}"
cd "${working_dir}"

# BreakpointOverlap
# ---------------------------------------------------------------------------------------------------------------------
BreakpointOverlap_prefix="${cohort_name}.breakpoint_overlap"
python /opt/sv-pipeline/04_variant_resolution/scripts/overlap_breakpoint_filter.py \
  "${cluster_vcfs}" \
  "${cluster_bothside_pass_lists}" \
  "${cluster_background_fail_lists}" \
  "${BreakpointOverlap_prefix}.dropped_records.vcf.gz" \
  | bgzip \
  > "${BreakpointOverlap_prefix}.vcf.gz"
tabix "${BreakpointOverlap_prefix}.vcf.gz"
tabix "${BreakpointOverlap_prefix}.dropped_records.vcf.gz"

BreakpointOverlap_out="$(realpath "${BreakpointOverlap_prefix}.vcf.gz")"
BreakpointOverlap_out_index="$(realpath "${BreakpointOverlap_prefix}.vcf.gz.tbi")"
BreakpointOverlap_dropped_record_vcf="$(realpath "${BreakpointOverlap_prefix}.dropped_records.vcf.gz")"
BreakpointOverlap_dropped_record_vcf_index="$(realpath "${BreakpointOverlap_prefix}.dropped_records.vcf.gz.tbi")"
