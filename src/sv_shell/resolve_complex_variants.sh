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
  output_dir=$(mktemp -d ${SV_SHELL_BASE_DIR}/output_resolve_complex_variants_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir="$(mktemp -d ${SV_SHELL_BASE_DIR}/wd_resolve_complex_variants_XXXXXXXX)"
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
ResolveCpxInv_wd=$(mktemp -d "${SV_SHELL_BASE_DIR}/wd_ResolveCpxInv_XXXXXXXX")
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

ResolveCpxInv_resolved_vcf_merged=$(jq -r ".resolved_vcf_merged" "${ResolveCpxInv_outputs_json}")

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


# ResolveCpxAll
# ---------------------------------------------------------------------------------------------------------------------
ResolveCpxAll_wd=$(mktemp -d "${SV_SHELL_BASE_DIR}/wd_ResolveCpxAll_XXXXXXXX")
ResolveCpxAll_wd="$(realpath ${ResolveCpxAll_wd})"
ResolveCpxAll_inputs_json="${ResolveCpxAll_wd}/inputs.json"
ResolveCpxAll_outputs_json="${ResolveCpxAll_wd}/outputs.json"

jq -n \
  --slurpfile inputs "${input_json}" \
  --arg vcf "${BreakpointOverlap_out}" \
  --arg prefix "${cohort_name}.all" \
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
  }' > "${ResolveCpxAll_inputs_json}"

bash /opt/sv_shell/resolve_complex_sv.sh "${ResolveCpxAll_inputs_json}" "${ResolveCpxAll_outputs_json}"
cd "${working_dir}"

ResolveCpxAll_resolved_vcf_merged=$(jq -r ".resolved_vcf_merged" "${ResolveCpxAll_outputs_json}")


# IntegrateResolvedVcfs
# Integrate inv-only and all-variants resolved VCFs
# ---------------------------------------------------------------------------------------------------------------------

IntegrateResolvedVcfs_integrated_vcf="$(realpath "${cohort_name}.resolved.vcf.gz")"
mkdir tmp
python /opt/sv-pipeline/04_variant_resolution/scripts/integrate_resolved_vcfs.py \
  --all-vcf "${ResolveCpxAll_resolved_vcf_merged}" \
  --inv-only-vcf "${ResolveCpxInv_resolved_vcf_merged}" \
  | bcftools sort --temp-dir ./tmp -Oz -o "${IntegrateResolvedVcfs_integrated_vcf}"
tabix "${IntegrateResolvedVcfs_integrated_vcf}"


# RenameVariants
# Apply consistent variant naming scheme to integrated VCF
# ---------------------------------------------------------------------------------------------------------------------

RenameVariants_renamed_vcf="$(realpath "${cohort_name}.renamed.vcf.gz")"

python /opt/sv-pipeline/04_variant_resolution/scripts/rename.py \
  --prefix "${cohort_name}" \
  "${IntegrateResolvedVcfs_integrated_vcf}" - \
  | bgzip \
  > "${RenameVariants_renamed_vcf}"
tabix "${RenameVariants_renamed_vcf}"



# UpdateBothsidePass
# Update SR background fail & bothside pass files
# ---------------------------------------------------------------------------------------------------------------------

UpdateBothsidePass_updated_list="$(realpath "${cohort_name}.sr_bothside_pass.updated3.txt")"

# append new ids to original list
svtk vcf2bed "${RenameVariants_renamed_vcf}" int.bed -i MEMBERS --no-samples --no-header

# match id one per line
# if an id is not found in the vcf, use previous id (in case vcf is a shard/subset)
# also sort by first column, which is support fraction for a bothside pass list
awk -F'[,\t]' -v OFS='\t' \
  '{ \
    if (ARGIND==1) for(i=6; i<=NF; ++i) MAP[$i]=$4; \
    else if ($NF in MAP) print $0,MAP[$NF]; \
    else print $0,$NF; \
  }' int.bed "${cluster_bothside_pass_lists}" \
  | sort -k1,1n \
  > "${UpdateBothsidePass_updated_list}"


# UpdateBackgroundFail
# Update SR background fail & bothside pass files
# ---------------------------------------------------------------------------------------------------------------------

UpdateBackgroundFail_updated_list="$(realpath "${cohort_name}.sr_background_fail.updated3.txt")"

# append new ids to original list
svtk vcf2bed "${RenameVariants_renamed_vcf}" int.bed -i MEMBERS --no-samples --no-header

# match id one per line
# if an id is not found in the vcf, use previous id (in case vcf is a shard/subset)
# also sort by first column, which is support fraction for a bothside pass list
awk -F'[,\t]' -v OFS='\t' \
  '{ \
    if (ARGIND==1) for(i=6; i<=NF; ++i) MAP[$i]=$4; \
    else if ($NF in MAP) print $0,MAP[$NF]; \
    else print $0,$NF; \
  }' int.bed "${cluster_background_fail_lists}" \
  | sort -k1,1n \
  > "${UpdateBackgroundFail_updated_list}"



# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------


complex_resolve_vcfs_output_dir="${output_dir}/$(basename "${RenameVariants_renamed_vcf}")"
mv "${RenameVariants_renamed_vcf}" "${complex_resolve_vcfs_output_dir}"
mv "${RenameVariants_renamed_vcf}.tbi" "${complex_resolve_vcfs_output_dir}.tbi"

complex_resolve_bothside_pass_list_output_dir="${output_dir}/$(basename "${UpdateBothsidePass_updated_list}")"
mv "${UpdateBothsidePass_updated_list}" "${complex_resolve_bothside_pass_list_output_dir}"

complex_resolve_background_fail_list_output_dir="${output_dir}/$(basename "${UpdateBackgroundFail_updated_list}")"
mv "${UpdateBackgroundFail_updated_list}" "${complex_resolve_background_fail_list_output_dir}"

breakpoint_overlap_dropped_record_vcfs_output_dir="${output_dir}/$(basename "${BreakpointOverlap_dropped_record_vcf}")"
mv "${BreakpointOverlap_dropped_record_vcf}" "${breakpoint_overlap_dropped_record_vcfs_output_dir}"
mv "${BreakpointOverlap_dropped_record_vcf}.tbi" "${breakpoint_overlap_dropped_record_vcfs_output_dir}.tbi"


jq -n \
  --arg complex_resolve_vcfs "${complex_resolve_vcfs_output_dir}" \
  --arg complex_resolve_vcf_indexes "${complex_resolve_vcfs_output_dir}.tbi" \
  --arg complex_resolve_bothside_pass_list "${complex_resolve_bothside_pass_list_output_dir}" \
  --arg complex_resolve_background_fail_list "${complex_resolve_background_fail_list_output_dir}" \
  --arg breakpoint_overlap_dropped_record_vcfs "${breakpoint_overlap_dropped_record_vcfs_output_dir}" \
  --arg breakpoint_overlap_dropped_record_vcf_indexes "${breakpoint_overlap_dropped_record_vcfs_output_dir}.tbi" \
  '{
      "complex_resolve_vcfs": $complex_resolve_vcfs,
      "complex_resolve_vcf_indexes": $complex_resolve_vcf_indexes,
      "complex_resolve_bothside_pass_list": $complex_resolve_bothside_pass_list,
      "complex_resolve_background_fail_list": $complex_resolve_background_fail_list,
      "breakpoint_overlap_dropped_record_vcfs": $breakpoint_overlap_dropped_record_vcfs,
      "breakpoint_overlap_dropped_record_vcf_indexes": $breakpoint_overlap_dropped_record_vcf_indexes
  }' > "${output_json_filename}"

echo "Successfully finished resolve complex variants, output json filename: ${output_json_filename}"
