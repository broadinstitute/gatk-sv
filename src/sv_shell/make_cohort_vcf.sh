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


# ResolveComplexVariants
# ---------------------------------------------------------------------------------------------------------------------
ResolveComplexVariants_output_dir=$(realpath $(mktemp -d "/output_ResolveComplexVariants_XXXXXXXX"))
ResolveComplexVariants_inputs_json_filename="${ResolveComplexVariants_output_dir}/inputs.json"
ResolveComplexVariants_outputs_json_filename="${ResolveComplexVariants_output_dir}/outputs.json"

jq -n \
  --slurpfile inputs "${input_json}" \
  --slurpfile cb "${CombineBatches_outputs_json_filename}" \
  '{
    "cohort_name": $inputs[0].cohort_name,
    "cluster_vcfs": $cb[0].combined_vcfs,
    "cluster_bothside_pass_lists": $cb[0].cluster_bothside_pass_lists,
    "cluster_background_fail_lists": $cb[0].cluster_background_fail_lists,
    "disc_files": $inputs[0].disc_files,
    "rf_cutoff_files": $inputs[0].rf_cutoff_files,
    "contig_list": $inputs[0].contig_list,
    "cytobands": $inputs[0].cytobands,
    "mei_bed": $inputs[0].mei_bed,
    "pe_exclude_list": $inputs[0].pe_exclude_list,
    "ref_dict": $inputs[0].reference_dict,
    "max_shard_size": $inputs[0].max_shard_size_resolve
  }' > "${ResolveComplexVariants_inputs_json_filename}"

bash /opt/sv_shell/resolve_complex_variants.sh \
  "${ResolveComplexVariants_inputs_json_filename}" \
  "${ResolveComplexVariants_outputs_json_filename}" \
  "${ResolveComplexVariants_output_dir}"


# GenotypeComplexVariants
# ---------------------------------------------------------------------------------------------------------------------
GenotypeComplexVariants_output_dir=$(realpath $(mktemp -d "/output_GenotypeComplexVariants_XXXXXXXX"))
GenotypeComplexVariants_inputs_json_filename="${GenotypeComplexVariants_output_dir}/inputs.json"
GenotypeComplexVariants_outputs_json_filename="${GenotypeComplexVariants_output_dir}/outputs.json"

jq -n \
  --slurpfile inputs "${input_json}" \
  --slurpfile resolve_cpx "${ResolveComplexVariants_outputs_json_filename}" \
  '{
    "cohort_name": $inputs[0].cohort_name,
    "batches": $inputs[0].batches,
    "merge_vcfs": $inputs[0].merge_complex_genotype_vcfs,
    "complex_resolve_vcfs": $resolve_cpx[0].complex_resolve_vcfs,
    "complex_resolve_vcf_indexes": $resolve_cpx[0].complex_resolve_vcf_indexes,
    "depth_vcfs": $inputs[0].depth_vcfs,
    "ped_file": $inputs[0].ped_file,
    "bincov_files": $inputs[0].bincov_files,
    "depth_gt_rd_sep_files": $inputs[0].depth_gt_rd_sep_files,
    "median_coverage_files": $inputs[0].median_coverage_files,
    "bin_exclude": $inputs[0].bin_exclude,
    "contig_list": $inputs[0].contig_list,
    "ref_dict": $inputs[0].reference_dict
  }' > "${GenotypeComplexVariants_inputs_json_filename}"

bash /opt/sv_shell/genotype_complex_variants.sh \
  "${GenotypeComplexVariants_inputs_json_filename}" \
  "${GenotypeComplexVariants_outputs_json_filename}" \
  "${GenotypeComplexVariants_output_dir}"


# CleanVcf
# ---------------------------------------------------------------------------------------------------------------------
CleanVcf_output_dir=$(realpath $(mktemp -d "/output_CleanVcf_XXXXXXXX"))
CleanVcf_inputs_json_filename="${CleanVcf_output_dir}/inputs.json"
CleanVcf_outputs_json_filename="${CleanVcf_output_dir}/outputs.json"

jq -n \
  --slurpfile inputs "${input_json}" \
  --slurpfile rcpx "${ResolveComplexVariants_outputs_json_filename}" \
  --slurpfile gcpx "${GenotypeComplexVariants_outputs_json_filename}" \
  '{
    "cohort_name": $inputs[0].cohort_name,
    "complex_genotype_vcf": $gcpx[0].complex_genotype_merged_vcf,
    "complex_resolve_bothside_pass_list": $rcpx[0].complex_resolve_bothside_pass_list,
    "complex_resolve_background_fail_list": $rcpx[0].complex_resolve_background_fail_list,
    "ped_file": $inputs[0].ped_file,
    "contig_list": $inputs[0].contig_list,
    "allosome_fai": $inputs[0].allosome_fai,
    "chr_x": $inputs[0].chr_x,
    "chr_y": $inputs[0].chr_y,
    "HERVK_reference": $inputs[0].HERVK_reference,
    "LINE1_reference": $inputs[0].LINE1_reference,
    "intron_reference": $inputs[0].intron_reference,
    "outlier_samples_list": $inputs[0].outlier_samples_list,
  }' > "${CleanVcf_inputs_json_filename}"

bash /opt/sv_shell/clean_vcf.sh \
  "${CleanVcf_inputs_json_filename}" \
  "${CleanVcf_outputs_json_filename}" \
  "${CleanVcf_output_dir}"


# -------------------------------------------------------
# ======================= Output ========================
# -------------------------------------------------------

jq -n \
  --slurpfile clean_vcf "${CleanVcf_outputs_json_filename}" \
  --slurpfile cb "${CombineBatches_outputs_json_filename}" \
  --slurpfile rcpx "${ResolveComplexVariants_outputs_json_filename}" \
  --slurpfile gcpx "${GenotypeComplexVariants_outputs_json_filename}" \
  '{
     "vcf": $clean_vcf[0].cleaned_vcf,
     "vcf_index": $clean_vcf[0].cleaned_vcf_index,
     "cluster_vcf": $cb[0].combine_batches_merged_vcf,
     "cluster_vcf_index": $cb[0].combine_batches_merged_vcf_index,
     "complex_resolve_vcf": $rcpx[0].complex_resolve_vcfs,
     "complex_resolve_vcf_index": $rcpx[0].complex_resolve_vcf_indexes,
     "complex_genotype_vcf": $gcpx[0].complex_genotype_merged_vcf,
     "complex_genotype_vcf_index": $gcpx[0].complex_genotype_merged_vcf_index,
     "combined_vcfs": $cb[0].combined_vcfs,
     "combined_vcf_indexes": $cb[0].combined_vcf_indexes,
     "cluster_bothside_pass_lists": $cb[0].cluster_bothside_pass_lists,
     "cluster_background_fail_lists": $cb[0].cluster_background_fail_lists,
     "complex_resolve_bothside_pass_list": $rcpx[0].complex_resolve_bothside_pass_list,
     "complex_resolve_background_fail_list": $rcpx[0].complex_resolve_background_fail_list,
     "breakpoint_overlap_dropped_record_vcfs": $rcpx[0].breakpoint_overlap_dropped_record_vcfs,
     "breakpoint_overlap_dropped_record_vcf_indexes": $rcpx[0].breakpoint_overlap_dropped_record_vcf_indexes,
     "complex_genotype_vcfs": $gcpx[0].complex_genotype_merged_vcf,
     "complex_genotype_vcf_indexes": $gcpx[0].complex_genotype_merged_vcf_index
   }' > "${output_json_filename}"

echo "Finished make cohort VCF, output json filename: ${output_json_filename}"