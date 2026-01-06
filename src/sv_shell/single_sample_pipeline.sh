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
  output_dir=$(mktemp -d /output_single_sample_XXXXXXXX)
else
  mkdir -p "${output_dir}"
fi
output_dir="$(realpath ${output_dir})"

if [ -z "${output_json_filename}" ]; then
  output_json_filename="${output_dir}/output.json"
else
  output_json_filename="$(realpath ${output_json_filename})"
fi

working_dir=$(mktemp -d /wd_single_sample_XXXXXXXX)
working_dir="$(realpath ${working_dir})"
cd "${working_dir}"
echo "Single-Sample Working directory: ${working_dir}"

batch=$(jq -r ".batch" "$input_json")
sample_id=$(jq -r ".sample_id" "$input_json")
run_vcf_qc=$(jq -r ".run_vcf_qc" "$input_json")
genome_file=$(jq -r ".genome_file" "$input_json")
wgd_scoring_mask=$(jq -r ".wgd_scoring_mask" "$input_json")
ref_panel_bincov_matrix=$(jq -r ".ref_panel_bincov_matrix" "$input_json")

ref_samples_list=$(jq -r ".ref_samples_list" "$input_json")
mapfile -t ref_samples < "${ref_samples_list}"
ref_samples_json_array=$(printf '%s\n' "${ref_samples[@]}" | jq -R . | jq -s -c .)

ref_pesr_disc_files_list=$(jq -r ".ref_pesr_disc_files_list" "$input_json")
mapfile -t ref_pesr_disc_files < "${ref_pesr_disc_files_list}"
ref_pesr_disc_files_json_array=$(printf '%s\n' "${ref_pesr_disc_files[@]}" | jq -R . | jq -s -c .)

ref_pesr_split_files_list=$(jq -r ".ref_pesr_split_files_list" "$input_json")
mapfile -t ref_pesr_split_files < "${ref_pesr_split_files_list}"
ref_pesr_split_files_json_array=$(printf '%s\n' "${ref_pesr_split_files[@]}" | jq -R . | jq -s -c .)

ref_pesr_sd_files_list=$(jq -r ".ref_pesr_sd_files_list" "$input_json")
mapfile -t ref_pesr_sd_files < "${ref_pesr_sd_files_list}"
ref_pesr_sd_files_json_array=$(printf '%s\n' "${ref_pesr_sd_files[@]}" | jq -R . | jq -s -c .)

gcnv_model_tars_list=$(jq -r ".gcnv_model_tars_list" "$input_json")
mapfile -t gcnv_model_tars < "${gcnv_model_tars_list}"
gcnv_model_tars_json_array=$(printf '%s\n' "${gcnv_model_tars[@]}" | jq -R . | jq -s -c .)


# -------------------------------------------------------
# ======================= Command =======================
# -------------------------------------------------------

# GatherSampleEvidence
# ---------------------------------------------------------------------------------------------------------------------

gather_sample_evidence_output_dir=$(mktemp -d "/output_GatherSampleEvidence_XXXXXXXX")
gather_sample_evidence_output_dir="$(realpath ${gather_sample_evidence_output_dir})"
gather_sample_evidence_inputs_json="${gather_sample_evidence_output_dir}/inputs.json"
gather_sample_evidence_outputs_json="${gather_sample_evidence_output_dir}/outputs.json"

jq -n \
  --slurpfile inputs "${input_json}" \
  '{
    "sample_id": $inputs[0].sample_id,
    "bam_or_cram_file": $inputs[0].bam_or_cram_file,
    "bam_or_cram_index": $inputs[0].bam_or_cram_index,
    "reference_dict": $inputs[0].reference_dict,
    "reference_fasta": $inputs[0].reference_fasta,
    "reference_index": $inputs[0].reference_index,
    "primary_contigs_list": $inputs[0].primary_contigs_list,
    "primary_contigs_fai": $inputs[0].primary_contigs_fai,
    "preprocessed_intervals": $inputs[0].preprocessed_intervals,
    "manta_region_bed": $inputs[0].manta_region_bed,
    "manta_region_bed_index": $inputs[0].manta_region_bed_index,
    "sd_locs_vcf": $inputs[0].sd_locs_vcf,
    "mei_bed": $inputs[0].mei_bed,
    "wham_include_list_bed_file": $inputs[0].wham_include_list_bed_file,
    "reference_bwa_alt": $inputs[0].reference_bwa_alt,
    "reference_bwa_amb": $inputs[0].reference_bwa_amb,
    "reference_bwa_ann": $inputs[0].reference_bwa_ann,
    "reference_bwa_bwt": $inputs[0].reference_bwa_bwt,
    "reference_bwa_pac": $inputs[0].reference_bwa_pac,
    "reference_bwa_sa": $inputs[0].reference_bwa_sa,
    "collect_coverage": true,
    "run_scramble": true,
    "run_manta": true,
    "run_wham": true,
    "collect_pesr": true,
    "scramble_alignment_score_cutoff": 90,
    "run_module_metrics": $inputs[0].run_sampleevidence_metrics
  }' > "${gather_sample_evidence_inputs_json}"

bash /opt/sv_shell/gather_sample_evidence.sh \
  "${gather_sample_evidence_inputs_json}" \
  "${gather_sample_evidence_outputs_json}" \
  "${gather_sample_evidence_output_dir}"



# TODO: TEMP --------- pipelining
# these inputs you get from gather sample evidence.
# WDL: # File case_counts_file_ = select_first([case_counts_file, GatherSampleEvidence.coverage_counts])
coverage_counts="/inputs/NA12878.counts.tsv.gz"



## TODO: the above file should be indexed.
## check if the output of the gather sample evidence pipeline already indexes it or not
## if not, then use the following solution
## The only solution I found that works is the following, but not sure if it is the best.
#SKIP_LINES=$(zcat "${coverage_counts}" | grep -c -E '^@|^CONTIG\s')
#tabix -S $SKIP_LINES -s 1 -b 2 -e 3 "${coverage_counts}"



# EvidenceQC
# ---------------------------------------------------------------------------------------------------------------------
evidence_qc_output_dir=$(mktemp -d "/output_evidence_qc_XXXXXXXX")
evidence_qc_output_dir="$(realpath ${evidence_qc_output_dir})"
evidence_qc_inputs_json_filename="${evidence_qc_output_dir}/inputs.json"
evidence_qc_outputs_json_filename="${evidence_qc_output_dir}/outputs.json"

jq -n \
  --slurpfile inputs "${input_json}" \
  --arg samples "${sample_id}" \
  --arg count_files "${coverage_counts}" \
  '{
      batch: $inputs[0].batch,
      samples: [$samples],
      run_vcf_qc: $inputs[0].run_vcf_qc,
      genome_file: $inputs[0].genome_file,
      count_files: [$count_files],
      run_ploidy: false,
      wgd_scoring_mask: $inputs[0].wgd_scoring_mask,
      reference_dict: $inputs[0].reference_dict,
      bincov_matrix: $inputs[0].ref_panel_bincov_matrix
  }' > "${evidence_qc_inputs_json_filename}"

bash /opt/sv_shell/evidence_qc.sh \
  "${evidence_qc_inputs_json_filename}" \
  "${evidence_qc_outputs_json_filename}" \
  "${evidence_qc_output_dir}"


# TODO: temp pipelining
gather_sample_evidence_outputs_json="/output_GatherSampleEvidence_mZnM8UoL/gather_sample_evidence_outputs.json"
evidence_qc_outputs_json_filename="/output_evidence_qc_SAHS0zt3/output.json"

# GatherBatchEvidence
# ---------------------------------------------------------------------------------------------------------------------
gather_batch_evidence_output_dir=$(mktemp -d "/output_gather_batch_evidence_XXXXXXXX")
gather_batch_evidence_output_dir="$(realpath ${gather_batch_evidence_output_dir})"
gather_batch_evidence_inputs_json_filename="${gather_batch_evidence_output_dir}/inputs.json"
gather_batch_evidence_outputs_json_filename="${gather_batch_evidence_output_dir}/outputs.json"

jq -n \
  --slurpfile inputs "${input_json}" \
  --slurpfile gse_outputs "${gather_sample_evidence_outputs_json}" \
  --slurpfile eqc_outputs "${evidence_qc_outputs_json_filename}" \
  --arg samples "${sample_id}" \
  --argjson ref_samples "${ref_samples_json_array}" \
  --argjson ref_pe_disc "${ref_pesr_disc_files_json_array}" \
  --argjson ref_pe_split "${ref_pesr_split_files_json_array}" \
  --argjson ref_pe_sd "${ref_pesr_sd_files_json_array}" \
  --argjson gcnv_model_tars "${gcnv_model_tars_json_array}" \
  '{
      batch: $inputs[0].batch,
      samples: [$samples],
      ref_panel_samples: $ref_samples,
      run_matrix_qc: false,
      ped_file: $inputs[0].ref_ped_file,
      genome_file: $inputs[0].genome_file,
      primary_contigs_fai: $inputs[0].primary_contigs_fai,
      reference_dict: $inputs[0].reference_dict,
      counts: [$gse_outputs[0].coverage_counts],
      sample_bincov_matrix: $inputs[0].sample_bincov_matrix,
      ref_panel_bincov_matrix: $inputs[0].ref_panel_bincov_matrix,
      bincov_matrix: $eqc_outputs[0].bincov_matrix,
      bincov_matrix_index: $eqc_outputs[0].bincov_matrix_index,
      PE_files: [$gse_outputs[0].pesr_disc],
      cytoband: $inputs[0].cytobands,
      mei_bed: $inputs[0].mei_bed,
      ref_panel_PE_files: $ref_pe_disc,
      SR_files: [$gse_outputs[0].pesr_split],
      ref_panel_SR_files: $ref_pe_split,
      SD_files: [$gse_outputs[0].pesr_sd],
      ref_panel_SD_files: $ref_pe_sd,
      sd_locs_vcf: $inputs[0].sd_locs_vcf,
      contig_ploidy_model_tar: $inputs[0].contig_ploidy_model_tar,
      gcnv_model_tars: $gcnv_model_tars,
      run_ploidy: true,
      append_first_sample_to_ped: true,
      gcnv_p_alt: $inputs[0].gcnv_p_alt,
      gcnv_cnv_coherence_length: $inputs[0].gcnv_cnv_coherence_length,
      gcnv_max_copy_number: $inputs[0].gcnv_max_copy_number,
      gcnv_mapping_error_rate: $inputs[0].gcnv_mapping_error_rate,
      gcnv_sample_psi_scale: $inputs[0].gcnv_sample_psi_scale,
      gcnv_depth_correction_tau: $inputs[0].gcnv_depth_correction_tau,
      gcnv_copy_number_posterior_expectation_mode: $inputs[0].gcnv_copy_number_posterior_expectation_mode,
      gcnv_active_class_padding_hybrid_mode: $inputs[0].gcnv_active_class_padding_hybrid_mode,
      gcnv_learning_rate: $inputs[0].gcnv_learning_rate,
      gcnv_adamax_beta_1: $inputs[0].gcnv_adamax_beta_1,
      gcnv_adamax_beta_2: $inputs[0].gcnv_adamax_beta_2,
      gcnv_log_emission_samples_per_round: $inputs[0].gcnv_log_emission_samples_per_round,
      gcnv_log_emission_sampling_median_rel_error: $inputs[0].gcnv_log_emission_sampling_median_rel_error,
      gcnv_log_emission_sampling_rounds: $inputs[0].gcnv_log_emission_sampling_rounds,
      gcnv_max_advi_iter_first_epoch: $inputs[0].gcnv_max_advi_iter_first_epoch,
      gcnv_max_advi_iter_subsequent_epochs: $inputs[0].gcnv_max_advi_iter_subsequent_epochs,
      gcnv_min_training_epochs: $inputs[0].gcnv_min_training_epochs,
      gcnv_max_training_epochs: $inputs[0].gcnv_max_training_epochs,
      gcnv_initial_temperature: $inputs[0].gcnv_initial_temperature,
      gcnv_num_thermal_advi_iters: $inputs[0].gcnv_num_thermal_advi_iters,
      gcnv_convergence_snr_averaging_window: $inputs[0].gcnv_convergence_snr_averaging_window,
      gcnv_convergence_snr_trigger_threshold: $inputs[0].gcnv_convergence_snr_trigger_threshold,
      gcnv_convergence_snr_countdown_window: $inputs[0].gcnv_convergence_snr_countdown_window,
      gcnv_max_calling_iters: $inputs[0].gcnv_max_calling_iters,
      gcnv_caller_update_convergence_threshold:  $inputs[0].gcnv_caller_update_convergence_threshold,
      gcnv_caller_internal_admixing_rate: $inputs[0].gcnv_caller_internal_admixing_rate,
      gcnv_caller_external_admixing_rate: $inputs[0].gcnv_caller_external_admixing_rate,
      gcnv_disable_annealing: $inputs[0].gcnv_disable_annealing,
      ref_copy_number_autosomal_contigs: $inputs[0].ref_copy_number_autosomal_contigs,
      allosomal_contigs: $inputs[0].allosomal_contigs,
      gcnv_qs_cutoff: $inputs[0].gcnv_qs_cutoff,
      dragen_vcfs: null,
      manta_vcfs: [$gse_outputs[0].manta_vcf],
      scramble_vcfs: [$gse_outputs[0].scramble_vcf],
      wham_vcfs: [$gse_outputs[0].wham_vcf],
      min_svsize: $inputs[0].min_svsize,
      cnmops_chrom_file: $inputs[0].autosome_file,
      cnmops_exclude_list: $inputs[0].cnmops_exclude_list,
      cnmops_allo_file: $inputs[0].allosome_file,
      cnmops_large_min_size: $inputs[0].cnmops_large_min_size,
      matrix_qc_distance: $inputs[0].matrix_qc_distance,
      run_module_metrics: $inputs[0].run_batchevidence_metrics,
      median_cov_mem_gb_per_sample: $inputs[0].median_cov_mem_gb_per_sample
  }' > "${gather_batch_evidence_inputs_json_filename}"

bash /opt/sv_shell/gather_batch_evidence.sh \
  "${gather_batch_evidence_inputs_json_filename}" \
  "${gather_batch_evidence_outputs_json_filename}" \
  "${gather_batch_evidence_output_dir}"

CombineMantaStd_working_dir=$(mktemp -d /wd_CombineMantaStd_XXXXXXXX)
tar xzf "${ref_std_manta_vcf_tar}" -C "${CombineMantaStd_working_dir}/"
tar xzf "${std_manta_vcf_tar}" -C "${CombineMantaStd_working_dir}/"
tar czf "${merged_manta_vcf_tar}" -C "${CombineMantaStd_working_dir}/" .

# CombineScrambleStd
# -----------------------
ref_std_scramble_vcf_tar=$(jq -r ".ref_std_scramble_vcf_tar" "$input_json")
std_scramble_vcf_tar=$(jq -r ".std_scramble_vcf_tar" "$gather_batch_evidence_outputs_json_filename")
merged_scramble_vcf_tar="${working_dir}/$(basename "${std_scramble_vcf_tar}")"

CombineScrambleStd_working_dir=$(mktemp -d /wd_CombineScrambleStd_XXXXXXXX)
tar xzf "${ref_std_scramble_vcf_tar}" -C "${CombineScrambleStd_working_dir}/"
tar xzf "${std_scramble_vcf_tar}" -C "${CombineScrambleStd_working_dir}/"
tar czf "${merged_scramble_vcf_tar}" -C "${CombineScrambleStd_working_dir}/" .

# CombineWhamStd
# -----------------------
ref_std_wham_vcf_tar=$(jq -r ".ref_std_wham_vcf_tar" "$input_json")
std_wham_vcf_tar=$(jq -r ".std_wham_vcf_tar" "$gather_batch_evidence_outputs_json_filename")
merged_wham_vcf_tar="${working_dir}/$(basename "${std_wham_vcf_tar}")"

CombineWhamStd_working_dir=$(mktemp -d /wd_CombineWhamStd_XXXXXXXX)
tar xzf "${ref_std_wham_vcf_tar}" -C "${CombineWhamStd_working_dir}/"
tar xzf "${std_wham_vcf_tar}" -C "${CombineWhamStd_working_dir}/"
tar czf "${merged_wham_vcf_tar}" -C "${CombineWhamStd_working_dir}/" .


# Merge depth
# ----------------------------------------------------------------------------------------------------------------------
# Note that the zcat command called in the following is implemented as a function in merge_depth.sh.
# However, for simplicity (merge_depth is using syntax that simplifies passing local arrays within the bash script),
# the pipe is copy-pasted here.
#
# MergeSetDel
# -----------------------
MergeSetDel_beds=(
    "$(jq -r ".merged_dels" "$gather_batch_evidence_outputs_json_filename")"
    "$(jq -r ".ref_panel_del_bed" "$input_json")"
)

zcat -f "${MergeSetDel_beds[@]}" \
  | sort -k1,1V -k2,2n \
  | awk -v OFS="\t" -v svtype="DEL" -v batch="${batch}" '{$4=batch"_"svtype"_"NR; print}' \
  | cat <(echo -e "#chr\\tstart\\tend\\tname\\tsample\\tsvtype\\tsources") - \
  | bgzip -c > "${batch}.DEL.bed.gz";
tabix -p bed "${batch}.DEL.bed.gz"


# MergeSetDup
# -----------------------
MergeSetDup_beds=(
    "$(jq -r ".merged_dups" "$gather_batch_evidence_outputs_json_filename")"
    "$(jq -r ".ref_panel_dup_bed" "$input_json")"
)

zcat -f "${MergeSetDup_beds[@]}" \
  | sort -k1,1V -k2,2n \
  | awk -v OFS="\t" -v svtype="DUP" -v batch="${batch}" '{$4=batch"_"svtype"_"NR; print}' \
  | cat <(echo -e "#chr\\tstart\\tend\\tname\\tsample\\tsvtype\\tsources") - \
  | bgzip -c > "${batch}.DUP.bed.gz";
tabix -p bed "${batch}.DUP.bed.gz"