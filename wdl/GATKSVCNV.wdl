version 1.0

import "Structs.wdl"
import "CNMOPS.wdl" as cnmops
import "CollectCoverage.wdl" as cov
import "MakeBincovMatrix.wdl" as mbm
import "GermlineCNVCase.wdl" as gcnv
import "MedianCov.wdl" as mc

workflow GATKSVCNV {
  input {
    # Batch info
    String batch
    Array[String] samples

    # Read depth evidence files
    Array[File] counts

    # Global files
    File ped_file
    File ref_fasta_dict

    # Condense read counts
    Int? condense_num_bins
    Int? condense_bin_size

    # gCNV inputs
    File contig_ploidy_model_tar
    Array[File] gcnv_model_tars

    File? gatk4_jar_override
    Float? gcnv_p_alt
    Float? gcnv_cnv_coherence_length
    Int? gcnv_max_copy_number

    Float? gcnv_mapping_error_rate
    Float? gcnv_sample_psi_scale
    Float? gcnv_depth_correction_tau
    String? gcnv_copy_number_posterior_expectation_mode
    Int? gcnv_active_class_padding_hybrid_mode

    Float? gcnv_learning_rate
    Float? gcnv_adamax_beta_1
    Float? gcnv_adamax_beta_2
    Int? gcnv_log_emission_samples_per_round
    Float? gcnv_log_emission_sampling_median_rel_error
    Int? gcnv_log_emission_sampling_rounds
    Int? gcnv_max_advi_iter_first_epoch
    Int? gcnv_max_advi_iter_subsequent_epochs
    Int? gcnv_min_training_epochs
    Int? gcnv_max_training_epochs
    Float? gcnv_initial_temperature
    Int? gcnv_num_thermal_advi_iters
    Int? gcnv_convergence_snr_averaging_window
    Float? gcnv_convergence_snr_trigger_threshold
    Int? gcnv_convergence_snr_countdown_window
    Int? gcnv_max_calling_iters
    Float? gcnv_caller_update_convergence_threshold
    Float? gcnv_caller_internal_admixing_rate
    Float? gcnv_caller_external_admixing_rate
    Boolean? gcnv_disable_annealing

    Int ref_copy_number_autosomal_contigs
    Array[String]? allosomal_contigs

    Int gcnv_qs_cutoff              # QS filtering cutoff

    # CNMops files
    File cnmops_chrom_file
    File cnmops_exclude_list
    File cnmops_allo_file
    Int? cnmops_large_min_size      # minimum size call to be detected by CNMOPS running in large mode

    # Runtime parameters
    String sv_base_mini_docker
    String sv_base_docker
    String sv_pipeline_docker
    String sv_pipeline_qc_docker
    String linux_docker
    String condense_counts_docker
    String gatk_docker
    String cnmops_docker

    RuntimeAttr? evidence_merging_bincov_runtime_attr
    RuntimeAttr? median_cov_runtime_attr
    Float? median_cov_mem_gb_per_sample
    RuntimeAttr? cnmops_sample10_runtime_attr   # Memory ignored if cnmops_mem_gb_override_sample10 given
    RuntimeAttr? cnmops_sample3_runtime_attr    # Memory ignored if cnmops_mem_gb_override_sample3 given
    Float? cnmops_mem_gb_override_sample10
    Float? cnmops_mem_gb_override_sample3

    RuntimeAttr? condense_counts_runtime_attr
    RuntimeAttr? runtime_attr_merge_sample
    RuntimeAttr? cnmops_ped_runtime_attr
    RuntimeAttr? cnmops_clean_runtime_attr

    RuntimeAttr? runtime_attr_ploidy
    RuntimeAttr? runtime_attr_case
    RuntimeAttr? runtime_attr_bundle
    RuntimeAttr? runtime_attr_postprocess
    RuntimeAttr? runtime_attr_explode
  }

  call mbm.MakeBincovMatrix {
    input:
      samples = samples,
      count_files = counts,
      batch = batch,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_base_docker = sv_base_docker,
      runtime_attr_override = evidence_merging_bincov_runtime_attr
  }

  Float median_cov_mem_gb_ = select_first([median_cov_mem_gb_per_sample, 0.5]) * length(samples) + 7.5
  call mc.MedianCov {
    input:
      bincov_matrix = MakeBincovMatrix.merged_bincov,
      cohort_id = batch,
      sv_pipeline_qc_docker = sv_pipeline_qc_docker,
      runtime_attr = median_cov_runtime_attr,
      mem_gb_override = median_cov_mem_gb_
  }

  call cnmops.CNMOPS as CNMOPSSmall {
    input:
      r1 = "3",
      r2 = "10",
      batch = batch,
      samples = samples,
      bincov_matrix = MakeBincovMatrix.merged_bincov,
      bincov_matrix_index = MakeBincovMatrix.merged_bincov_idx,
      chrom_file = cnmops_chrom_file,
      ped_file = ped_file,
      exclude_list = cnmops_exclude_list,
      allo_file = cnmops_allo_file,
      ref_dict = ref_fasta_dict,
      prefix = "small",
      mem_gb_override_sample10 = cnmops_mem_gb_override_sample10,
      mem_gb_override_sample3 = cnmops_mem_gb_override_sample3,
      linux_docker = linux_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      cnmops_docker = cnmops_docker,
      runtime_attr_sample10 = cnmops_sample10_runtime_attr,
      runtime_attr_sample3 = cnmops_sample3_runtime_attr,
      runtime_attr_ped = cnmops_ped_runtime_attr,
      runtime_attr_clean = cnmops_clean_runtime_attr
  }

  call cnmops.CNMOPS as CNMOPSLarge {
    input:
      r1 = "1000",
      r2 = "100",
      batch = batch,
      samples = samples,
      bincov_matrix = MakeBincovMatrix.merged_bincov,
      bincov_matrix_index = MakeBincovMatrix.merged_bincov_idx,
      chrom_file = cnmops_chrom_file,
      ped_file = ped_file,
      exclude_list = cnmops_exclude_list,
      allo_file = cnmops_allo_file,
      ref_dict = ref_fasta_dict,
      prefix = "large",
      min_size=cnmops_large_min_size,
      mem_gb_override_sample10 = cnmops_mem_gb_override_sample10,
      mem_gb_override_sample3 = cnmops_mem_gb_override_sample3,
      linux_docker = linux_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      cnmops_docker = cnmops_docker,
      runtime_attr_sample10 = cnmops_sample10_runtime_attr,
      runtime_attr_sample3 = cnmops_sample3_runtime_attr,
      runtime_attr_ped = cnmops_ped_runtime_attr,
      runtime_attr_clean = cnmops_clean_runtime_attr
  }

  scatter (i in range(length(samples))) {
    call MergeSampleCnmops {
      input:
        sample_id = samples[i],
        cnmops_del_beds = [CNMOPSSmall.Del, CNMOPSLarge.Del],
        cnmops_dup_beds = [CNMOPSSmall.Dup, CNMOPSLarge.Dup],
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_merge_sample
    }
  }

  scatter (i in range(length(samples))) {
    call cov.CondenseReadCounts {
      input:
        counts = counts[i],
        sample = samples[i],
        num_bins = condense_num_bins,
        expected_bin_size = condense_bin_size,
        condense_counts_docker = condense_counts_docker,
        runtime_attr_override=condense_counts_runtime_attr
    }
  }

  call gcnv.CNVGermlineCaseWorkflow {
    input:
      counts = CondenseReadCounts.out,
      count_entity_ids = samples,
      contig_ploidy_model_tar = contig_ploidy_model_tar,
      gcnv_model_tars = gcnv_model_tars,
      gatk_docker = gatk_docker,
      linux_docker = linux_docker,
      sv_base_mini_docker = sv_base_mini_docker,
      gatk4_jar_override = gatk4_jar_override,
      gcnv_p_alt = gcnv_p_alt,
      gcnv_cnv_coherence_length = gcnv_cnv_coherence_length,
      gcnv_max_copy_number = gcnv_max_copy_number,
      gcnv_mapping_error_rate = gcnv_mapping_error_rate,
      gcnv_sample_psi_scale = gcnv_sample_psi_scale,
      gcnv_depth_correction_tau = gcnv_depth_correction_tau,
      gcnv_copy_number_posterior_expectation_mode = gcnv_copy_number_posterior_expectation_mode,
      gcnv_active_class_padding_hybrid_mode = gcnv_active_class_padding_hybrid_mode,
      gcnv_learning_rate = gcnv_learning_rate,
      gcnv_adamax_beta_1 = gcnv_adamax_beta_1,
      gcnv_adamax_beta_2 = gcnv_adamax_beta_2,
      gcnv_log_emission_samples_per_round = gcnv_log_emission_samples_per_round,
      gcnv_log_emission_sampling_median_rel_error = gcnv_log_emission_sampling_median_rel_error,
      gcnv_log_emission_sampling_rounds = gcnv_log_emission_sampling_rounds,
      gcnv_max_advi_iter_first_epoch = gcnv_max_advi_iter_first_epoch,
      gcnv_max_advi_iter_subsequent_epochs = gcnv_max_advi_iter_subsequent_epochs,
      gcnv_min_training_epochs = gcnv_min_training_epochs,
      gcnv_max_training_epochs = gcnv_max_training_epochs,
      gcnv_initial_temperature = gcnv_initial_temperature,
      gcnv_num_thermal_advi_iters = gcnv_num_thermal_advi_iters,
      gcnv_convergence_snr_averaging_window = gcnv_convergence_snr_averaging_window,
      gcnv_convergence_snr_trigger_threshold = gcnv_convergence_snr_trigger_threshold,
      gcnv_convergence_snr_countdown_window = gcnv_convergence_snr_countdown_window,
      gcnv_max_calling_iters = gcnv_max_calling_iters,
      gcnv_caller_update_convergence_threshold = gcnv_caller_update_convergence_threshold,
      gcnv_caller_internal_admixing_rate = gcnv_caller_internal_admixing_rate,
      gcnv_caller_external_admixing_rate = gcnv_caller_external_admixing_rate,
      gcnv_disable_annealing = gcnv_disable_annealing,
      ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
      allosomal_contigs = allosomal_contigs,
      runtime_attr_ploidy = runtime_attr_ploidy,
      runtime_attr_case = runtime_attr_case,
      runtime_attr_bundle = runtime_attr_bundle,
      runtime_attr_postprocess = runtime_attr_postprocess,
      runtime_attr_explode = runtime_attr_explode
  }

  output {
    File medcov = MedianCov.medianCov
    Array[File] cnmops_beds = MergeSampleCnmops.out
    Array[File] cnmops_bed_indexes = MergeSampleCnmops.out_index
    Array[File] gcnv_segments_vcfs = CNVGermlineCaseWorkflow.genotyped_segments_vcf
    Array[File] gcnv_ploidy_call_tars = CNVGermlineCaseWorkflow.sample_contig_ploidy_calls_tars
  }
}

task MergeSampleCnmops {
  input {
    Array[File] cnmops_del_beds
    Array[File] cnmops_dup_beds
    String sample_id
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{sample_id}.cnmops.bed.gz"
    File out_index = "~{sample_id}.cnmops.bed.gz.tbi"
  }
  command <<<

    set -euo pipefail

    # Merge DEL and DUP separately
    zcat ~{sep=" " cnmops_del_beds} | awk -F "\t" -v OFS="\t" '{if ($5=="~{sample_id}") print}' \
      | sort -k1,1V -k2,2n \
      | bedtools merge -i stdin -d 0 -c 4,5,6,7 -o distinct \
      > ~{sample_id}.del.bed
    zcat ~{sep=" " cnmops_dup_beds} | awk -F "\t" -v OFS="\t" '{if ($5=="~{sample_id}") print}' \
      | sort -k1,1V -k2,2n \
      | bedtools merge -i stdin -d 0 -c 4,5,6,7 -o distinct \
      > ~{sample_id}.dup.bed

    cat ~{sample_id}.del.bed ~{sample_id}.dup.bed \
      | sort -k1,1V -k2,2n \
      > ~{sample_id}.cnmops.bed

    bgzip ~{sample_id}.cnmops.bed
    tabix -p bed ~{sample_id}.cnmops.bed.gz

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}
