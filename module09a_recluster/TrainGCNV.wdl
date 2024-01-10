version 1.0

import "Structs.wdl"
import "CollectCoverage.wdl" as cov
import "GermlineCNVCohort.wdl" as gcnv_cohort
import "Utils.wdl" as util

# Trains gCNV model on a cohort with counts already collected
workflow TrainGCNV {
  input {
    Array[String] samples
    Array[File] count_files

    # Common parameters
    String cohort
    File reference_fasta
    File reference_index    # Index (.fai), must be in same dir as fasta
    File reference_dict     # Dictionary (.dict), must be in same dir as fasta

    # Options for subsetting samples for training. Both options require providing sv_pipeline_base_docker
    # Assumes all other inputs correspond to the full sample list. Intended for Terra
    Int? n_samples_subsample # Number of samples to subsample from provided sample list for trainGCNV (rec: ~100)
    Int subsample_seed = 42
    # Subset of full sample list on which to train the gCNV model. Overrides n_samples_subsample if both provided
    Array[String]? sample_ids_training_subset

    # Condense read counts
    Int? min_interval_size
    Int? max_interval_size

    # gCNV common inputs
    Int ref_copy_number_autosomal_contigs
    Array[String]? allosomal_contigs

    # Interval filtering inputs
    File? exclude_intervals_for_filter_intervals_ploidy
    File? exclude_intervals_for_filter_intervals_cnv

    # gCNV cohort mode inputs
    Boolean? filter_intervals
    File contig_ploidy_priors
    Int num_intervals_per_scatter
    Boolean? do_explicit_gc_correction
    Boolean? gcnv_enable_bias_factors

    # gCNV additional arguments
    Float? gcnv_learning_rate
    Int? gcnv_num_thermal_advi_iters
    Int? gcnv_max_advi_iter_first_epoch
    Int? gcnv_max_advi_iter_subsequent_epochs
    Int? gcnv_max_training_epochs
    Int? gcnv_min_training_epochs
    Int? gcnv_convergence_snr_averaging_window
    Int? gcnv_convergence_snr_countdown_window
    Int? gcnv_cnv_coherence_length
    String? gcnv_copy_number_posterior_expectation_mode
    Int? gcnv_log_emission_sampling_rounds
    Float? gcnv_p_alt
    Float? gcnv_sample_psi_scale
    Float? ploidy_sample_psi_scale

    Float? gcnv_caller_update_convergence_threshold
    Float? gcnv_class_coherence_length
    Float? gcnv_convergence_snr_trigger_threshold
    Float? gcnv_interval_psi_scale
    Float? gcnv_log_emission_sampling_median_rel_error
    Float? gcnv_log_mean_bias_standard_deviation
    Int? gcnv_max_bias_factors
    Int? gcnv_max_calling_iters
    Float? ploidy_global_psi_scale
    Float? ploidy_mean_bias_standard_deviation
    Float? gcnv_depth_correction_tau

    # gCNV model building arguments
    Float? gcnv_model_learning_rate
    Int? gcnv_model_num_thermal_advi_iters
    Int? gcnv_model_max_advi_iter_first_epoch
    Int? gcnv_model_max_advi_iter_subsequent_epochs
    Int? gcnv_model_max_training_epochs
    Int? gcnv_model_min_training_epochs
    Int? gcnv_model_convergence_snr_averaging_window
    Int? gcnv_model_convergence_snr_countdown_window
    Int? gcnv_model_cnv_coherence_length
    Int? gcnv_model_log_emission_sampling_rounds

    # Docker
    String sv_base_mini_docker
    String linux_docker
    String gatk_docker
    String condense_counts_docker
    String? sv_pipeline_base_docker # required if using n_samples_subsample or sample_ids_training_subset to subset samples

    # Runtime configuration overrides
    RuntimeAttr? condense_counts_runtime_attr
    RuntimeAttr? counts_to_intervals_runtime_attr
    RuntimeAttr? runtime_attr_annotate
    RuntimeAttr? runtime_attr_filter
    RuntimeAttr? runtime_attr_scatter
    RuntimeAttr? runtime_attr_ploidy
    RuntimeAttr? runtime_attr_cohort
    RuntimeAttr? runtime_attr_postprocess
    RuntimeAttr? runtime_attr_explode
  }

  if (defined(sample_ids_training_subset)) {
    call util.GetSubsampledIndices {
      input:
        all_strings = write_lines(samples),
        subset_strings = write_lines(select_first([sample_ids_training_subset])),
        prefix = cohort,
        sv_pipeline_base_docker = select_first([sv_pipeline_base_docker])
    }
  }

  if (defined(n_samples_subsample) && !defined(sample_ids_training_subset)) {
    call util.RandomSubsampleStringArray {
      input:
        strings = write_lines(samples),
        seed = subsample_seed,
        subset_size = select_first([n_samples_subsample]),
        prefix = cohort,
        sv_pipeline_base_docker = select_first([sv_pipeline_base_docker])
    }
  }

  Array[Int] sample_indices = select_first([GetSubsampledIndices.subsample_indices_array, RandomSubsampleStringArray.subsample_indices_array, range(length(samples))])

  scatter (i in sample_indices) {
    String sample_ids_ = samples[i]
    call cov.CondenseReadCounts as CondenseReadCounts {
      input:
        counts = count_files[i],
        sample = samples[i],
        min_interval_size = min_interval_size,
        max_interval_size = max_interval_size,
        condense_counts_docker = condense_counts_docker,
        runtime_attr_override=condense_counts_runtime_attr
    }
  }

  call cov.CountsToIntervals {
    input:
      counts = CondenseReadCounts.out[0],
      output_name = "condensed_intervals",
      linux_docker = linux_docker,
      runtime_attr_override = counts_to_intervals_runtime_attr
  }

  call gcnv_cohort.CNVGermlineCohortWorkflow {
    input:
      preprocessed_intervals = CountsToIntervals.out,
      filter_intervals = filter_intervals,
      counts = CondenseReadCounts.out,
      count_entity_ids = sample_ids_,
      cohort_entity_id = cohort,
      contig_ploidy_priors = contig_ploidy_priors,
      num_intervals_per_scatter = num_intervals_per_scatter,
      ref_fasta_dict = reference_dict,
      ref_fasta_fai = reference_index,
      ref_fasta = reference_fasta,
      exclude_intervals_for_filter_intervals_ploidy=exclude_intervals_for_filter_intervals_ploidy,
      exclude_intervals_for_filter_intervals_cnv=exclude_intervals_for_filter_intervals_cnv,
      do_explicit_gc_correction = do_explicit_gc_correction,
      gcnv_enable_bias_factors = gcnv_enable_bias_factors,
      ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
      allosomal_contigs = allosomal_contigs,
      gatk_docker = gatk_docker,
      linux_docker = linux_docker,
      sv_base_mini_docker = sv_base_mini_docker,
      gcnv_learning_rate = gcnv_learning_rate,
      gcnv_max_advi_iter_first_epoch = gcnv_max_advi_iter_first_epoch,
      gcnv_num_thermal_advi_iters = gcnv_model_num_thermal_advi_iters,
      gcnv_max_advi_iter_subsequent_epochs = gcnv_model_max_advi_iter_subsequent_epochs,
      gcnv_max_training_epochs = gcnv_model_max_training_epochs,
      gcnv_min_training_epochs = gcnv_model_min_training_epochs,
      gcnv_convergence_snr_averaging_window = gcnv_model_convergence_snr_averaging_window,
      gcnv_convergence_snr_countdown_window = gcnv_model_convergence_snr_countdown_window,
      gcnv_cnv_coherence_length = gcnv_model_cnv_coherence_length,
      gcnv_class_coherence_length = gcnv_class_coherence_length,
      gcnv_copy_number_posterior_expectation_mode = gcnv_copy_number_posterior_expectation_mode,
      gcnv_log_emission_sampling_rounds = gcnv_model_log_emission_sampling_rounds,
      gcnv_p_alt = gcnv_p_alt,
      gcnv_sample_psi_scale = gcnv_sample_psi_scale,
      ploidy_sample_psi_scale = ploidy_sample_psi_scale,
      gcnv_caller_update_convergence_threshold = gcnv_caller_update_convergence_threshold,
      gcnv_convergence_snr_trigger_threshold = gcnv_convergence_snr_trigger_threshold,
      gcnv_interval_psi_scale = gcnv_interval_psi_scale,
      gcnv_log_emission_sampling_median_rel_error = gcnv_log_emission_sampling_median_rel_error,
      gcnv_log_mean_bias_standard_deviation = gcnv_log_mean_bias_standard_deviation,
      gcnv_max_bias_factors = gcnv_max_bias_factors,
      gcnv_max_calling_iters = gcnv_max_calling_iters,
      ploidy_global_psi_scale = ploidy_global_psi_scale,
      ploidy_mean_bias_standard_deviation = ploidy_mean_bias_standard_deviation,
      gcnv_depth_correction_tau = gcnv_depth_correction_tau,
      runtime_attr_annotate = runtime_attr_annotate,
      runtime_attr_filter = runtime_attr_filter,
      runtime_attr_scatter = runtime_attr_scatter,
      runtime_attr_ploidy = runtime_attr_ploidy,
      runtime_attr_cohort = runtime_attr_cohort,
      runtime_attr_postprocess = runtime_attr_postprocess,
      runtime_attr_explode = runtime_attr_explode
  }

  output {
    File? annotated_intervals = CNVGermlineCohortWorkflow.annotated_intervals
    File? filtered_intervals_cnv = CNVGermlineCohortWorkflow.filtered_intervals_cnv
    File? filtered_intervals_ploidy = CNVGermlineCohortWorkflow.filtered_intervals_ploidy
    File? cohort_contig_ploidy_model_tar = CNVGermlineCohortWorkflow.contig_ploidy_model_tar
    File? cohort_contig_ploidy_calls_tar = CNVGermlineCohortWorkflow.contig_ploidy_calls_tar
    Array[File]? cohort_gcnv_model_tars = CNVGermlineCohortWorkflow.gcnv_model_tars
    Array[Array[File]]? cohort_gcnv_calls_tars = CNVGermlineCohortWorkflow.gcnv_calls_tars
    Array[File]? cohort_gcnv_tracking_tars = CNVGermlineCohortWorkflow.gcnv_tracking_tars
    Array[File]? cohort_genotyped_intervals_vcfs = CNVGermlineCohortWorkflow.genotyped_intervals_vcfs
    Array[File]? cohort_genotyped_segments_vcfs = CNVGermlineCohortWorkflow.genotyped_segments_vcfs
    Array[File]? cohort_denoised_copy_ratios = CNVGermlineCohortWorkflow.denoised_copy_ratios
  }
}
