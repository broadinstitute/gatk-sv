version 1.0

import "GatherBatchEvidence.wdl" as batchevidence
import "ClusterBatch.wdl" as clusterbatch
import "GenerateBatchMetrics.wdl" as batchmetrics
import "FilterBatch.wdl" as filterbatch
import "Structs.wdl"

# One mighty WDL to rule them all...
# Runs GatherBatchEvidence, ClusterBatch, GenerateBatchMetrics, FilterBatch

workflow GATKSVPipelinePhase1 {
  input {

    # Batch info
    String batch
    Array[String] samples

    # Global files
    File ped_file
    File genome_file
    File contigs          # .fai file of included contigs
    File reference_fasta
    File reference_index
    File reference_dict     # Dictionary (.dict), must be in same dir as fasta

    String sv_base_mini_docker
    String sv_base_docker
    String sv_pipeline_base_docker
    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker
    String sv_pipeline_qc_docker
    String linux_docker
    String cnmops_docker
    String gatk_docker
    String? gcnv_gatk_docker
    String condense_counts_docker

    ############################################################
    ## GatherBatchEvidence
    ############################################################

    # PE/SR/BAF/RD files
    # Supply either BAF_files or (SD_files and sd_locs_vcf)
    Array[File?]? BAF_files
    Array[File] PE_files
    Array[File] SR_files
    Array[File]? SD_files
    File? sd_locs_vcf
    Array[File] counts
    File? bincov_matrix
    File? bincov_matrix_index

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

    Float? ploidy_sample_psi_scale
    Int ref_copy_number_autosomal_contigs
    Array[String]? allosomal_contigs

    Int gcnv_qs_cutoff              # QS filtering cutoff

    # SV tool calls
    Array[File]? manta_vcfs        # Manta VCF
    Array[File]? melt_vcfs         # Melt VCF
    Array[File]? scramble_vcfs     # Scramble VCF
    Array[File]? wham_vcfs         # Wham VCF
    Int min_svsize                 # Minimum SV length to include

    # QC options
    Boolean run_matrix_qc

    # CNMops files
    File cnmops_chrom_file
    File cnmops_exclude_list
    File cnmops_allo_file
    Int? cnmops_large_min_size

    # Resolve files
    File cytoband
    File mei_bed

    # QC files
    Int matrix_qc_distance

    # Runtime parameters
    RuntimeAttr? runtime_attr_merge_baf
    RuntimeAttr? runtime_attr_bem

    RuntimeAttr? evidence_merging_bincov_runtime_attr # Disk space ignored, use evidence_merging_bincov_size_mb

    RuntimeAttr? cnmops_sample10_runtime_attr
    RuntimeAttr? cnmops_sample3_runtime_attr

    RuntimeAttr? median_cov_runtime_attr        # Memory ignored, use median_cov_mem_gb_per_sample
    Float? median_cov_mem_gb_per_sample

    RuntimeAttr? preprocess_calls_runtime_attr
    RuntimeAttr? depth_merge_set_runtime_attr
    RuntimeAttr? depth_merge_sample_runtime_attr
    RuntimeAttr? cnmops_ped_runtime_attr
    RuntimeAttr? cnmops_clean_runtime_attr
    RuntimeAttr? matrix_qc_pesrbaf_runtime_attr
    RuntimeAttr? matrix_qc_rd_runtime_attr

    RuntimeAttr? runtime_attr_ploidy
    RuntimeAttr? runtime_attr_case
    RuntimeAttr? runtime_attr_postprocess
    RuntimeAttr? runtime_attr_explode

    ############################################################
    ## ClusterBatch
    ############################################################
    String? chr_x
    String? chr_y

    Int? depth_records_per_bed_shard_cluster_batch
    File depth_exclude_intervals
    Float depth_exclude_overlap_fraction
    Float depth_interval_overlap
    String? depth_clustering_algorithm

    Int? pesr_min_size
    File pesr_exclude_intervals
    Float pesr_interval_overlap
    Int pesr_breakend_window
    String? pesr_clustering_algorithm

    Int? N_IQR_cutoff_plotting

    File? baseline_depth_vcf_cluster_batch
    File? baseline_manta_vcf_cluster_batch
    File? baseline_wham_vcf_cluster_batch
    File? baseline_melt_vcf_cluster_batch

    Float? java_mem_fraction_cluster_batch

    RuntimeAttr? runtime_attr_ids_from_vcf_list_cluster_batch
    RuntimeAttr? runtime_attr_create_ploidy_cluster_batch
    RuntimeAttr? runtime_attr_prepare_pesr_vcfs_cluster_batch
    RuntimeAttr? runtime_attr_svcluster_manta_cluster_batch
    RuntimeAttr? runtime_attr_svcluster_melt_cluster_batch
    RuntimeAttr? runtime_attr_svcluster_scramble_cluster_batch
    RuntimeAttr? runtime_attr_svcluster_wham_cluster_batch
    RuntimeAttr? runtime_override_concat_vcfs_pesr_cluster_batch
    RuntimeAttr? runtime_attr_gatk_to_svtk_vcf_pesr_cluster_batch
    RuntimeAttr? runtime_attr_scatter_bed_cluster_batch
    RuntimeAttr? runtime_attr_cnv_bed_to_gatk_vcf_cluster_batch
    RuntimeAttr? runtime_attr_exclude_intervals_depth_cluster_batch
    RuntimeAttr? runtime_attr_svcluster_depth_cluster_batch
    RuntimeAttr? runtime_attr_gatk_to_svtk_vcf_depth_cluster_batch
    RuntimeAttr? runtime_override_concat_vcfs_depth_cluster_batch
    RuntimeAttr? runtime_attr_exclude_intervals_pesr_cluster_batch
    RuntimeAttr? runtime_attr_count_svs
    RuntimeAttr? runtime_attr_plot_svcounts
    RuntimeAttr? runtime_attr_cat_outliers_preview

    ############################################################
    ## GenerateBatchMetrics
    ############################################################

    Int BAF_split_size
    Int RD_split_size
    Int PE_split_size
    Int SR_split_size
    Int common_cnv_size_cutoff

    File rmsk
    File segdups
    File autosome_contigs
    File allosome_contigs

    RuntimeAttr? runtime_attr_sample_list
    RuntimeAttr? runtime_attr_petest
    RuntimeAttr? runtime_attr_srtest
    RuntimeAttr? runtime_attr_rdtest
    RuntimeAttr? runtime_attr_baftest
    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_split_rd_vcf
    RuntimeAttr? runtime_attr_split_baf_vcf
    RuntimeAttr? runtime_attr_merge_allo
    RuntimeAttr? runtime_attr_merge_stats

    ############################################################
    ## FilterBatch
    ############################################################

    File? outlier_cutoff_table
    Int outlier_cutoff_nIQR

    RuntimeAttr? runtime_attr_adjudicate
    RuntimeAttr? runtime_attr_filter_annotate_vcf
    RuntimeAttr? runtime_attr_merge_pesr_vcfs
    RuntimeAttr? runtime_attr_identify_outliers
    RuntimeAttr? runtime_attr_subset_vcf
    RuntimeAttr? runtime_attr_cat_outliers
    RuntimeAttr? runtime_attr_filter_samples
    RuntimeAttr? runtime_attr_get_male_only

    ############################################################
    ## Module metrics parameters for GatherBatchEvidence, ClusterBatch, GenerateBatchMetrics, FilterBatch metrics
    ############################################################

    # Run module metrics workflow at the end - by default on except for GatherBatchEvidence because of runtime/expense
    Boolean? run_batchevidence_metrics
    Boolean? run_clusterbatch_metrics
    Boolean? run_batchmetrics_metrics
    Boolean? run_filterbatch_metrics
    File? primary_contigs_list  # required if run_module_metrics = true

  }

  call batchevidence.GatherBatchEvidence as GatherBatchEvidence {
    input:
      batch = batch,
      samples = samples,
      ped_file = ped_file,
      genome_file = genome_file,
      primary_contigs_fai = contigs,
      BAF_files = BAF_files,
      PE_files = PE_files,
      SR_files = SR_files,
      SD_files = SD_files,
      sd_locs_vcf = sd_locs_vcf,
      ref_dict = reference_dict,
      cytoband = cytoband,
      mei_bed = mei_bed,
      counts = counts,
      bincov_matrix = bincov_matrix,
      bincov_matrix_index = bincov_matrix_index,
      contig_ploidy_model_tar = contig_ploidy_model_tar,
      gcnv_model_tars = gcnv_model_tars,
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
      gcnv_qs_cutoff=gcnv_qs_cutoff,
      manta_vcfs=manta_vcfs,
      melt_vcfs=melt_vcfs,
      scramble_vcfs=scramble_vcfs,
      wham_vcfs=wham_vcfs,
      min_svsize=min_svsize,
      run_matrix_qc=run_matrix_qc,
      cnmops_chrom_file=cnmops_chrom_file,
      cnmops_exclude_list=cnmops_exclude_list,
      cnmops_allo_file=cnmops_allo_file,
      cnmops_large_min_size=cnmops_large_min_size,
      matrix_qc_distance=matrix_qc_distance,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_base_docker=sv_base_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_qc_docker=sv_pipeline_qc_docker,
      linux_docker=linux_docker,
      cnmops_docker=cnmops_docker,
      gatk_docker=gatk_docker,
      gcnv_gatk_docker=gcnv_gatk_docker,
      condense_counts_docker=condense_counts_docker,
      evidence_merging_bincov_runtime_attr=evidence_merging_bincov_runtime_attr,
      cnmops_sample10_runtime_attr=cnmops_sample10_runtime_attr,
      cnmops_sample3_runtime_attr=cnmops_sample3_runtime_attr,
      median_cov_runtime_attr=median_cov_runtime_attr,
      median_cov_mem_gb_per_sample=median_cov_mem_gb_per_sample,
      preprocess_calls_runtime_attr=preprocess_calls_runtime_attr,
      depth_merge_set_runtime_attr=depth_merge_set_runtime_attr,
      depth_merge_sample_runtime_attr=depth_merge_sample_runtime_attr,
      cnmops_ped_runtime_attr=cnmops_ped_runtime_attr,
      cnmops_clean_runtime_attr=cnmops_clean_runtime_attr,
      matrix_qc_pesrbaf_runtime_attr=matrix_qc_pesrbaf_runtime_attr,
      matrix_qc_rd_runtime_attr=matrix_qc_rd_runtime_attr,
      runtime_attr_ploidy = runtime_attr_ploidy,
      runtime_attr_case = runtime_attr_case,
      runtime_attr_postprocess = runtime_attr_postprocess,
      runtime_attr_explode = runtime_attr_explode,
      run_module_metrics = run_batchevidence_metrics,
      primary_contigs_list = primary_contigs_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call clusterbatch.ClusterBatch as ClusterBatch {
    input:
      manta_vcf_tar=GatherBatchEvidence.std_manta_vcf_tar,
      wham_vcf_tar=GatherBatchEvidence.std_wham_vcf_tar,
      scramble_vcf_tar=GatherBatchEvidence.std_scramble_vcf_tar,
      melt_vcf_tar=GatherBatchEvidence.std_melt_vcf_tar,
      del_bed=GatherBatchEvidence.merged_dels,
      dup_bed=GatherBatchEvidence.merged_dups,
      batch=batch,
      ped_file=ped_file,
      contig_list=contigs,
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_index,
      reference_dict=reference_dict,
      chr_x=chr_x,
      chr_y=chr_y,
      depth_records_per_bed_shard=depth_records_per_bed_shard_cluster_batch,
      depth_exclude_intervals=depth_exclude_intervals,
      depth_exclude_overlap_fraction=depth_exclude_overlap_fraction,
      depth_interval_overlap=depth_interval_overlap,
      depth_clustering_algorithm=depth_clustering_algorithm,
      pesr_min_size=pesr_min_size,
      pesr_exclude_intervals=pesr_exclude_intervals,
      pesr_interval_overlap=pesr_interval_overlap,
      pesr_breakend_window=pesr_breakend_window,
      pesr_clustering_algorithm=pesr_clustering_algorithm,
      N_IQR_cutoff_plotting = N_IQR_cutoff_plotting,
      run_module_metrics=run_clusterbatch_metrics,
      linux_docker=linux_docker,
      sv_pipeline_base_docker=sv_pipeline_base_docker,
      baseline_depth_vcf=baseline_depth_vcf_cluster_batch,
      baseline_manta_vcf=baseline_manta_vcf_cluster_batch,
      baseline_wham_vcf=baseline_wham_vcf_cluster_batch,
      baseline_melt_vcf=baseline_melt_vcf_cluster_batch,
      gatk_docker=gatk_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      java_mem_fraction=java_mem_fraction_cluster_batch,
      runtime_attr_ids_from_vcf_list=runtime_attr_ids_from_vcf_list_cluster_batch,
      runtime_attr_create_ploidy=runtime_attr_create_ploidy_cluster_batch,
      runtime_attr_prepare_pesr_vcfs=runtime_attr_prepare_pesr_vcfs_cluster_batch,
      runtime_attr_svcluster_manta=runtime_attr_svcluster_manta_cluster_batch,
      runtime_attr_svcluster_melt=runtime_attr_svcluster_melt_cluster_batch,
      runtime_attr_svcluster_scramble=runtime_attr_svcluster_scramble_cluster_batch,
      runtime_attr_svcluster_wham=runtime_attr_svcluster_wham_cluster_batch,
      runtime_override_concat_vcfs_pesr=runtime_override_concat_vcfs_pesr_cluster_batch,
      runtime_attr_gatk_to_svtk_vcf_pesr=runtime_attr_gatk_to_svtk_vcf_pesr_cluster_batch,
      runtime_attr_scatter_bed=runtime_attr_scatter_bed_cluster_batch,
      runtime_attr_cnv_bed_to_gatk_vcf=runtime_attr_cnv_bed_to_gatk_vcf_cluster_batch,
      runtime_attr_exclude_intervals_depth=runtime_attr_exclude_intervals_depth_cluster_batch,
      runtime_attr_svcluster_depth=runtime_attr_svcluster_depth_cluster_batch,
      runtime_attr_gatk_to_svtk_vcf_depth=runtime_attr_gatk_to_svtk_vcf_depth_cluster_batch,
      runtime_override_concat_vcfs_depth=runtime_override_concat_vcfs_depth_cluster_batch,
      runtime_attr_exclude_intervals_pesr=runtime_attr_exclude_intervals_pesr_cluster_batch,
      runtime_attr_count_svs = runtime_attr_count_svs,
      runtime_attr_plot_svcounts = runtime_attr_plot_svcounts,
      runtime_attr_cat_outliers_preview = runtime_attr_cat_outliers_preview
  }

  call batchmetrics.GenerateBatchMetrics as GenerateBatchMetrics {
    input:
      batch=batch,
      depth_vcf=ClusterBatch.clustered_depth_vcf,
      melt_vcf=ClusterBatch.clustered_melt_vcf,
      scramble_vcf=ClusterBatch.clustered_scramble_vcf,
      wham_vcf=ClusterBatch.clustered_wham_vcf,
      manta_vcf=ClusterBatch.clustered_manta_vcf,
      baf_metrics=GatherBatchEvidence.merged_BAF,
      discfile=GatherBatchEvidence.merged_PE,
      coveragefile=GatherBatchEvidence.merged_bincov,
      splitfile=GatherBatchEvidence.merged_SR,
      medianfile=GatherBatchEvidence.median_cov,
      BAF_split_size=BAF_split_size,
      RD_split_size=RD_split_size,
      PE_split_size=PE_split_size,
      SR_split_size=SR_split_size,
      common_cnv_size_cutoff=common_cnv_size_cutoff,
      ref_dict=reference_dict,
      rmsk=rmsk,
      segdups=segdups,
      ped_file=ped_file,
      autosome_contigs=autosome_contigs,
      allosome_contigs=allosome_contigs,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_base_docker=sv_base_docker,
      sv_pipeline_base_docker=sv_pipeline_base_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
      linux_docker=linux_docker,
      runtime_attr_sample_list=runtime_attr_sample_list,
      runtime_attr_petest=runtime_attr_petest,
      runtime_attr_srtest=runtime_attr_srtest,
      runtime_attr_rdtest=runtime_attr_rdtest,
      runtime_attr_baftest=runtime_attr_baftest,
      runtime_attr_split_vcf=runtime_attr_split_vcf,
      runtime_attr_split_rd_vcf=runtime_attr_split_rd_vcf,
      runtime_attr_split_baf_vcf=runtime_attr_split_baf_vcf,
      runtime_attr_merge_allo=runtime_attr_merge_allo,
      runtime_attr_merge_baf=runtime_attr_merge_baf,
      runtime_attr_merge_stats=runtime_attr_merge_stats,
      run_module_metrics = run_batchmetrics_metrics,
      primary_contigs_list = primary_contigs_list
  }

  call filterbatch.FilterBatch as FilterBatch {
    input:
      batch=batch,
      manta_vcf=ClusterBatch.clustered_manta_vcf,
      wham_vcf=ClusterBatch.clustered_wham_vcf,
      melt_vcf=ClusterBatch.clustered_melt_vcf,
      scramble_vcf=ClusterBatch.clustered_scramble_vcf,
      depth_vcf=ClusterBatch.clustered_depth_vcf,
      outlier_cutoff_table=outlier_cutoff_table,
      evidence_metrics=GenerateBatchMetrics.metrics,
      evidence_metrics_common=GenerateBatchMetrics.metrics_common,
      outlier_cutoff_nIQR=outlier_cutoff_nIQR,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      linux_docker=linux_docker,
      runtime_attr_adjudicate=runtime_attr_adjudicate,
      runtime_attr_filter_annotate_vcf=runtime_attr_filter_annotate_vcf,
      runtime_attr_merge_pesr_vcfs=runtime_attr_merge_pesr_vcfs,
      runtime_attr_identify_outliers=runtime_attr_identify_outliers,
      runtime_attr_subset_vcf=runtime_attr_subset_vcf,
      runtime_attr_cat_outliers=runtime_attr_cat_outliers,
      runtime_attr_filter_samples=runtime_attr_filter_samples,
      run_module_metrics = run_filterbatch_metrics,
      primary_contigs_list = primary_contigs_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      ped_file = ped_file
  }

  output {
    # Module 00
    File merged_BAF = GatherBatchEvidence.merged_BAF
    File merged_BAF_index = GatherBatchEvidence.merged_BAF_index
    File merged_SR = GatherBatchEvidence.merged_SR
    File merged_SR_index = GatherBatchEvidence.merged_SR_index
    File merged_PE = GatherBatchEvidence.merged_PE
    File merged_PE_index = GatherBatchEvidence.merged_PE_index
    File merged_bincov = GatherBatchEvidence.merged_bincov
    File merged_bincov_index = GatherBatchEvidence.merged_bincov_index

    File median_cov = GatherBatchEvidence.median_cov

    File? PE_stats = GatherBatchEvidence.PE_stats
    File? RD_stats = GatherBatchEvidence.RD_stats
    File? SR_stats = GatherBatchEvidence.SR_stats
    File? BAF_stats = GatherBatchEvidence.BAF_stats
    File? Matrix_QC_plot=GatherBatchEvidence.Matrix_QC_plot

    File merged_dels = GatherBatchEvidence.merged_dels
    File merged_dups = GatherBatchEvidence.merged_dups

    File? std_manta_vcf_tar = GatherBatchEvidence.std_manta_vcf_tar
    File? std_melt_vcf_tar = GatherBatchEvidence.std_melt_vcf_tar
    File? std_scramble_vcf_tar = GatherBatchEvidence.std_scramble_vcf_tar
    File? std_wham_vcf_tar = GatherBatchEvidence.std_wham_vcf_tar

    File? metrics_file_batchevidence = GatherBatchEvidence.metrics_file_batchevidence

    # ClusterBatch
    File depth_vcf = ClusterBatch.clustered_depth_vcf
    File depth_vcf_index = ClusterBatch.clustered_depth_vcf_index
    File? manta_vcf = ClusterBatch.clustered_manta_vcf
    File? manta_vcf_index = ClusterBatch.clustered_manta_vcf_index
    File? wham_vcf = ClusterBatch.clustered_wham_vcf
    File? wham_vcf_index = ClusterBatch.clustered_wham_vcf_index
    File? melt_vcf = ClusterBatch.clustered_melt_vcf
    File? melt_vcf_index = ClusterBatch.clustered_melt_vcf_index
    File? scramble_vcf = ClusterBatch.clustered_scramble_vcf
    File? scramble_vcf_index = ClusterBatch.clustered_scramble_vcf_index
    Array[File]? clustered_sv_counts = ClusterBatch.clustered_sv_counts
    Array[File]? clustered_sv_count_plots = ClusterBatch.clustered_sv_count_plots
    File? clustered_outlier_samples_preview = ClusterBatch.clustered_outlier_samples_preview
    File? clustered_outlier_samples_with_reason = ClusterBatch.clustered_outlier_samples_with_reason
    Int? clustered_num_outlier_samples = ClusterBatch.clustered_num_outlier_samples

    File? metrics_file_clusterbatch = ClusterBatch.metrics_file_clusterbatch

    # GenerateBatchMetrics
    File evidence_metrics = GenerateBatchMetrics.metrics
    File evidence_metrics_common = GenerateBatchMetrics.metrics_common

    File? metrics_file_batchmetrics = GenerateBatchMetrics.metrics_file_batchmetrics

    # FilterBatch
    File? filtered_manta_vcf = FilterBatch.filtered_manta_vcf
    File? filtered_wham_vcf = FilterBatch.filtered_wham_vcf
    File? filtered_melt_vcf = FilterBatch.filtered_melt_vcf
    File? filtered_scramble_vcf = FilterBatch.filtered_scramble_vcf
    File? filtered_depth_vcf = FilterBatch.filtered_depth_vcf
    File? filtered_pesr_vcf = FilterBatch.filtered_pesr_vcf
    File cutoffs = FilterBatch.cutoffs
    File scores = FilterBatch.scores
    File RF_intermediate_files = FilterBatch.RF_intermediate_files
    Array[String] outlier_samples_excluded = FilterBatch.outlier_samples_excluded
    Array[String] batch_samples_postOutlierExclusion = FilterBatch.batch_samples_postOutlierExclusion
    File outlier_samples_excluded_file = FilterBatch.outlier_samples_excluded_file
    File batch_samples_postOutlierExclusion_file = FilterBatch.batch_samples_postOutlierExclusion_file

    File? sites_filtered_manta_vcf = FilterBatch.sites_filtered_manta_vcf
    File? sites_filtered_wham_vcf = FilterBatch.sites_filtered_wham_vcf
    File? sites_filtered_melt_vcf = FilterBatch.sites_filtered_melt_vcf
    File? sites_filtered_depth_vcf = FilterBatch.sites_filtered_depth_vcf

    File? metrics_file_filterbatch = FilterBatch.metrics_file_filterbatch
  }
}
