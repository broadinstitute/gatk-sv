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
    File reference_index    # Index (.fai), must be in same dir as fasta
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
    Array[File?]? BAF_files
    Array[File] PE_files
    Array[File] SR_files
    Array[File] counts
    File? bincov_matrix
    File? bincov_matrix_index
    File inclusion_bed

    # BAF generation if BAF_files unavailable
    # BAF Option #1, gVCFs
    Array[File]? gvcfs
    File? unpadded_intervals_file
    File? dbsnp_vcf
    File? dbsnp_vcf_index
    File? gvcf_gcs_project_for_requester_pays

    # BAF Option #2, position-sharded VCFs
    Array[File]? snp_vcfs
    File? snp_vcf_header  # Only use if snp vcfs are unheadered

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
    Array[File]? delly_vcfs        # Delly VCF
    Array[File]? melt_vcfs         # Melt VCF
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
    RuntimeAttr? runtime_attr_shard_baf
    RuntimeAttr? runtime_attr_merge_baf
    RuntimeAttr? runtime_attr_shard_pe
    RuntimeAttr? runtime_attr_merge_pe
    RuntimeAttr? runtime_attr_shard_sr
    RuntimeAttr? runtime_attr_merge_sr

    RuntimeAttr? runtime_attr_set_sample
    RuntimeAttr? evidence_merging_bincov_runtime_attr # Disk space ignored, use evidence_merging_bincov_size_mb

    RuntimeAttr? cnmops_sample10_runtime_attr   # Memory ignored if cnmops_mem_gb_override_sample10 given
    RuntimeAttr? cnmops_sample3_runtime_attr    # Memory ignored if cnmops_mem_gb_override_sample3 given
    Float? cnmops_mem_gb_override_sample10
    Float? cnmops_mem_gb_override_sample3

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
    RuntimeAttr? runtime_attr_bundle
    RuntimeAttr? runtime_attr_postprocess
    RuntimeAttr? runtime_attr_explode

    ############################################################
    ## ClusterBatch
    ############################################################

    Int pesr_svsize
    Float pesr_frac
    String pesr_flags
    Int pesr_distance
    File pesr_exclude_list
    String depth_flags
    Float depth_frac

    File? depth_exclude_list
    Float? depth_exclude_list_frac_max

    RuntimeAttr? runtime_attr_pesr_cluster
    RuntimeAttr? runtime_attr_pesr_concat
    RuntimeAttr? runtime_attr_depth_cluster
    RuntimeAttr? runtime_attr_depth_concat
    RuntimeAttr? runtime_attr_depth_vcf
    RuntimeAttr? runtime_attr_rdtest_bed

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
    RuntimeAttr? runtime_attr_merge_baf
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
    RuntimeAttr? runtime_attr_exclude_outliers
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
      gvcfs = gvcfs,
      unpadded_intervals_file = unpadded_intervals_file,
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_index = dbsnp_vcf_index,
      gvcf_gcs_project_for_requester_pays = gvcf_gcs_project_for_requester_pays,
      ref_fasta = reference_fasta,
      ref_fasta_index = reference_index,
      ref_dict = reference_dict,
      snp_vcfs = snp_vcfs,
      snp_vcf_header = snp_vcf_header,
      cytoband = cytoband,
      mei_bed = mei_bed,
      inclusion_bed = inclusion_bed,
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
      delly_vcfs=delly_vcfs,
      melt_vcfs=melt_vcfs,
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
      runtime_attr_set_sample = runtime_attr_set_sample,
      runtime_attr_shard_baf = runtime_attr_shard_baf,
      runtime_attr_merge_baf = runtime_attr_merge_baf,
      runtime_attr_shard_pe = runtime_attr_shard_pe,
      runtime_attr_merge_pe = runtime_attr_merge_pe,
      runtime_attr_shard_sr = runtime_attr_shard_sr,
      runtime_attr_merge_sr = runtime_attr_merge_sr,
      evidence_merging_bincov_runtime_attr=evidence_merging_bincov_runtime_attr,
      cnmops_sample10_runtime_attr=cnmops_sample10_runtime_attr,
      cnmops_sample3_runtime_attr=cnmops_sample3_runtime_attr,
      cnmops_mem_gb_override_sample10=cnmops_mem_gb_override_sample10,
      cnmops_mem_gb_override_sample3=cnmops_mem_gb_override_sample3,
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
      runtime_attr_bundle = runtime_attr_bundle,
      runtime_attr_postprocess = runtime_attr_postprocess,
      runtime_attr_explode = runtime_attr_explode,
      run_module_metrics = run_batchevidence_metrics,
      primary_contigs_list = primary_contigs_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call clusterbatch.ClusterBatch as ClusterBatch {
    input:
      manta_vcfs=GatherBatchEvidence.std_manta_vcf,
      delly_vcfs=GatherBatchEvidence.std_delly_vcf,
      wham_vcfs=GatherBatchEvidence.std_wham_vcf,
      melt_vcfs=GatherBatchEvidence.std_melt_vcf,
      del_bed=GatherBatchEvidence.merged_dels,
      dup_bed=GatherBatchEvidence.merged_dups,
      batch=batch,
      pesr_svsize=pesr_svsize,
      pesr_frac=pesr_frac,
      pesr_flags=pesr_flags,
      pesr_distance=pesr_distance,
      pesr_exclude_list=pesr_exclude_list,
      depth_flags=depth_flags,
      depth_frac=depth_frac,
      contigs=contigs,
      depth_exclude_list=depth_exclude_list,
      depth_exclude_list_frac_max=depth_exclude_list_frac_max,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_pesr_cluster=runtime_attr_pesr_cluster,
      runtime_attr_pesr_concat=runtime_attr_pesr_concat,
      runtime_attr_depth_cluster=runtime_attr_depth_cluster,
      runtime_attr_depth_concat=runtime_attr_depth_concat,
      runtime_attr_depth_vcf=runtime_attr_depth_vcf,
      runtime_attr_rdtest_bed=runtime_attr_rdtest_bed,
      run_module_metrics = run_clusterbatch_metrics,
      primary_contigs_list = primary_contigs_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker, 
      linux_docker = linux_docker
  }

  call batchmetrics.GenerateBatchMetrics as GenerateBatchMetrics {
    input:
      batch=batch,
      depth_vcf=ClusterBatch.depth_vcf,
      melt_vcf=ClusterBatch.melt_vcf,
      delly_vcf=ClusterBatch.delly_vcf,
      wham_vcf=ClusterBatch.wham_vcf,
      manta_vcf=ClusterBatch.manta_vcf,
      baf_metrics=select_first([GatherBatchEvidence.merged_BAF]),
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
      manta_vcf=ClusterBatch.manta_vcf,
      delly_vcf=ClusterBatch.delly_vcf,
      wham_vcf=ClusterBatch.wham_vcf,
      melt_vcf=ClusterBatch.melt_vcf,
      depth_vcf=ClusterBatch.depth_vcf,
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
      runtime_attr_exclude_outliers=runtime_attr_exclude_outliers,
      runtime_attr_cat_outliers=runtime_attr_cat_outliers,
      runtime_attr_filter_samples=runtime_attr_filter_samples,
      run_module_metrics = run_filterbatch_metrics,
      primary_contigs_list = primary_contigs_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      ped_file = ped_file
  }

  output {
    # Module 00
    File merged_BAF = select_first([GatherBatchEvidence.merged_BAF])
    File merged_BAF_index = select_first([GatherBatchEvidence.merged_BAF_index])
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

    Array[File]? std_manta_vcf = GatherBatchEvidence.std_manta_vcf
    Array[File]? std_delly_vcf = GatherBatchEvidence.std_delly_vcf
    Array[File]? std_melt_vcf = GatherBatchEvidence.std_melt_vcf
    Array[File]? std_wham_vcf = GatherBatchEvidence.std_wham_vcf

    File? metrics_file_batchevidence = GatherBatchEvidence.metrics_file_batchevidence

    # ClusterBatch
    File depth_vcf = ClusterBatch.depth_vcf
    File depth_vcf_index = ClusterBatch.depth_vcf_index
    File? manta_vcf = ClusterBatch.manta_vcf
    File? manta_vcf_index = ClusterBatch.manta_vcf_index
    File? delly_vcf = ClusterBatch.delly_vcf
    File? delly_vcf_index = ClusterBatch.delly_vcf_index
    File? wham_vcf = ClusterBatch.wham_vcf
    File? wham_vcf_index = ClusterBatch.wham_vcf_index
    File? melt_vcf = ClusterBatch.melt_vcf
    File? melt_vcf_index = ClusterBatch.melt_vcf_index

    File? metrics_file_clusterbatch = ClusterBatch.metrics_file_clusterbatch

    # GenerateBatchMetrics
    File evidence_metrics = GenerateBatchMetrics.metrics
    File evidence_metrics_common = GenerateBatchMetrics.metrics_common

    File? metrics_file_batchmetrics = GenerateBatchMetrics.metrics_file_batchmetrics

    # FilterBatch
    File? filtered_manta_vcf = FilterBatch.filtered_manta_vcf
    File? filtered_delly_vcf = FilterBatch.filtered_delly_vcf
    File? filtered_wham_vcf = FilterBatch.filtered_wham_vcf
    File? filtered_melt_vcf = FilterBatch.filtered_melt_vcf
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
