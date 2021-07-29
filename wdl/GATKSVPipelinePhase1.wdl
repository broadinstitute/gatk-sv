version 1.0

import "Module00c.wdl" as m00c
import "Module01.wdl" as m01
import "Module02.wdl" as m02
import "Module03.wdl" as m03
import "Structs.wdl"

# One mighty WDL to rule them all...
# Runs Modules 00c, 01, 02, and 03

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
    ## Module 00c
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
    ## Module 01
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
    ## Module 02
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
    ## Module 03
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

    ############################################################
    ## Module metrics parameters for 00c, 01, 02, and 03 metrics
    ############################################################

    # Run module metrics workflow at the end - by default on except for Module00c because of runtime/expense
    Boolean? run_00c_metrics
    Boolean? run_01_metrics
    Boolean? run_02_metrics
    Boolean? run_03_metrics
    File? primary_contigs_list  # required if run_module_metrics = true

  }

  call m00c.Module00c as Module00c {
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
      run_module_metrics = run_00c_metrics,
      primary_contigs_list = primary_contigs_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call m01.Module01 as Module01 {
    input:
      manta_vcfs=Module00c.std_manta_vcf,
      delly_vcfs=Module00c.std_delly_vcf,
      wham_vcfs=Module00c.std_wham_vcf,
      melt_vcfs=Module00c.std_melt_vcf,
      del_bed=Module00c.merged_dels,
      dup_bed=Module00c.merged_dups,
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
      run_module_metrics = run_01_metrics,
      primary_contigs_list = primary_contigs_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call m02.Module02 as Module02 {
    input:
      batch=batch,
      depth_vcf=Module01.depth_vcf,
      melt_vcf=Module01.melt_vcf,
      delly_vcf=Module01.delly_vcf,
      wham_vcf=Module01.wham_vcf,
      manta_vcf=Module01.manta_vcf,
      baf_metrics=select_first([Module00c.merged_BAF]),
      discfile=Module00c.merged_PE,
      coveragefile=Module00c.merged_bincov,
      splitfile=Module00c.merged_SR,
      medianfile=Module00c.median_cov,
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
      run_module_metrics = run_02_metrics,
      primary_contigs_list = primary_contigs_list
  }

  call m03.Module03 as Module03 {
    input:
      batch=batch,
      manta_vcf=Module01.manta_vcf,
      delly_vcf=Module01.delly_vcf,
      wham_vcf=Module01.wham_vcf,
      melt_vcf=Module01.melt_vcf,
      depth_vcf=Module01.depth_vcf,
      outlier_cutoff_table=outlier_cutoff_table,
      evidence_metrics=Module02.metrics,
      evidence_metrics_common=Module02.metrics_common,
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
      run_module_metrics = run_03_metrics,
      primary_contigs_list = primary_contigs_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      ped_file = ped_file
  }

  output {
    # Module 00
    File merged_BAF = select_first([Module00c.merged_BAF])
    File merged_BAF_index = select_first([Module00c.merged_BAF_index])
    File merged_SR = Module00c.merged_SR
    File merged_SR_index = Module00c.merged_SR_index
    File merged_PE = Module00c.merged_PE
    File merged_PE_index = Module00c.merged_PE_index
    File merged_bincov = Module00c.merged_bincov
    File merged_bincov_index = Module00c.merged_bincov_index

    File median_cov = Module00c.median_cov

    File? PE_stats = Module00c.PE_stats
    File? RD_stats = Module00c.RD_stats
    File? SR_stats = Module00c.SR_stats
    File? BAF_stats = Module00c.BAF_stats
    File? Matrix_QC_plot=Module00c.Matrix_QC_plot

    File merged_dels = Module00c.merged_dels
    File merged_dups = Module00c.merged_dups

    Array[File]? std_manta_vcf = Module00c.std_manta_vcf
    Array[File]? std_delly_vcf = Module00c.std_delly_vcf
    Array[File]? std_melt_vcf = Module00c.std_melt_vcf
    Array[File]? std_wham_vcf = Module00c.std_wham_vcf

    File? metrics_file_00c = Module00c.metrics_file_00c

    # Module 01
    File? depth_vcf = Module01.depth_vcf
    File? manta_vcf = Module01.manta_vcf
    File? delly_vcf = Module01.delly_vcf
    File? wham_vcf = Module01.wham_vcf
    File? melt_vcf = Module01.melt_vcf

    File? metrics_file_01 = Module01.metrics_file_01

    # Module 02
    File evidence_metrics = Module02.metrics
    File evidence_metrics_common = Module02.metrics_common

    File? metrics_file_02 = Module02.metrics_file_02

    # Module 03
    File? filtered_manta_vcf = Module03.filtered_manta_vcf
    File? filtered_delly_vcf = Module03.filtered_delly_vcf
    File? filtered_wham_vcf = Module03.filtered_wham_vcf
    File? filtered_melt_vcf = Module03.filtered_melt_vcf
    File? filtered_depth_vcf = Module03.filtered_depth_vcf
    File filtered_pesr_vcf = Module03.filtered_pesr_vcf
    File cutoffs = Module03.cutoffs
    File scores = Module03.scores
    File RF_intermediate_files = Module03.RF_intermediate_files
    Array[String] outlier_samples_excluded = Module03.outlier_samples_excluded
    Array[String] batch_samples_postOutlierExclusion = Module03.batch_samples_postOutlierExclusion
    File outlier_samples_excluded_file = Module03.outlier_samples_excluded_file
    File batch_samples_postOutlierExclusion_file = Module03.batch_samples_postOutlierExclusion_file

    File? metrics_file_03 = Module03.metrics_file_03
  }
}
