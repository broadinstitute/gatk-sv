version 1.0

import "GatherSampleEvidence.wdl" as sampleevidence
import "EvidenceQC.wdl" as evidenceqc
import "PloidyEstimation.wdl" as pe
import "GatherBatchEvidence.wdl" as batchevidence
import "DepthPreprocessing.wdl" as dpn
import "ClusterBatch.wdl" as clusterbatch
import "GenerateBatchMetrics.wdl" as batchmetrics
import "SRTest.wdl" as SRTest
import "FilterBatchSamples.wdl" as filterbatch
import "GenotypeBatch.wdl" as genotypebatch
import "MakeCohortVcf.wdl" as makecohortvcf
import "TasksMakeCohortVcf.wdl" as tasks_makecohortvcf
import "AnnotateVcf.wdl" as annotate
import "GermlineCNVCase.wdl" as gcnv
import "SingleSampleFiltering.wdl" as SingleSampleFiltering
import "GATKSVPipelineSingleSampleMetrics.wdl" as SingleSampleMetrics
import "Utils.wdl" as utils
import "Structs.wdl"

# GATK SV Pipeline single sample mode
# Runs GatherSampleEvidence, EvidenceQC, GatherBatchEvidence, ClusterBatch, FilterBatch.MergePesrVcfs, GenotypeBatch, 
# MakeCohortVcf (CombineBatches, ResolveComplexVariants, GenotypeComplexVariants, GenotypeComplexVariants), and AnnotateVcf

workflow GATKSVPipelineSingleSample {
  meta {
    allowNestedInputs: true
  }

  input {
    # Batch info
    String batch
    String sample_id

    # Define raw callers to use
    # Overrides presence of case_*_vcf parameters below
    Boolean use_manta = true
    Boolean use_melt = true
    Boolean use_scramble = false
    Boolean use_wham = true

    # If GatherSampleEvidence outputs already prepared
    File? case_manta_vcf
    File? case_melt_vcf
    File? case_scramble_vcf
    File? case_wham_vcf
    File? case_counts_file
    File? case_pe_file
    File? case_sr_file
    File? case_sd_file

    # Global files
    File ref_ped_file
    File ref_samples_list
    File genome_file
    File primary_contigs_list
    File primary_contigs_fai
    File reference_fasta
    File reference_index    # Index (.fai), must be in same dir as fasta
    File reference_dict     # Dictionary (.dict), must be in same dir as fasta
    File ref_panel_vcf
    File autosome_file      # fai of autosomal contigs
    File allosome_file      # fai of allosomal contigs

    String sv_base_mini_docker
    String sv_base_docker
    String sv_pipeline_docker
    String sv_pipeline_hail_docker
    String sv_pipeline_updates_docker
    String sv_pipeline_rdtest_docker
    String sv_pipeline_base_docker
    String sv_pipeline_qc_docker
    String linux_docker
    String cnmops_docker
    String gatk_docker
    String? gcnv_gatk_docker
    String? gatk_docker_pesr_override
    String condense_counts_docker
    String genomes_in_the_cloud_docker
    String samtools_cloud_docker
    String cloud_sdk_docker

    # Must be provided if corresponding use_* is true and case_*_vcf is not provided
    String? manta_docker
    String? melt_docker
    String? scramble_docker
    String? wham_docker

    ############################################################
    ## GatherSampleEvidence
    ############################################################

    # Required if any GatherSampleEvidence outputs need to be generated (vcfs, counts, pe/sr files)
    # (When "If GatherSampleEvidence outputs already prepared" section above is used)
    File? bam_or_cram_file
    File? bam_or_cram_index

    # Common parameters
    String? reference_version   # Either "38" or "19"

    # Coverage collection inputs
    File preprocessed_intervals

    # Manta inputs
    File manta_region_bed
    File? manta_region_bed_index
    Float? manta_jobs_per_cpu
    Int? manta_mem_gb_per_job

    # PESR inputs
    File sd_locs_vcf

    # Melt inputs
    File? melt_standard_vcf_header # required if use_melt True
    File? melt_metrics_intervals
    Float? insert_size
    Int? read_length
    Float? coverage
    File? metrics_intervals
    Int? pf_reads_improper_pairs
    Float? pct_chimeras
    Float? total_reads

    # Wham inputs
    File wham_include_list_bed_file

    # Run GatherSampleEvidence metrics - default is off for single sample pipeline
    Boolean? run_sampleevidence_metrics = false

    # Runtime configuration overrides
    RuntimeAttr? runtime_attr_localize_reads
    RuntimeAttr? runtime_attr_manta
    RuntimeAttr? runtime_attr_melt_coverage
    RuntimeAttr? runtime_attr_melt_metrics
    RuntimeAttr? runtime_attr_melt
    RuntimeAttr? runtime_attr_scramble_part1
    RuntimeAttr? runtime_attr_scramble_part2
    RuntimeAttr? runtime_attr_pesr
    RuntimeAttr? runtime_attr_wham

    ############################################################
    ## EvidenceQC
    ############################################################

    # Optional QC tasks
    Boolean run_vcf_qc

    # WGD files
    File wgd_scoring_mask

    RuntimeAttr? runtime_attr_qc
    RuntimeAttr? runtime_attr_qc_outlier
    RuntimeAttr? ploidy_score_runtime_attr
    RuntimeAttr? ploidy_build_runtime_attr
    RuntimeAttr? wgd_build_runtime_attr
    RuntimeAttr? wgd_score_runtime_attr

    ############################################################
    ## GatherBatchEvidence
    ############################################################

    # Parameters
    Int min_svsize                  # Minimum SV length to include

    # gCNV inputs
    File contig_ploidy_model_tar
    File gcnv_model_tars_list  # list of files, one per line

    # bincov counts files (for cn.mops)
    File ref_panel_bincov_matrix

    File ref_pesr_disc_files_list  # list of files, one per line
    File ref_pesr_split_files_list  # list of files, one per line
    File ref_pesr_sd_files_list  # list of files, one per line

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
    File cnmops_exclude_list

    Int? cnmops_large_min_size

    # QC files
    Int matrix_qc_distance

    # Run GatherBatchEvidence metrics - default is off for single sample pipeline
    Boolean? run_batchevidence_metrics = false

    RuntimeAttr? median_cov_runtime_attr        # Memory ignored, use median_cov_mem_gb_per_sample
    Float? median_cov_mem_gb_per_sample

    RuntimeAttr? evidence_merging_bincov_runtime_attr # Disk space ignored, use evidence_merging_bincov_size_mb

    RuntimeAttr? cnmops_sample10_runtime_attr
    RuntimeAttr? cnmops_sample3_runtime_attr

    RuntimeAttr? add_sample_to_ped_runtime_attr
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

    # Depth merging parameters
    RuntimeAttr? runtime_attr_depth_merge_pre_clusterbatch

    # VCF merging parameters
    RuntimeAttr? runtime_attr_combine_manta_std
    RuntimeAttr? runtime_attr_combine_wham_std
    RuntimeAttr? runtime_attr_combine_scramble_std
    RuntimeAttr? runtime_attr_combine_melt_std

    # Reference panel standardized caller VCFs
    File ref_std_manta_vcf_tar
    File ref_std_wham_vcf_tar
    File? ref_std_scramble_vcf_tar
    File? ref_std_melt_vcf_tar

    File ref_panel_del_bed
    File ref_panel_dup_bed

    Int? depth_records_per_bed_shard_cluster_batch
    Float depth_exclude_overlap_fraction
    Float depth_interval_overlap
    String? depth_clustering_algorithm

    Int? pesr_min_size
    File pesr_exclude_intervals
    Float pesr_interval_overlap
    Int pesr_breakend_window
    String? pesr_clustering_algorithm

    File? baseline_depth_vcf_cluster_batch
    File? baseline_manta_vcf_cluster_batch
    File? baseline_wham_vcf_cluster_batch
    File? baseline_scramble_vcf_cluster_batch
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

    # Run ClusterBatch metrics - default is off for single sample pipeline
    Boolean? run_clusterbatch_metrics = false

    RuntimeAttr? runtime_attr_filter_vcf_by_id

    ############################################################
    ## GenerateBatchMetrics/FilterBatch
    ############################################################

    File rmsk
    File segdups

    Int? min_large_pesr_call_size_for_filtering
    Float? min_large_pesr_depth_overlap_fraction

    RuntimeAttr? runtime_attr_filter_large_pesr
    RuntimeAttr? runtime_attr_srtest
    RuntimeAttr? runtime_attr_split_vcf_srtest
    RuntimeAttr? runtime_attr_merge_allo
    RuntimeAttr? runtime_attr_merge_stats
    RuntimeAttr? runtime_attr_rewritesrcoords

    RuntimeAttr? runtime_attr_merge_pesr_vcfs
    RuntimeAttr? runtime_attr_get_male_only

    ############################################################
    ## GenotypeBatch
    ############################################################

    Int genotyping_n_per_split
    Int n_RD_genotype_bins  # number of RdTest bins

    File cutoffs

    File genotype_pesr_pesr_sepcutoff
    File genotype_pesr_depth_sepcutoff
    File genotype_depth_pesr_sepcutoff
    File genotype_depth_depth_sepcutoff

    File SR_metrics
    File PE_metrics

    File bin_exclude

    # Run GenotypeBatch metrics - default is off for single sample pipeline
    Boolean? run_genotypebatch_metrics = false

    # Common
    RuntimeAttr? runtime_attr_merge_counts
    RuntimeAttr? runtime_attr_split_variants
    RuntimeAttr? runtime_attr_make_subset_vcf
    RuntimeAttr? runtime_attr_rdtest_genotype
    RuntimeAttr? runtime_attr_add_genotypes
    RuntimeAttr? runtime_attr_genotype_depths_concat_vcfs
    RuntimeAttr? runtime_attr_genotype_pesr_concat_vcfs
    RuntimeAttr? runtime_attr_split_vcf_genotypebatch


    # Master
    RuntimeAttr? runtime_attr_add_batch
    RuntimeAttr? runtime_attr_index_vcf

    # PESR part 2
    RuntimeAttr? runtime_attr_count_pe
    RuntimeAttr? runtime_attr_genotype_pe
    RuntimeAttr? runtime_attr_count_sr
    RuntimeAttr? runtime_attr_genotype_sr
    RuntimeAttr? runtime_attr_integrate_gq
    RuntimeAttr? runtime_attr_integrate_pesr_gq
    RuntimeAttr? runtime_attr_triple_stream_cat

    # Depth part 2
    RuntimeAttr? runtime_attr_integrate_depth_gq

    ############################################################
    ## MakeCohortVcf
    ############################################################

    Float clean_vcf_min_sr_background_fail_batches

    File cytobands

    File mei_bed
    File depth_exclude_list
    File empty_file

    Int max_shard_size_resolve
    Int clean_vcf_max_shards_per_chrom_clean_vcf_step1
    Int clean_vcf_min_records_per_shard_clean_vcf_step1
    Int clean_vcf_samples_per_clean_vcf_step2_shard
    Int clean_vcf5_records_per_shard
    Int clean_vcf1b_records_per_shard

    String? chr_x
    String? chr_y

    Int? clean_vcf_random_seed

    # Run MakeCohortVcf metrics - default is off for single sample pipeline
    Boolean? run_makecohortvcf_metrics = false

    # overrides for local tasks
    RuntimeAttr? runtime_overide_get_discfile_size
    RuntimeAttr? runtime_override_update_sr_list_cluster
    RuntimeAttr? runtime_override_merge_pesr_depth
    RuntimeAttr? runtime_override_integrate_resolved_vcfs
    RuntimeAttr? runtime_override_rename_variants
    RuntimeAttr? runtime_override_rename_cleaned_samples

    # overrides for mini tasks
    RuntimeAttr? runtime_override_clean_background_fail
    RuntimeAttr? runtime_override_make_cpx_cnv_input_file
    RuntimeAttr? runtime_override_subset_inversions
    RuntimeAttr? runtime_override_concat_merged_vcfs
    RuntimeAttr? runtime_override_concat_cpx_vcfs
    RuntimeAttr? runtime_override_concat_cleaned_vcfs

    # overrides for VcfClusterContig
    RuntimeAttr? runtime_override_join_vcfs
    RuntimeAttr? runtime_override_subset_bothside_pass
    RuntimeAttr? runtime_override_subset_background_fail
    RuntimeAttr? runtime_override_subset_sv_type
    RuntimeAttr? runtime_override_shard_clusters
    RuntimeAttr? runtime_override_shard_vids
    RuntimeAttr? runtime_override_pull_vcf_shard
    RuntimeAttr? runtime_override_svtk_vcf_cluster
    RuntimeAttr? runtime_override_get_vcf_header_with_members_info_line
    RuntimeAttr? runtime_override_cluster_merge
    RuntimeAttr? runtime_override_concat_vcf_cluster
    RuntimeAttr? runtime_override_concat_svtypes
    RuntimeAttr? runtime_override_concat_sharded_cluster
    RuntimeAttr? runtime_override_make_sites_only
    RuntimeAttr? runtime_override_preconcat_sharded_cluster
    RuntimeAttr? runtime_override_hail_merge_sharded_cluster
    RuntimeAttr? runtime_override_fix_header_sharded_cluster
    RuntimeAttr? runtime_override_concat_large_pesr_depth

    # overrides for ResolveComplexVariants
    RuntimeAttr? runtime_override_update_sr_list_pass
    RuntimeAttr? runtime_override_update_sr_list_fail
    RuntimeAttr? runtime_override_breakpoint_overlap_filter
    RuntimeAttr? runtime_override_concat_resolve

    RuntimeAttr? runtime_override_get_se_cutoff
    RuntimeAttr? runtime_override_shard_vcf_cpx
    RuntimeAttr? runtime_override_shard_vids_resolve
    RuntimeAttr? runtime_override_resolve_prep
    RuntimeAttr? runtime_override_resolve_cpx_per_shard
    RuntimeAttr? runtime_override_restore_unresolved_cnv_per_shard
    RuntimeAttr? runtime_override_concat_resolved_per_shard
    RuntimeAttr? runtime_override_preconcat_resolve
    RuntimeAttr? runtime_override_hail_merge_resolve
    RuntimeAttr? runtime_override_fix_header_resolve

    RuntimeAttr? runtime_override_get_se_cutoff_inv
    RuntimeAttr? runtime_override_shard_vcf_cpx_inv
    RuntimeAttr? runtime_override_shard_vids_resolve_inv
    RuntimeAttr? runtime_override_resolve_prep_inv
    RuntimeAttr? runtime_override_resolve_cpx_per_shard_inv
    RuntimeAttr? runtime_override_restore_unresolved_cnv_per_shard_inv
    RuntimeAttr? runtime_override_concat_resolved_per_shard_inv
    RuntimeAttr? runtime_override_pull_vcf_shard_inv
    RuntimeAttr? runtime_override_preconcat_resolve_inv
    RuntimeAttr? runtime_override_hail_merge_resolve_inv
    RuntimeAttr? runtime_override_fix_header_resolve_inv

    # overrides for GenotypeComplexContig
    RuntimeAttr? runtime_override_ids_from_median
    RuntimeAttr? runtime_override_split_vcf_to_genotype
    RuntimeAttr? runtime_override_concat_cpx_cnv_vcfs
    RuntimeAttr? runtime_override_get_cpx_cnv_intervals
    RuntimeAttr? runtime_override_parse_genotypes
    RuntimeAttr? runtime_override_merge_melted_gts
    RuntimeAttr? runtime_override_split_bed_by_size
    RuntimeAttr? runtime_override_rd_genotype
    RuntimeAttr? runtime_override_concat_melted_genotypes
    RuntimeAttr? runtime_attr_ids_from_vcf_regeno
    RuntimeAttr? runtime_attr_subset_ped_regeno
    RuntimeAttr? runtime_override_preconcat_regeno
    RuntimeAttr? runtime_override_hail_merge_regeno
    RuntimeAttr? runtime_override_fix_header_regeno

    # overrides for CleanVcfContig
    RuntimeAttr? runtime_override_preconcat_clean_final
    RuntimeAttr? runtime_override_hail_merge_clean_final
    RuntimeAttr? runtime_override_fix_header_clean_final
    RuntimeAttr? runtime_attr_format_clean

    RuntimeAttr? runtime_override_clean_vcf_1a
    RuntimeAttr? runtime_override_clean_vcf_2
    RuntimeAttr? runtime_override_clean_vcf_3
    RuntimeAttr? runtime_override_clean_vcf_4
    RuntimeAttr? runtime_override_clean_vcf_5_scatter
    RuntimeAttr? runtime_override_clean_vcf_5_make_cleangq
    RuntimeAttr? runtime_override_clean_vcf_5_find_redundant_multiallelics
    RuntimeAttr? runtime_override_clean_vcf_5_polish
    RuntimeAttr? runtime_override_stitch_fragmented_cnvs
    RuntimeAttr? runtime_override_final_cleanup

    RuntimeAttr? runtime_attr_override_subset_large_cnvs_1b
    RuntimeAttr? runtime_attr_override_sort_bed_1b
    RuntimeAttr? runtime_attr_override_intersect_bed_1b
    RuntimeAttr? runtime_attr_override_build_dict_1b
    RuntimeAttr? runtime_attr_override_scatter_1b
    RuntimeAttr? runtime_attr_override_filter_vcf_1b
    RuntimeAttr? runtime_override_concat_vcfs_1b
    RuntimeAttr? runtime_override_cat_multi_cnvs_1b

    RuntimeAttr? runtime_override_preconcat_step1
    RuntimeAttr? runtime_override_hail_merge_step1
    RuntimeAttr? runtime_override_fix_header_step1

    RuntimeAttr? runtime_override_preconcat_drc
    RuntimeAttr? runtime_override_hail_merge_drc
    RuntimeAttr? runtime_override_fix_header_drc

    RuntimeAttr? runtime_override_split_vcf_to_clean
    RuntimeAttr? runtime_override_combine_step_1_sex_chr_revisions
    RuntimeAttr? runtime_override_split_include_list
    RuntimeAttr? runtime_override_combine_clean_vcf_2
    RuntimeAttr? runtime_override_combine_revised_4
    RuntimeAttr? runtime_override_combine_multi_ids_4
    RuntimeAttr? runtime_override_drop_redundant_cnvs
    RuntimeAttr? runtime_override_combine_step_1_vcfs
    RuntimeAttr? runtime_override_sort_drop_redundant_cnvs

    # overrides for VcfQc
    RuntimeAttr? runtime_override_site_level_benchmark_plot
    RuntimeAttr? runtime_override_per_sample_benchmark_plot
    RuntimeAttr? runtime_override_subset_vcf
    RuntimeAttr? runtime_override_preprocess_vcf
    RuntimeAttr? runtime_override_site_level_benchmark
    RuntimeAttr? runtime_override_merge_site_level_benchmark
    RuntimeAttr? runtime_override_merge_sharded_per_sample_vid_lists
    RuntimeAttr? runtime_override_plot_qc_vcf_wide
    RuntimeAttr? runtime_override_plot_qc_per_sample
    RuntimeAttr? runtime_override_plot_qc_per_family
    RuntimeAttr? runtime_override_sanitize_outputs
    RuntimeAttr? runtime_override_merge_vcfwide_stat_shards
    RuntimeAttr? runtime_override_merge_vcf_2_bed
    RuntimeAttr? runtime_override_collect_sharded_vcf_stats
    RuntimeAttr? runtime_override_svtk_vcf_2_bed
    RuntimeAttr? runtime_override_scatter_vcf
    RuntimeAttr? runtime_override_merge_subvcf_stat_shards
    RuntimeAttr? runtime_override_collect_vids_per_sample
    RuntimeAttr? runtime_override_split_samples_list
    RuntimeAttr? runtime_override_tar_shard_vid_lists
    RuntimeAttr? runtime_override_benchmark_samples
    RuntimeAttr? runtime_override_split_shuffled_list
    RuntimeAttr? runtime_override_merge_and_tar_shard_benchmarks

    ############################################################
    ## AnnotateVcf
    ############################################################

    File protein_coding_gtf
    File noncoding_bed
    Int? promoter_window
    Int? max_breakend_as_cnv_length
    Int annotation_sv_per_shard

    File? external_af_ref_bed             # bed file with population AFs for annotation
    String? external_af_ref_bed_prefix    # name of external AF bed file call set
    Array[String]? external_af_population # populations to annotate external AFs (required if ref_bed set, use "ALL" for all)

    RuntimeAttr? runtime_attr_svannotate

    ############################################################
    ## Single sample filtering
    ############################################################

    Float? max_ref_panel_carrier_freq

    ############################################################
    ## Single sample metrics
    ############################################################

    File? baseline_cleaned_vcf
    File? baseline_final_vcf
    File? baseline_genotyped_pesr_vcf
    File? baseline_genotyped_depth_vcf
    File? baseline_non_genotyped_unique_depth_calls_vcf

    ############################################################
    ## QC
    ############################################################

    File qc_definitions

    # Do not use
    String? NONE_STRING_

  }

  String? manta_docker_ = if (!defined(case_manta_vcf) && use_manta) then manta_docker else NONE_STRING_
  String? melt_docker_ = if (!defined(case_melt_vcf) && use_melt) then melt_docker else NONE_STRING_
  String? scramble_docker_ = if (!defined(case_scramble_vcf) && use_scramble) then scramble_docker else NONE_STRING_
  String? wham_docker_ = if (!defined(case_wham_vcf) && use_wham) then wham_docker else NONE_STRING_

  Boolean collect_coverage = !defined(case_counts_file)
  Boolean collect_pesr = !defined(case_pe_file) || !defined(case_sr_file) || !defined(case_sd_file)

  Boolean run_sampleevidence = defined(manta_docker_) || defined(melt_docker_) || defined(scramble_docker_) || defined(wham_docker_) || collect_coverage || collect_pesr

  if (run_sampleevidence) {
    call sampleevidence.GatherSampleEvidence as GatherSampleEvidence {
      input:
        bam_or_cram_file=select_first([bam_or_cram_file]),
        bam_or_cram_index=bam_or_cram_index,
        sample_id=sample_id,
        collect_coverage = collect_coverage,
        collect_pesr = collect_pesr,
        primary_contigs_list=primary_contigs_list,
        reference_fasta=reference_fasta,
        reference_index=reference_index,
        reference_dict=reference_dict,
        reference_version=reference_version,
        preprocessed_intervals=preprocessed_intervals,
        manta_region_bed=manta_region_bed,
        manta_region_bed_index=manta_region_bed_index,
        manta_jobs_per_cpu=manta_jobs_per_cpu,
        manta_mem_gb_per_job=manta_mem_gb_per_job,
        sd_locs_vcf=sd_locs_vcf,
        melt_standard_vcf_header=melt_standard_vcf_header,
        melt_metrics_intervals=melt_metrics_intervals,
        insert_size=insert_size,
        read_length=read_length,
        coverage=coverage,
        metrics_intervals=metrics_intervals,
        pf_reads_improper_pairs=pf_reads_improper_pairs,
        pct_chimeras=pct_chimeras,
        total_reads=total_reads,
        wham_include_list_bed_file=wham_include_list_bed_file,
        run_module_metrics = run_sampleevidence_metrics,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        manta_docker=manta_docker_,
        melt_docker=melt_docker_,
        scramble_docker=scramble_docker_,
        wham_docker=wham_docker_,
        gatk_docker=gatk_docker,
        gatk_docker_pesr_override = gatk_docker_pesr_override,
        genomes_in_the_cloud_docker=genomes_in_the_cloud_docker,
        samtools_cloud_docker=samtools_cloud_docker,
        cloud_sdk_docker = cloud_sdk_docker,
        runtime_attr_localize_reads=runtime_attr_localize_reads,
        runtime_attr_manta=runtime_attr_manta,
        runtime_attr_melt_coverage=runtime_attr_melt_coverage,
        runtime_attr_melt_metrics=runtime_attr_melt_metrics,
        runtime_attr_melt=runtime_attr_melt,
        runtime_attr_scramble_part1=runtime_attr_scramble_part1,
        runtime_attr_scramble_part2=runtime_attr_scramble_part2,
        runtime_attr_pesr=runtime_attr_pesr,
        runtime_attr_wham=runtime_attr_wham
    }
  }

  File case_counts_file_ = select_first([case_counts_file, GatherSampleEvidence.coverage_counts])
  File case_pe_file_ = select_first([case_pe_file, GatherSampleEvidence.pesr_disc])
  File case_sr_file_ = select_first([case_sr_file, GatherSampleEvidence.pesr_split])
  File case_sd_file_ = select_first([case_sd_file, GatherSampleEvidence.pesr_sd])

  call evidenceqc.EvidenceQC as EvidenceQC {
    input:
      batch=batch,
      samples=[sample_id],
      run_vcf_qc=run_vcf_qc,
      genome_file=genome_file,
      counts=[case_counts_file_],
      run_ploidy = false,
      wgd_scoring_mask=wgd_scoring_mask,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_qc_docker=sv_pipeline_qc_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_base_docker=sv_base_docker,
      runtime_attr_qc=runtime_attr_qc,
      runtime_attr_qc_outlier=runtime_attr_qc_outlier,
      wgd_build_runtime_attr=wgd_build_runtime_attr,
      wgd_score_runtime_attr=wgd_score_runtime_attr
  }

  if (use_manta) {
    Array[File] manta_vcfs_ = [select_first([case_manta_vcf, GatherSampleEvidence.manta_vcf])]
  }
  if (use_melt) {
    Array[File] melt_vcfs_ = [select_first([case_melt_vcf, GatherSampleEvidence.melt_vcf])]
  }
  if (use_scramble) {
    Array[File] scramble_vcfs_ = [select_first([case_scramble_vcf, GatherSampleEvidence.scramble_vcf])]
  }
  if (use_wham) {
    Array[File] wham_vcfs_ = [select_first([case_wham_vcf, GatherSampleEvidence.wham_vcf])]
  }

  Array[String] ref_samples = read_lines(ref_samples_list)

  call batchevidence.GatherBatchEvidence as GatherBatchEvidence {
    input:
      batch=batch,
      samples=[sample_id],
      ref_panel_samples=ref_samples,
      run_matrix_qc=false,
      ped_file=ref_ped_file,
      genome_file=genome_file,
      primary_contigs_fai=primary_contigs_fai,
      ref_dict=reference_dict,
      counts=[case_counts_file_],
      ref_panel_bincov_matrix=ref_panel_bincov_matrix,
      bincov_matrix=EvidenceQC.bincov_matrix,
      bincov_matrix_index=EvidenceQC.bincov_matrix_index,
      PE_files=[case_pe_file_],
      cytoband=cytobands,
      mei_bed=mei_bed,
      ref_panel_PE_files=read_lines(ref_pesr_disc_files_list),
      SR_files=[case_sr_file_],
      ref_panel_SR_files=read_lines(ref_pesr_split_files_list),
      SD_files=[case_sd_file_],
      ref_panel_SD_files=read_lines(ref_pesr_sd_files_list),
      sd_locs_vcf=sd_locs_vcf,
      contig_ploidy_model_tar = contig_ploidy_model_tar,
      gcnv_model_tars = read_lines(gcnv_model_tars_list),
      gatk4_jar_override = gatk4_jar_override,
      run_ploidy = true,
      append_first_sample_to_ped = true,
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
      manta_vcfs=manta_vcfs_,
      melt_vcfs=melt_vcfs_,
      scramble_vcfs=scramble_vcfs_,
      wham_vcfs=wham_vcfs_,
      min_svsize=min_svsize,
      cnmops_chrom_file=autosome_file,
      cnmops_exclude_list=cnmops_exclude_list,
      cnmops_allo_file=allosome_file,
      cnmops_large_min_size=cnmops_large_min_size,
      matrix_qc_distance=matrix_qc_distance,
      run_module_metrics = run_batchevidence_metrics,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_base_docker=sv_base_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_qc_docker=sv_pipeline_qc_docker,
      linux_docker=linux_docker,
      cnmops_docker=cnmops_docker,
      gatk_docker = gatk_docker,
      gcnv_gatk_docker=gcnv_gatk_docker,
      condense_counts_docker = condense_counts_docker,
      median_cov_runtime_attr=median_cov_runtime_attr,
      median_cov_mem_gb_per_sample=median_cov_mem_gb_per_sample,
      evidence_merging_bincov_runtime_attr=evidence_merging_bincov_runtime_attr,
      cnmops_sample10_runtime_attr=cnmops_sample10_runtime_attr,
      cnmops_sample3_runtime_attr=cnmops_sample3_runtime_attr,
      preprocess_calls_runtime_attr=preprocess_calls_runtime_attr,
      depth_merge_set_runtime_attr=depth_merge_set_runtime_attr,
      depth_merge_sample_runtime_attr=depth_merge_sample_runtime_attr,
      cnmops_ped_runtime_attr=cnmops_ped_runtime_attr,
      cnmops_clean_runtime_attr=cnmops_clean_runtime_attr,
      matrix_qc_pesrbaf_runtime_attr=matrix_qc_pesrbaf_runtime_attr,
      matrix_qc_rd_runtime_attr=matrix_qc_rd_runtime_attr,
      ploidy_score_runtime_attr=ploidy_score_runtime_attr,
      ploidy_build_runtime_attr=ploidy_build_runtime_attr,
      add_sample_to_ped_runtime_attr=add_sample_to_ped_runtime_attr,
      runtime_attr_ploidy = runtime_attr_ploidy,
      runtime_attr_case = runtime_attr_case,
      runtime_attr_postprocess = runtime_attr_postprocess,
      runtime_attr_explode = runtime_attr_explode
  }

  File combined_ped_file = select_first([GatherBatchEvidence.combined_ped_file])

  # Merge calls with reference panel
  if (defined(GatherBatchEvidence.std_manta_vcf_tar)) {
    call utils.CombineTars as CombineMantaStd {
      input:
        tar1=ref_std_manta_vcf_tar,
        tar2=select_first([GatherBatchEvidence.std_manta_vcf_tar]),
        linux_docker=linux_docker,
        runtime_attr_override=runtime_attr_combine_manta_std
    }
  }
  if (defined(GatherBatchEvidence.std_wham_vcf_tar)) {
    call utils.CombineTars as CombineWhamStd {
      input:
        tar1=ref_std_wham_vcf_tar,
        tar2=select_first([GatherBatchEvidence.std_wham_vcf_tar]),
        linux_docker=linux_docker,
        runtime_attr_override=runtime_attr_combine_wham_std
    }
  }
  if (defined(GatherBatchEvidence.std_scramble_vcf_tar) && defined(ref_std_scramble_vcf_tar)) {
    call utils.CombineTars as CombineScrambleStd {
      input:
        tar1=select_first([ref_std_scramble_vcf_tar]),
        tar2=select_first([GatherBatchEvidence.std_scramble_vcf_tar]),
        linux_docker=linux_docker,
        runtime_attr_override=runtime_attr_combine_scramble_std
    }
  }
  if (defined(GatherBatchEvidence.std_melt_vcf_tar) && defined(ref_std_melt_vcf_tar)) {
    call utils.CombineTars as CombineMeltStd {
      input:
        tar1=select_first([ref_std_melt_vcf_tar]),
        tar2=select_first([GatherBatchEvidence.std_melt_vcf_tar]),
        linux_docker=linux_docker,
        runtime_attr_override=runtime_attr_combine_melt_std
    }
  }
  File merged_manta_vcf_tar = select_first([CombineMantaStd.out, ref_std_manta_vcf_tar])
  File merged_wham_vcf_tar = select_first([CombineWhamStd.out, ref_std_wham_vcf_tar])
  if (defined(CombineScrambleStd.out) || defined(ref_std_scramble_vcf_tar)) {
    File merged_scramble_vcf_tar = select_first([CombineScrambleStd.out, ref_std_scramble_vcf_tar])
  }
  if (defined(CombineMeltStd.out) || defined(ref_std_melt_vcf_tar)) {
    File merged_melt_vcf_tar = select_first([CombineMeltStd.out, ref_std_melt_vcf_tar])
  }

  call dpn.MergeSet as MergeSetDel {
    input:
      beds=[GatherBatchEvidence.merged_dels, ref_panel_del_bed],
      svtype="DEL",
      batch=batch,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_depth_merge_pre_clusterbatch
  }
  call dpn.MergeSet as MergeSetDup {
    input:
      beds=[GatherBatchEvidence.merged_dups, ref_panel_dup_bed],
      svtype="DUP",
      batch=batch,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_depth_merge_pre_clusterbatch
  }

  call clusterbatch.ClusterBatch {
    input:
      manta_vcf_tar=merged_manta_vcf_tar,
      wham_vcf_tar=merged_wham_vcf_tar,
      scramble_vcf_tar=merged_scramble_vcf_tar,
      melt_vcf_tar=merged_melt_vcf_tar,
      del_bed=MergeSetDel.out,
      dup_bed=MergeSetDup.out,
      batch=batch,
      ped_file=combined_ped_file,
      contig_list=primary_contigs_list,
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_index,
      reference_dict=reference_dict,
      chr_x=chr_x,
      chr_y=chr_y,
      depth_records_per_bed_shard=depth_records_per_bed_shard_cluster_batch,
      depth_exclude_intervals=depth_exclude_list,
      depth_exclude_overlap_fraction=depth_exclude_overlap_fraction,
      depth_interval_overlap=depth_interval_overlap,
      depth_clustering_algorithm=depth_clustering_algorithm,
      pesr_min_size=pesr_min_size,
      pesr_exclude_intervals=pesr_exclude_intervals,
      pesr_interval_overlap=pesr_interval_overlap,
      pesr_breakend_window=pesr_breakend_window,
      pesr_clustering_algorithm=pesr_clustering_algorithm,
      run_module_metrics=run_clusterbatch_metrics,
      linux_docker=linux_docker,
      sv_pipeline_base_docker=sv_pipeline_base_docker,
      baseline_depth_vcf=baseline_depth_vcf_cluster_batch,
      baseline_manta_vcf=baseline_manta_vcf_cluster_batch,
      baseline_wham_vcf=baseline_wham_vcf_cluster_batch,
      baseline_scramble_vcf=baseline_scramble_vcf_cluster_batch,
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
      runtime_attr_svcluster_wham=runtime_attr_svcluster_wham_cluster_batch,
      runtime_attr_svcluster_scramble=runtime_attr_svcluster_scramble_cluster_batch,
      runtime_override_concat_vcfs_pesr=runtime_override_concat_vcfs_pesr_cluster_batch,
      runtime_attr_gatk_to_svtk_vcf_pesr=runtime_attr_gatk_to_svtk_vcf_pesr_cluster_batch,
      runtime_attr_scatter_bed=runtime_attr_scatter_bed_cluster_batch,
      runtime_attr_cnv_bed_to_gatk_vcf=runtime_attr_cnv_bed_to_gatk_vcf_cluster_batch,
      runtime_attr_exclude_intervals_depth=runtime_attr_exclude_intervals_depth_cluster_batch,
      runtime_attr_svcluster_depth=runtime_attr_svcluster_depth_cluster_batch,
      runtime_attr_gatk_to_svtk_vcf_depth=runtime_attr_gatk_to_svtk_vcf_depth_cluster_batch,
      runtime_override_concat_vcfs_depth=runtime_override_concat_vcfs_depth_cluster_batch,
      runtime_attr_exclude_intervals_pesr=runtime_attr_exclude_intervals_pesr_cluster_batch
  }

  # Pull out clustered calls from this sample only
  if (use_manta) {
    call SingleSampleFiltering.FilterVcfBySampleGenotypeAndAddEvidenceAnnotation as FilterManta {
        input :
            vcf_gz=select_first([ClusterBatch.clustered_manta_vcf]),
            sample_id=sample_id,
            evidence="RD,PE,SR",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_filter_vcf_by_id
    }
  }
  if (use_wham) {
    call SingleSampleFiltering.FilterVcfBySampleGenotypeAndAddEvidenceAnnotation as FilterWham {
        input :
            vcf_gz=select_first([ClusterBatch.clustered_wham_vcf]),
            sample_id=sample_id,
            evidence="RD,PE,SR",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_filter_vcf_by_id
    }
  }
  if (use_melt) {
    call SingleSampleFiltering.FilterVcfBySampleGenotypeAndAddEvidenceAnnotation as FilterMelt {
        input :
            vcf_gz=select_first([ClusterBatch.clustered_melt_vcf]),
            sample_id=sample_id,
            evidence="RD,PE,SR",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_filter_vcf_by_id
    }
  }
  if (use_scramble) {
    call SingleSampleFiltering.FilterVcfBySampleGenotypeAndAddEvidenceAnnotation as FilterScramble {
        input :
            vcf_gz=select_first([ClusterBatch.clustered_scramble_vcf]),
            sample_id=sample_id,
            evidence="RD,PE,SR",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_filter_vcf_by_id
    }
  }

  call SingleSampleFiltering.FilterVcfBySampleGenotypeAndAddEvidenceAnnotation as FilterDepth {
    input :
      vcf_gz=ClusterBatch.clustered_depth_vcf,
      sample_id=sample_id,
      evidence="RD",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_filter_vcf_by_id
  }

  call tasks_makecohortvcf.ConcatVcfs as MergePesrVcfs {
    input:
      vcfs=select_all([FilterManta.out, FilterWham.out, FilterScramble.out, FilterMelt.out]),
      vcfs_idx=select_all([FilterManta.out_index, FilterWham.out_index, FilterScramble.out_index, FilterMelt.out_index]),
      allow_overlaps=true,
      outfile_prefix="~{batch}.filtered_pesr_merged",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_merge_pesr_vcfs
  }

  call SingleSampleFiltering.FilterLargePESRCallsWithoutRawDepthSupport as FilterLargePESRCallsWithoutRawDepthSupport {
    input:
      pesr_vcf=MergePesrVcfs.concat_vcf,
      raw_dels=GatherBatchEvidence.merged_dels,
      raw_dups=GatherBatchEvidence.merged_dups,
      min_large_pesr_call_size_for_filtering=min_large_pesr_call_size_for_filtering,
      min_large_pesr_depth_overlap_fraction=min_large_pesr_depth_overlap_fraction,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override=runtime_attr_filter_large_pesr
  }

  call batchmetrics.GetSampleLists as SamplesList {
    input:
      ped_file = combined_ped_file,
      samples = flatten([[sample_id], ref_samples]),
      sv_base_docker = sv_base_docker
  }

  call batchmetrics.GetMaleOnlyVariantIDs {
    input:
      vcf = ClusterBatch.clustered_depth_vcf,
      female_samples = SamplesList.female_samples,
      male_samples = SamplesList.male_samples,
      contig = select_first([chr_x, "chrX"]),
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_get_male_only
  }

  call SRTest.SRTest as SRTest {
    input:
      splitfile = GatherBatchEvidence.merged_SR,
      medianfile = GatherBatchEvidence.median_cov,
      ped_file = combined_ped_file,
      vcf = FilterLargePESRCallsWithoutRawDepthSupport.out,
      autosome_contigs = autosome_file,
      ref_dict=reference_dict,
      split_size = genotyping_n_per_split,
      algorithm = "PESR",
      allosome_contigs = allosome_file,
      batch = batch,
      samples = SamplesList.samples_file,
      male_samples = SamplesList.male_samples,
      female_samples = SamplesList.female_samples,
      male_only_variant_ids = GetMaleOnlyVariantIDs.male_only_variant_ids,
      run_common = false,
      sv_base_mini_docker = sv_base_mini_docker,
      linux_docker = linux_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_srtest = runtime_attr_srtest,
      runtime_attr_split_vcf = runtime_attr_split_vcf_srtest,
      runtime_attr_merge_allo = runtime_attr_merge_allo,
      runtime_attr_merge_stats = runtime_attr_merge_stats
  }

  call batchmetrics.AggregateTests as AggregateTests {
    input:
      vcf=FilterLargePESRCallsWithoutRawDepthSupport.out,
      srtest=SRTest.srtest,
      rmsk=rmsk,
      segdups=segdups,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call SingleSampleFiltering.RewriteSRCoords as RewriteSRCoords {
    input:
      vcf = FilterLargePESRCallsWithoutRawDepthSupport.out,
      metrics = AggregateTests.metrics,
      cutoffs = cutoffs,
      prefix = batch,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_rewritesrcoords
  }

  call genotypebatch.GenotypeBatch as GenotypeBatch {
    input:
      batch_pesr_vcf=RewriteSRCoords.annotated_vcf,
      batch_depth_vcf=FilterDepth.out,
      cohort_pesr_vcf=RewriteSRCoords.annotated_vcf,
      cohort_depth_vcf=FilterDepth.out,
      batch=batch,
      n_per_split=genotyping_n_per_split,
      medianfile=GatherBatchEvidence.median_cov,
      coveragefile=GatherBatchEvidence.merged_bincov,
      coveragefile_index=GatherBatchEvidence.merged_bincov_index,
      discfile=GatherBatchEvidence.merged_PE,
      discfile_index=GatherBatchEvidence.merged_PE_index,
      splitfile=GatherBatchEvidence.merged_SR,
      splitfile_index=GatherBatchEvidence.merged_SR_index,
      ref_dict=reference_dict,
      n_RD_genotype_bins=n_RD_genotype_bins,
      genotype_pesr_pesr_sepcutoff=genotype_pesr_pesr_sepcutoff,
      genotype_pesr_depth_sepcutoff=genotype_pesr_depth_sepcutoff,
      genotype_depth_pesr_sepcutoff=genotype_depth_pesr_sepcutoff,
      genotype_depth_depth_sepcutoff=genotype_depth_depth_sepcutoff,
      SR_metrics=SR_metrics,
      PE_metrics=PE_metrics,
      bin_exclude=bin_exclude,
      run_module_metrics = run_genotypebatch_metrics,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
      linux_docker=linux_docker,
      runtime_attr_split_vcf=runtime_attr_split_vcf_genotypebatch,
      runtime_attr_merge_counts=runtime_attr_merge_counts,
      runtime_attr_split_variants=runtime_attr_split_variants,
      runtime_attr_make_subset_vcf=runtime_attr_make_subset_vcf,
      runtime_attr_rdtest_genotype=runtime_attr_rdtest_genotype,
      runtime_attr_add_genotypes=runtime_attr_add_genotypes,
      runtime_attr_genotype_depths_concat_vcfs=runtime_attr_genotype_depths_concat_vcfs,
      runtime_attr_genotype_pesr_concat_vcfs=runtime_attr_genotype_pesr_concat_vcfs,
      runtime_attr_add_batch=runtime_attr_add_batch,
      runtime_attr_index_vcf=runtime_attr_index_vcf,
      runtime_attr_count_pe=runtime_attr_count_pe,
      runtime_attr_genotype_pe=runtime_attr_genotype_pe,
      runtime_attr_count_sr=runtime_attr_count_sr,
      runtime_attr_genotype_sr=runtime_attr_genotype_sr,
      runtime_attr_integrate_gq=runtime_attr_integrate_gq,
      runtime_attr_integrate_pesr_gq=runtime_attr_integrate_pesr_gq,
      runtime_attr_triple_stream_cat=runtime_attr_triple_stream_cat,
      runtime_attr_integrate_depth_gq=runtime_attr_integrate_depth_gq
  }

  call SingleSampleFiltering.ConvertCNVsWithoutDepthSupportToBNDs as ConvertCNVsWithoutDepthSupportToBNDs {
    input:
      genotyped_pesr_vcf=GenotypeBatch.genotyped_pesr_vcf,
      allosome_file=allosome_file,
      merged_famfile=combined_ped_file,
      case_sample=sample_id,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call makecohortvcf.MakeCohortVcf as MakeCohortVcf {
    input:
      raw_sr_bothside_pass_files=[GenotypeBatch.sr_bothside_pass],
      raw_sr_background_fail_files=[GenotypeBatch.sr_background_fail],
      min_sr_background_fail_batches=clean_vcf_min_sr_background_fail_batches,
      ped_file=combined_ped_file,
      pesr_vcfs=[ConvertCNVsWithoutDepthSupportToBNDs.out_vcf],
      depth_vcfs=[GenotypeBatch.genotyped_depth_vcf],
      contig_list=primary_contigs_fai,
      allosome_fai=allosome_file,
      ref_dict=reference_dict,

      merge_complex_genotype_vcfs = true,
      cytobands=cytobands,

      bin_exclude=bin_exclude,

      disc_files=[GatherBatchEvidence.merged_PE],
      bincov_files=[GatherBatchEvidence.merged_bincov],

      mei_bed=mei_bed,
      pe_exclude_list=pesr_exclude_intervals,
      depth_exclude_list=depth_exclude_list,
      empty_file=empty_file,

      cohort_name=batch,

      rf_cutoff_files=[cutoffs],
      batches=[batch],
      depth_gt_rd_sep_files=[genotype_depth_depth_sepcutoff],
      median_coverage_files=[GatherBatchEvidence.median_cov],

      max_shard_size_resolve=max_shard_size_resolve,
      max_shards_per_chrom_clean_vcf_step1=clean_vcf_max_shards_per_chrom_clean_vcf_step1,
      min_records_per_shard_clean_vcf_step1=clean_vcf_min_records_per_shard_clean_vcf_step1,
      samples_per_clean_vcf_step2_shard=clean_vcf_samples_per_clean_vcf_step2_shard,
      clean_vcf5_records_per_shard=clean_vcf5_records_per_shard,
      clean_vcf1b_records_per_shard=clean_vcf1b_records_per_shard,

      chr_x=select_first([chr_x, "chrX"]),
      chr_y=select_first([chr_y, "chrY"]),

      random_seed=clean_vcf_random_seed,

      run_module_metrics = run_makecohortvcf_metrics,

      primary_contigs_list=primary_contigs_list,
      sv_pipeline_base_docker=sv_pipeline_base_docker,

      linux_docker=linux_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_hail_docker=sv_pipeline_hail_docker,
      sv_pipeline_updates_docker=sv_pipeline_updates_docker,
      sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
      sv_pipeline_qc_docker=sv_pipeline_qc_docker,
      sv_base_mini_docker=sv_base_mini_docker,

      runtime_overide_get_discfile_size=runtime_overide_get_discfile_size,
      runtime_override_update_sr_list_cluster=runtime_override_update_sr_list_cluster,
      runtime_override_merge_pesr_depth=runtime_override_merge_pesr_depth,
      runtime_override_integrate_resolved_vcfs=runtime_override_integrate_resolved_vcfs,
      runtime_override_rename_variants=runtime_override_rename_variants,
      runtime_override_rename_cleaned_samples=runtime_override_rename_cleaned_samples,
      runtime_override_breakpoint_overlap_filter=runtime_override_breakpoint_overlap_filter,
      runtime_override_clean_background_fail=runtime_override_clean_background_fail,
      runtime_override_make_cpx_cnv_input_file=runtime_override_make_cpx_cnv_input_file,
      runtime_override_concat_merged_vcfs=runtime_override_concat_merged_vcfs,
      runtime_override_concat_cpx_vcfs=runtime_override_concat_cpx_vcfs,
      runtime_override_concat_cleaned_vcfs=runtime_override_concat_cleaned_vcfs,
      runtime_override_join_vcfs=runtime_override_join_vcfs,
      runtime_override_subset_bothside_pass=runtime_override_subset_bothside_pass,
      runtime_override_subset_background_fail=runtime_override_subset_background_fail,
      runtime_override_subset_sv_type=runtime_override_subset_sv_type,
      runtime_override_shard_clusters=runtime_override_shard_clusters,
      runtime_override_shard_vids=runtime_override_shard_vids,
      runtime_override_pull_vcf_shard=runtime_override_pull_vcf_shard,
      runtime_override_svtk_vcf_cluster=runtime_override_svtk_vcf_cluster,
      runtime_override_get_vcf_header_with_members_info_line=runtime_override_get_vcf_header_with_members_info_line,
      runtime_override_cluster_merge=runtime_override_cluster_merge,
      runtime_override_concat_vcf_cluster=runtime_override_concat_vcf_cluster,
      runtime_override_concat_svtypes=runtime_override_concat_svtypes,
      runtime_override_concat_sharded_cluster=runtime_override_concat_sharded_cluster,
      runtime_override_make_sites_only=runtime_override_make_sites_only,
      runtime_override_preconcat_sharded_cluster=runtime_override_preconcat_sharded_cluster,
      runtime_override_hail_merge_sharded_cluster=runtime_override_hail_merge_sharded_cluster,
      runtime_override_fix_header_sharded_cluster=runtime_override_fix_header_sharded_cluster,
      runtime_override_concat_large_pesr_depth=runtime_override_concat_large_pesr_depth,
      runtime_override_update_sr_list_pass=runtime_override_update_sr_list_pass,
      runtime_override_update_sr_list_fail=runtime_override_update_sr_list_fail,
      runtime_override_subset_inversions=runtime_override_subset_inversions,
      runtime_override_concat_resolve=runtime_override_concat_resolve,
      runtime_override_get_se_cutoff=runtime_override_get_se_cutoff,
      runtime_override_shard_vcf_cpx=runtime_override_shard_vcf_cpx,
      runtime_override_shard_vids_resolve=runtime_override_shard_vids_resolve,
      runtime_override_resolve_prep=runtime_override_resolve_prep,
      runtime_override_resolve_cpx_per_shard=runtime_override_resolve_cpx_per_shard,
      runtime_override_restore_unresolved_cnv_per_shard=runtime_override_restore_unresolved_cnv_per_shard,
      runtime_override_concat_resolved_per_shard=runtime_override_concat_resolved_per_shard,
      runtime_override_preconcat_resolve=runtime_override_preconcat_resolve,
      runtime_override_hail_merge_resolve=runtime_override_hail_merge_resolve,
      runtime_override_fix_header_resolve=runtime_override_fix_header_resolve,
      runtime_override_get_se_cutoff_inv=runtime_override_get_se_cutoff_inv,
      runtime_override_shard_vcf_cpx_inv=runtime_override_shard_vcf_cpx_inv,
      runtime_override_shard_vids_resolve_inv=runtime_override_shard_vids_resolve_inv,
      runtime_override_resolve_prep_inv=runtime_override_resolve_prep_inv,
      runtime_override_resolve_cpx_per_shard_inv=runtime_override_resolve_cpx_per_shard_inv,
      runtime_override_restore_unresolved_cnv_per_shard_inv=runtime_override_restore_unresolved_cnv_per_shard_inv,
      runtime_override_concat_resolved_per_shard_inv=runtime_override_concat_resolved_per_shard_inv,
      runtime_override_pull_vcf_shard_inv=runtime_override_pull_vcf_shard_inv,
      runtime_override_preconcat_resolve_inv=runtime_override_preconcat_resolve_inv,
      runtime_override_hail_merge_resolve_inv=runtime_override_hail_merge_resolve_inv,
      runtime_override_fix_header_resolve_inv=runtime_override_fix_header_resolve_inv,
      runtime_override_ids_from_median=runtime_override_ids_from_median,
      runtime_override_split_vcf_to_genotype=runtime_override_split_vcf_to_genotype,
      runtime_override_concat_cpx_cnv_vcfs=runtime_override_concat_cpx_cnv_vcfs,
      runtime_override_get_cpx_cnv_intervals=runtime_override_get_cpx_cnv_intervals,
      runtime_override_parse_genotypes=runtime_override_parse_genotypes,
      runtime_override_merge_melted_gts=runtime_override_merge_melted_gts,
      runtime_override_split_bed_by_size=runtime_override_split_bed_by_size,
      runtime_override_rd_genotype=runtime_override_rd_genotype,
      runtime_override_concat_melted_genotypes=runtime_override_concat_melted_genotypes,
      runtime_attr_ids_from_vcf_regeno=runtime_attr_ids_from_vcf_regeno,
      runtime_attr_subset_ped_regeno=runtime_attr_subset_ped_regeno,
      runtime_override_preconcat_regeno=runtime_override_preconcat_regeno,
      runtime_override_hail_merge_regeno=runtime_override_hail_merge_regeno,
      runtime_override_fix_header_regeno=runtime_override_fix_header_regeno,
      runtime_override_preconcat_clean_final=runtime_override_preconcat_clean_final,
      runtime_override_hail_merge_clean_final=runtime_override_hail_merge_clean_final,
      runtime_override_fix_header_clean_final=runtime_override_fix_header_clean_final,
      runtime_override_clean_vcf_1a=runtime_override_clean_vcf_1a,
      runtime_override_clean_vcf_2=runtime_override_clean_vcf_2,
      runtime_override_clean_vcf_3=runtime_override_clean_vcf_3,
      runtime_override_clean_vcf_4=runtime_override_clean_vcf_4,
      runtime_override_clean_vcf_5_scatter=runtime_override_clean_vcf_5_scatter,
      runtime_override_clean_vcf_5_make_cleangq=runtime_override_clean_vcf_5_make_cleangq,
      runtime_override_clean_vcf_5_find_redundant_multiallelics=runtime_override_clean_vcf_5_find_redundant_multiallelics,
      runtime_override_clean_vcf_5_polish=runtime_override_clean_vcf_5_polish,
      runtime_override_stitch_fragmented_cnvs=runtime_override_stitch_fragmented_cnvs,
      runtime_override_final_cleanup=runtime_override_final_cleanup,
      runtime_attr_override_subset_large_cnvs_1b=runtime_attr_override_subset_large_cnvs_1b,
      runtime_attr_override_sort_bed_1b=runtime_attr_override_sort_bed_1b,
      runtime_attr_override_intersect_bed_1b=runtime_attr_override_intersect_bed_1b,
      runtime_attr_override_build_dict_1b=runtime_attr_override_build_dict_1b,
      runtime_attr_override_scatter_1b=runtime_attr_override_scatter_1b,
      runtime_attr_override_filter_vcf_1b=runtime_attr_override_filter_vcf_1b,
      runtime_override_concat_vcfs_1b=runtime_override_concat_vcfs_1b,
      runtime_override_cat_multi_cnvs_1b=runtime_override_cat_multi_cnvs_1b,
      runtime_override_preconcat_step1=runtime_override_preconcat_step1,
      runtime_override_hail_merge_step1=runtime_override_hail_merge_step1,
      runtime_override_fix_header_step1=runtime_override_fix_header_step1,
      runtime_override_preconcat_drc=runtime_override_preconcat_drc,
      runtime_override_hail_merge_drc=runtime_override_hail_merge_drc,
      runtime_override_fix_header_drc=runtime_override_fix_header_drc,
      runtime_override_split_vcf_to_clean=runtime_override_split_vcf_to_clean,
      runtime_override_combine_step_1_sex_chr_revisions=runtime_override_combine_step_1_sex_chr_revisions,
      runtime_override_split_include_list=runtime_override_split_include_list,
      runtime_override_combine_clean_vcf_2=runtime_override_combine_clean_vcf_2,
      runtime_override_combine_revised_4=runtime_override_combine_revised_4,
      runtime_override_combine_multi_ids_4=runtime_override_combine_multi_ids_4,
      runtime_override_drop_redundant_cnvs=runtime_override_drop_redundant_cnvs,
      runtime_override_combine_step_1_vcfs=runtime_override_combine_step_1_vcfs,
      runtime_override_sort_drop_redundant_cnvs=runtime_override_sort_drop_redundant_cnvs,
      runtime_override_site_level_benchmark_plot=runtime_override_site_level_benchmark_plot,
      runtime_override_per_sample_benchmark_plot=runtime_override_per_sample_benchmark_plot,
      runtime_override_subset_vcf=runtime_override_subset_vcf,
      runtime_override_preprocess_vcf=runtime_override_preprocess_vcf,
      runtime_override_site_level_benchmark=runtime_override_site_level_benchmark,
      runtime_override_merge_site_level_benchmark=runtime_override_merge_site_level_benchmark,
      runtime_override_merge_sharded_per_sample_vid_lists=runtime_override_merge_sharded_per_sample_vid_lists,
      runtime_override_plot_qc_vcf_wide=runtime_override_plot_qc_vcf_wide,
      runtime_override_plot_qc_per_sample=runtime_override_plot_qc_per_sample,
      runtime_override_plot_qc_per_family=runtime_override_plot_qc_per_family,
      runtime_override_sanitize_outputs=runtime_override_sanitize_outputs,
      runtime_override_merge_vcfwide_stat_shards=runtime_override_merge_vcfwide_stat_shards,
      runtime_override_merge_vcf_2_bed=runtime_override_merge_vcf_2_bed,
      runtime_override_collect_sharded_vcf_stats=runtime_override_collect_sharded_vcf_stats,
      runtime_override_svtk_vcf_2_bed=runtime_override_svtk_vcf_2_bed,
      runtime_override_scatter_vcf=runtime_override_scatter_vcf,
      runtime_override_merge_subvcf_stat_shards=runtime_override_merge_subvcf_stat_shards,
      runtime_override_collect_vids_per_sample=runtime_override_collect_vids_per_sample,
      runtime_override_split_samples_list=runtime_override_split_samples_list,
      runtime_override_tar_shard_vid_lists=runtime_override_tar_shard_vid_lists,
      runtime_override_benchmark_samples=runtime_override_benchmark_samples,
      runtime_override_split_shuffled_list=runtime_override_split_shuffled_list,
      runtime_override_merge_and_tar_shard_benchmarks=runtime_override_merge_and_tar_shard_benchmarks,
      runtime_attr_format_clean=runtime_attr_format_clean

  }

  call SingleSampleFiltering.FilterVcfForShortDepthCalls as FilterVcfDepthLt5kb {
    input:
      vcf_gz=MakeCohortVcf.vcf,
      min_length=5000,
      filter_name="DEPTH_LT_5KB",
      sv_base_mini_docker=sv_base_mini_docker
  }

  call SingleSampleFiltering.GetUniqueNonGenotypedDepthCalls as GetUniqueNonGenotypedDepthCalls {
    input:
      vcf_gz=select_first([MakeCohortVcf.complex_genotype_vcf]),
      sample_id=sample_id,
      ref_panel_dels=ref_panel_del_bed,
      ref_panel_dups=ref_panel_dup_bed,
      sv_base_mini_docker=sv_base_mini_docker
  }

  call SingleSampleFiltering.FilterVcfForCaseSampleGenotype as FilterVcfForCaseSampleGenotype {
    input:
      vcf_gz=FilterVcfDepthLt5kb.out,
      sample_id=sample_id,
      sv_base_mini_docker=sv_base_mini_docker
  }

  call SingleSampleFiltering.FilterVcfWithReferencePanelCalls as FilterVcfWithReferencePanelCalls {
    input:
      single_sample_vcf=FilterVcfForCaseSampleGenotype.out,
      cohort_vcf=ref_panel_vcf,
      case_sample_id=sample_id,
      max_ref_panel_carrier_freq=max_ref_panel_carrier_freq,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call SingleSampleMetrics.SingleSampleMetrics as SampleFilterMetrics {
    input:
      name = batch,
      ref_samples = ref_samples,
      case_sample = sample_id,
      wgd_scores = EvidenceQC.WGD_scores,
      sample_counts = case_counts_file_,
      contig_list = primary_contigs_list,
      linux_docker = linux_docker,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call utils.RunQC as SampleFilterQC {
    input:
      name=batch,
      metrics=SampleFilterMetrics.metrics_file,
      qc_definitions = qc_definitions,
      sv_pipeline_base_docker=sv_pipeline_base_docker
  }

  call SingleSampleFiltering.SampleQC as FilterSample {
    input:
      vcf=FilterVcfWithReferencePanelCalls.out,
      sample_filtering_qc_file=SampleFilterQC.out,
      sv_pipeline_base_docker=sv_pipeline_base_docker,
  }

  call annotate.AnnotateVcf {
       input:
        vcf = FilterSample.out,
        prefix = batch,
        contig_list = primary_contigs_list,
        protein_coding_gtf = protein_coding_gtf,
        noncoding_bed = noncoding_bed,
        promoter_window = promoter_window,
        max_breakend_as_cnv_length = max_breakend_as_cnv_length,
        ref_bed = external_af_ref_bed,
        ref_prefix = external_af_ref_bed_prefix,
        population = external_af_population,
        use_hail = false,
        sv_per_shard = annotation_sv_per_shard,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        gatk_docker = gatk_docker,
        runtime_attr_svannotate = runtime_attr_svannotate
  }

  call SingleSampleFiltering.VcfToBed as VcfToBed {
    input:
      vcf = AnnotateVcf.annotated_vcf,
      prefix = batch,
      sv_pipeline_docker = sv_pipeline_docker
  }

  call SingleSampleFiltering.UpdateBreakendRepresentation {
    input:
      vcf=AnnotateVcf.annotated_vcf,
      vcf_idx=AnnotateVcf.annotated_vcf_index,
      ref_fasta=reference_fasta,
      ref_fasta_idx=reference_index,
      prefix=basename(AnnotateVcf.annotated_vcf, ".vcf.gz") + ".final_cleanup",
      sv_pipeline_docker=sv_pipeline_docker
  }

  call SingleSampleMetrics.SingleSampleMetrics {
    input:
      name = batch,
      ref_samples = ref_samples,
      case_sample = sample_id,
      wgd_scores = EvidenceQC.WGD_scores,
      sample_pe = case_pe_file_,
      sample_sr = case_sr_file_,
      sample_counts = case_counts_file_,
      cleaned_vcf = MakeCohortVcf.vcf,
      final_vcf = UpdateBreakendRepresentation.out,
      genotyped_pesr_vcf = ConvertCNVsWithoutDepthSupportToBNDs.out_vcf,
      genotyped_depth_vcf = GenotypeBatch.genotyped_depth_vcf,
      non_genotyped_unique_depth_calls_vcf = GetUniqueNonGenotypedDepthCalls.out,
      contig_list = primary_contigs_list,
      linux_docker = linux_docker,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call utils.RunQC as SingleSampleQC {
    input:
      name = batch,
      metrics = SingleSampleMetrics.metrics_file,
      qc_definitions = qc_definitions,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  output {
    File final_vcf = UpdateBreakendRepresentation.out
    File final_vcf_idx = UpdateBreakendRepresentation.out_idx

    File final_bed = VcfToBed.bed

    # These files contain events reported in the internal VCF representation
    # They are less VCF-spec compliant but may be useful if components of the pipeline need to be re-run
    # on the output.
    File pre_cleanup_vcf = AnnotateVcf.annotated_vcf
    File pre_cleanup_vcf_idx = AnnotateVcf.annotated_vcf_index

    File ploidy_matrix = select_first([GatherBatchEvidence.batch_ploidy_matrix])
    File ploidy_plots = select_first([GatherBatchEvidence.batch_ploidy_plots])
    File metrics_file = SingleSampleMetrics.metrics_file
    File qc_file = SingleSampleQC.out

    # These files contain any depth based calls made in the case sample that did not pass genotyping
    # in the case sample and do not match a depth-based call from the reference panel.
    File non_genotyped_unique_depth_calls = GetUniqueNonGenotypedDepthCalls.out
    File non_genotyped_unique_depth_calls_idx = GetUniqueNonGenotypedDepthCalls.out_idx
  }
}
