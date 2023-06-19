version 1.0

import "CombineBatches.wdl" as Cluster
import "ResolveComplexVariants.wdl" as ComplexResolve
import "GenotypeComplexVariants.wdl" as ComplexGenotype
import "CleanVcf.wdl" as Clean
import "MainVcfQc.wdl" as VcfQc

workflow MakeCohortVcf {
  input {
    String cohort_name
    Array[String] batches
    File ped_file # cohort ped file

    # Merge contig vcfs at each stage for QC
    # Not recommended for very large cohorts
    Boolean merge_cluster_vcfs = false
    Boolean merge_complex_resolve_vcfs = false
    Boolean merge_complex_genotype_vcfs = false

    Array[File] pesr_vcfs
    Array[File] depth_vcfs
    Array[File] disc_files
    Array[File] bincov_files

    Array[File] raw_sr_bothside_pass_files
    Array[File] raw_sr_background_fail_files
    Array[File] depth_gt_rd_sep_files
    Array[File] median_coverage_files
    Array[File] rf_cutoff_files

    # Enables use of Hail for merging and sorting VCFs
    # Recommended for cohorts of 10,000 samples or more
    # Requires that DataProc be enabled in the GCP project
    Boolean use_hail = false
    String? gcs_project

    File bin_exclude
    File contig_list
    File allosome_fai
    Int? localize_shard_size
    File cytobands
    File mei_bed
    File pe_exclude_list
    File depth_exclude_list
    File ref_dict
    Int max_shard_size_resolve
    Int max_shards_per_chrom_clean_vcf_step1
    Int min_records_per_shard_clean_vcf_step1
    Int clean_vcf1b_records_per_shard
    Int samples_per_clean_vcf_step2_shard
    Int clean_vcf5_records_per_shard
    Int? clean_vcf5_threads_per_task
    Float min_sr_background_fail_batches
    Int? max_samples_per_shard_clean_vcf_step3

    String chr_x
    String chr_y

    File empty_file
    File? outlier_samples_list
    Int? random_seed
    Int? max_gq  # Max GQ for plotting. Default = 99, ie. GQ is on a scale of [0,99]. Prior to CleanVcf, use 999

    Array[Array[String]]? site_level_comparison_datasets    # Array of two-element arrays, one per dataset, each of format [prefix, gs:// path to directory with one BED per population]
    Array[Array[String]]? sample_level_comparison_datasets  # Array of two-element arrays, one per dataset, each of format [prefix, gs:// path to per-sample tarballs]
    File? sample_renaming_tsv # File with mapping to rename sample IDs for compatibility with sample_level_comparison_datasets

    # Module metrics parameters
    # Run module metrics workflow at the end - on by default
    Boolean? run_module_metrics
    String? sv_pipeline_base_docker  # required if run_module_metrics = true
    File? primary_contigs_list  # required if run_module_metrics = true
    File? baseline_cluster_vcf  # baseline files are optional for metrics workflow
    File? baseline_complex_resolve_vcf
    File? baseline_complex_genotype_vcf
    File? baseline_cleaned_vcf

    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_hail_docker
    String sv_pipeline_updates_docker
    String sv_pipeline_rdtest_docker
    String sv_pipeline_qc_docker

    # overrides for local tasks
    RuntimeAttr? runtime_overide_get_discfile_size
    RuntimeAttr? runtime_override_update_sr_list_cluster
    RuntimeAttr? runtime_override_merge_pesr_depth
    RuntimeAttr? runtime_override_integrate_resolved_vcfs
    RuntimeAttr? runtime_override_rename_variants
    RuntimeAttr? runtime_override_rename_cleaned_samples

    RuntimeAttr? runtime_override_breakpoint_overlap_filter

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
    RuntimeAttr? runtime_attr_format_clean

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
  }

  call Cluster.CombineBatches {
    input:
      cohort_name=cohort_name,
      batches=batches,
      merge_vcfs=merge_cluster_vcfs,
      pesr_vcfs=pesr_vcfs,
      depth_vcfs=depth_vcfs,
      raw_sr_bothside_pass_files=raw_sr_bothside_pass_files,
      raw_sr_background_fail_files=raw_sr_background_fail_files,
      contig_list=contig_list,
      localize_shard_size=localize_shard_size,
      pe_exclude_list=pe_exclude_list,
      depth_exclude_list=depth_exclude_list,
      min_sr_background_fail_batches=min_sr_background_fail_batches,
      empty_file=empty_file,
      use_hail=use_hail,
      gcs_project=gcs_project,
      sv_pipeline_hail_docker=sv_pipeline_hail_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_override_update_sr_list=runtime_override_update_sr_list_cluster,
      runtime_override_merge_pesr_depth=runtime_override_merge_pesr_depth,
      runtime_override_clean_background_fail=runtime_override_clean_background_fail,
      runtime_override_concat=runtime_override_cluster_merge,
      runtime_override_join_vcfs=runtime_override_join_vcfs,
      runtime_override_subset_bothside_pass=runtime_override_subset_bothside_pass,
      runtime_override_subset_background_fail=runtime_override_subset_background_fail,
      runtime_override_subset_sv_type=runtime_override_subset_sv_type,
      runtime_override_shard_clusters=runtime_override_shard_clusters,
      runtime_override_shard_vids=runtime_override_shard_vids,
      runtime_override_pull_vcf_shard=runtime_override_pull_vcf_shard,
      runtime_override_svtk_vcf_cluster=runtime_override_svtk_vcf_cluster,
      runtime_override_get_vcf_header_with_members_info_line=runtime_override_get_vcf_header_with_members_info_line,
      runtime_override_concat_vcf_cluster=runtime_override_concat_vcf_cluster,
      runtime_override_concat_svtypes=runtime_override_concat_svtypes,
      runtime_override_concat_sharded_cluster=runtime_override_concat_sharded_cluster,
      runtime_override_make_sites_only=runtime_override_make_sites_only,
      runtime_override_preconcat_sharded_cluster=runtime_override_preconcat_sharded_cluster,
      runtime_override_hail_merge_sharded_cluster=runtime_override_hail_merge_sharded_cluster,
      runtime_override_fix_header_sharded_cluster=runtime_override_fix_header_sharded_cluster,
      runtime_override_concat_large_pesr_depth=runtime_override_concat_large_pesr_depth
  }

  call ComplexResolve.ResolveComplexVariants {
    input:
      cohort_name=cohort_name,
      merge_vcfs=merge_complex_resolve_vcfs,
      cluster_vcfs=CombineBatches.combined_vcfs,
      cluster_bothside_pass_lists=CombineBatches.cluster_bothside_pass_lists,
      cluster_background_fail_lists=CombineBatches.cluster_background_fail_lists,
      disc_files=disc_files,
      rf_cutoff_files=rf_cutoff_files,
      contig_list=contig_list,
      cytobands=cytobands,
      mei_bed=mei_bed,
      pe_exclude_list=pe_exclude_list,
      ref_dict=ref_dict,
      use_hail=use_hail,
      gcs_project=gcs_project,
      sv_pipeline_hail_docker=sv_pipeline_hail_docker,
      max_shard_size=max_shard_size_resolve,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_override_update_sr_list_pass=runtime_override_update_sr_list_pass,
      runtime_override_update_sr_list_fail=runtime_override_update_sr_list_fail,
      runtime_override_integrate_resolved_vcfs=runtime_override_integrate_resolved_vcfs,
      runtime_override_rename_variants=runtime_override_rename_variants,
      runtime_override_breakpoint_overlap_filter=runtime_override_breakpoint_overlap_filter,
      runtime_override_subset_inversions=runtime_override_subset_inversions,
      runtime_override_concat=runtime_override_concat_resolve,

      runtime_override_get_se_cutoff=runtime_override_get_se_cutoff,
      runtime_override_shard_vcf_cpx=runtime_override_shard_vcf_cpx,
      runtime_override_shard_vids=runtime_override_shard_vids_resolve,
      runtime_override_resolve_prep=runtime_override_resolve_prep,
      runtime_override_resolve_cpx_per_shard=runtime_override_resolve_cpx_per_shard,
      runtime_override_restore_unresolved_cnv_per_shard=runtime_override_restore_unresolved_cnv_per_shard,
      runtime_override_concat_resolved_per_shard=runtime_override_concat_resolved_per_shard,
      runtime_override_pull_vcf_shard=runtime_override_pull_vcf_shard,
      runtime_override_preconcat=runtime_override_preconcat_resolve,
      runtime_override_hail_merge=runtime_override_hail_merge_resolve,
      runtime_override_fix_header=runtime_override_fix_header_resolve,

      runtime_override_get_se_cutoff_inv=runtime_override_get_se_cutoff_inv,
      runtime_override_shard_vcf_cpx_inv=runtime_override_shard_vcf_cpx_inv,
      runtime_override_shard_vids_inv=runtime_override_shard_vids_resolve_inv,
      runtime_override_resolve_prep_inv=runtime_override_resolve_prep_inv,
      runtime_override_resolve_cpx_per_shard_inv=runtime_override_resolve_cpx_per_shard_inv,
      runtime_override_restore_unresolved_cnv_per_shard_inv=runtime_override_restore_unresolved_cnv_per_shard_inv,
      runtime_override_concat_resolved_per_shard_inv=runtime_override_concat_resolved_per_shard_inv,
      runtime_override_pull_vcf_shard_inv=runtime_override_pull_vcf_shard_inv,
      runtime_override_preconcat_inv=runtime_override_preconcat_resolve_inv,
      runtime_override_hail_merge_inv=runtime_override_hail_merge_resolve_inv,
      runtime_override_fix_header_inv=runtime_override_fix_header_resolve_inv
  }

  call ComplexGenotype.GenotypeComplexVariants {
    input:
      cohort_name=cohort_name,
      batches=batches,
      merge_vcfs=merge_complex_genotype_vcfs,
      complex_resolve_vcfs=ResolveComplexVariants.complex_resolve_vcfs,
      complex_resolve_vcf_indexes=ResolveComplexVariants.complex_resolve_vcf_indexes,
      depth_vcfs=depth_vcfs,
      ped_file=ped_file,
      bincov_files=bincov_files,
      depth_gt_rd_sep_files=depth_gt_rd_sep_files,
      median_coverage_files=median_coverage_files,
      bin_exclude=bin_exclude,
      contig_list=contig_list,
      ref_dict=ref_dict,
      linux_docker=linux_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_hail_docker=sv_pipeline_hail_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_updates_docker=sv_pipeline_updates_docker,
      sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
      runtime_override_ids_from_median=runtime_override_ids_from_median,
      runtime_override_split_vcf_to_genotype=runtime_override_split_vcf_to_genotype,
      runtime_override_concat_cpx_cnv_vcfs=runtime_override_concat_cpx_cnv_vcfs,
      runtime_override_get_cpx_cnv_intervals=runtime_override_get_cpx_cnv_intervals,
      runtime_override_parse_genotypes=runtime_override_parse_genotypes,
      runtime_override_merge_melted_gts=runtime_override_merge_melted_gts,
      runtime_override_split_bed_by_size=runtime_override_split_bed_by_size,
      runtime_override_rd_genotype=runtime_override_rd_genotype,
      runtime_override_concat_melted_genotypes=runtime_override_concat_melted_genotypes,
      runtime_attr_ids_from_vcf=runtime_attr_ids_from_vcf_regeno,
      runtime_attr_subset_ped=runtime_attr_subset_ped_regeno,
      runtime_override_preconcat=runtime_override_preconcat_regeno,
      runtime_override_hail_merge=runtime_override_hail_merge_regeno,
      runtime_override_fix_header=runtime_override_fix_header_regeno
  }

  call Clean.CleanVcf {
    input:
      cohort_name=cohort_name,
      complex_genotype_vcfs=GenotypeComplexVariants.complex_genotype_vcfs,
      complex_resolve_bothside_pass_lists=ResolveComplexVariants.complex_resolve_bothside_pass_lists,
      complex_resolve_background_fail_lists=ResolveComplexVariants.complex_resolve_background_fail_lists,
      ped_file=ped_file,
      contig_list=contig_list,
      allosome_fai=allosome_fai,
      chr_x=chr_x,
      chr_y=chr_y,
      max_shards_per_chrom_step1=max_shards_per_chrom_clean_vcf_step1,
      min_records_per_shard_step1=min_records_per_shard_clean_vcf_step1,
      clean_vcf1b_records_per_shard=clean_vcf1b_records_per_shard,
      samples_per_step2_shard=samples_per_clean_vcf_step2_shard,
      max_samples_per_shard_step3=max_samples_per_shard_clean_vcf_step3,
      clean_vcf5_records_per_shard=clean_vcf5_records_per_shard,
      outlier_samples_list=outlier_samples_list,
      use_hail=use_hail,
      gcs_project=gcs_project,
      combine_batches_merged_vcf=CombineBatches.combine_batches_merged_vcf,
      resolve_complex_merged_vcf=ResolveComplexVariants.complex_resolve_merged_vcf,
      genotype_complex_merged_vcf=GenotypeComplexVariants.complex_genotype_merged_vcf,
      baseline_cluster_vcf = baseline_cluster_vcf,
      baseline_complex_resolve_vcf = baseline_complex_resolve_vcf,
      baseline_complex_genotype_vcf = baseline_complex_genotype_vcf,
      baseline_cleaned_vcf = baseline_cleaned_vcf,
      primary_contigs_list=primary_contigs_list,
      run_module_metrics=run_module_metrics,
      sv_pipeline_base_docker=sv_pipeline_base_docker,
      linux_docker=linux_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_hail_docker=sv_pipeline_hail_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_updates_docker=sv_pipeline_updates_docker,
      runtime_override_preconcat_clean_final=runtime_override_preconcat_clean_final,
      runtime_override_hail_merge_clean_final=runtime_override_hail_merge_clean_final,
      runtime_override_fix_header_clean_final=runtime_override_fix_header_clean_final,
      runtime_override_concat_cleaned_vcfs=runtime_override_concat_cleaned_vcfs,
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
      runtime_attr_format=runtime_attr_format_clean
  }

  call VcfQc.MainVcfQc {
    input:
      vcfs=[CleanVcf.cleaned_vcf],
      ped_file=ped_file,
      prefix="~{cohort_name}.cleaned",
      sv_per_shard=2500,
      samples_per_shard=600,
      site_level_comparison_datasets=site_level_comparison_datasets,
      sample_level_comparison_datasets=sample_level_comparison_datasets,
      sample_renaming_tsv=sample_renaming_tsv,
      primary_contigs_fai=contig_list,
      random_seed=random_seed,
      sv_pipeline_qc_docker=sv_pipeline_qc_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
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
      runtime_override_merge_and_tar_shard_benchmarks=runtime_override_merge_and_tar_shard_benchmarks
  }


  output {
    File vcf = CleanVcf.cleaned_vcf
    File vcf_index = CleanVcf.cleaned_vcf_index
    File vcf_qc = MainVcfQc.sv_vcf_qc_output

    # If merge_intermediate_vcfs enabled
    File? cluster_vcf = CombineBatches.combine_batches_merged_vcf
    File? cluster_vcf_index = CombineBatches.combine_batches_merged_vcf_index
    File? complex_resolve_vcf = ResolveComplexVariants.complex_resolve_merged_vcf
    File? complex_resolve_vcf_index = ResolveComplexVariants.complex_resolve_merged_vcf_index
    File? complex_genotype_vcf = GenotypeComplexVariants.complex_genotype_merged_vcf
    File? complex_genotype_vcf_index = GenotypeComplexVariants.complex_genotype_merged_vcf_index

    # CombineBatches
    Array[File] combined_vcfs = CombineBatches.combined_vcfs
    Array[File] combined_vcf_indexes = CombineBatches.combined_vcf_indexes
    Array[File] cluster_bothside_pass_lists = CombineBatches.cluster_bothside_pass_lists
    Array[File] cluster_background_fail_lists = CombineBatches.cluster_background_fail_lists

    # ResolveComplexVariants
    Array[File] complex_resolve_vcfs = ResolveComplexVariants.complex_resolve_vcfs
    Array[File] complex_resolve_vcf_indexes = ResolveComplexVariants.complex_resolve_vcf_indexes
    Array[File] complex_resolve_bothside_pass_lists = ResolveComplexVariants.complex_resolve_bothside_pass_lists
    Array[File] complex_resolve_background_fail_lists = ResolveComplexVariants.complex_resolve_background_fail_lists
    Array[File] breakpoint_overlap_dropped_record_vcfs = ResolveComplexVariants.breakpoint_overlap_dropped_record_vcfs
    Array[File] breakpoint_overlap_dropped_record_vcf_indexes = ResolveComplexVariants.breakpoint_overlap_dropped_record_vcf_indexes

    # GenotypeComplexVariants
    Array[File] complex_genotype_vcfs = GenotypeComplexVariants.complex_genotype_vcfs
    Array[File] complex_genotype_vcf_indexes = GenotypeComplexVariants.complex_genotype_vcfs

    File? metrics_file_makecohortvcf = CleanVcf.metrics_file_makecohortvcf
  }
}
