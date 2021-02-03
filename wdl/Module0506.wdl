version 1.0

import "Module0506Cluster.wdl" as Cluster
import "Module0506ComplexResolve.wdl" as ComplexResolve
import "Module0506ComplexGenotype.wdl" as ComplexGenotype
import "Module0506Clean.wdl" as Clean
import "MasterVcfQc.wdl" as VcfQc

workflow Module0506 {
  input {
    String cohort_name
    Array[String] batches
    Array[File] ped_files

    # Merge contig vcfs at each stage for QC
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

    File bin_exclude
    File contig_list
    Int max_shards_per_chrom
    Int min_variants_per_shard
    File cytobands
    File mei_bed
    File pe_exclude_list
    File depth_exclude_list
    File ref_dict
    Int max_shards_per_chrom_clean_vcf_step1
    Int min_records_per_shard_clean_vcf_step1
    Int samples_per_clean_vcf_step2_shard
    Float min_sr_background_fail_batches

    File empty_file
    File? outlier_samples_list
    Int? random_seed

    Array[File]? thousand_genomes_tarballs
    Array[File]? hgsv_tarballs
    Array[File]? asc_tarballs
    File? sanders_2015_tarball
    File? collins_2017_tarball
    File? werling_2018_tarball

    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker
    String sv_pipeline_qc_docker

    # overrides for local tasks
    RuntimeAttr? runtime_overide_get_discfile_size
    RuntimeAttr? runtime_override_update_sr_list
    RuntimeAttr? runtime_override_merge_pesr_depth
    RuntimeAttr? runtime_override_breakpoint_overlap_filter
    RuntimeAttr? runtime_override_integrate_resolved_vcfs
    RuntimeAttr? runtime_override_rename_variants
    RuntimeAttr? runtime_override_rename_cleaned_samples

    # overrides for mini tasks
    RuntimeAttr? runtime_override_ids_from_vcf
    RuntimeAttr? runtime_override_clean_bothside_pass
    RuntimeAttr? runtime_override_clean_background_fail
    RuntimeAttr? runtime_override_merge_fam_file_list
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
    RuntimeAttr? runtime_override_concat_sv_types
    RuntimeAttr? runtime_override_shard_vcf_precluster
    RuntimeAttr? runtime_override_svtk_vcf_cluster
    RuntimeAttr? runtime_override_get_vcf_header_with_members_info_line
    RuntimeAttr? runtime_override_cluster_merge
    RuntimeAttr? runtime_override_concat_shards

    # overrides for ResolveComplexContig
    RuntimeAttr? runtime_override_get_se_cutoff
    RuntimeAttr? runtime_override_shard_vcf_cpx
    RuntimeAttr? runtime_override_resolve_prep
    RuntimeAttr? runtime_override_resolve_cpx_per_shard
    RuntimeAttr? runtime_override_restore_unresolved_cnv_per_shard
    RuntimeAttr? runtime_override_concat_resolved_per_shard
    RuntimeAttr? runtime_override_complex_resolve_merge

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
    RuntimeAttr? runtime_override_complex_genotype_merge

    # overrides for CleanVcfContig
    RuntimeAttr? runtime_override_clean_vcf_1a
    RuntimeAttr? runtime_override_clean_vcf_1b
    RuntimeAttr? runtime_override_clean_vcf_2
    RuntimeAttr? runtime_override_clean_vcf_3
    RuntimeAttr? runtime_override_clean_vcf_4
    RuntimeAttr? runtime_override_clean_vcf_5
    RuntimeAttr? runtime_override_drop_redundant_cnvs
    RuntimeAttr? runtime_override_stitch_fragmented_cnvs
    RuntimeAttr? runtime_override_final_cleanup
    RuntimeAttr? runtime_override_split_vcf_to_clean
    RuntimeAttr? runtime_override_combine_step_1_vcfs
    RuntimeAttr? runtime_override_combine_step_1_sex_chr_revisions
    RuntimeAttr? runtime_override_split_include_list
    RuntimeAttr? runtime_override_combine_clean_vcf_2
    RuntimeAttr? runtime_override_combine_revised_4
    RuntimeAttr? runtime_override_combine_multi_ids_4

    # overrides for VcfQc
    RuntimeAttr? runtime_override_plot_qc_vcf_wide
    RuntimeAttr? runtime_override_thousand_g_benchmark
    RuntimeAttr? runtime_override_thousand_g_plot
    RuntimeAttr? runtime_override_asc_benchmark
    RuntimeAttr? runtime_override_asc_plot
    RuntimeAttr? runtime_override_hgsv_benchmark
    RuntimeAttr? runtime_override_hgsv_plot
    RuntimeAttr? runtime_override_plot_qc_per_sample
    RuntimeAttr? runtime_override_plot_qc_per_family
    RuntimeAttr? runtime_override_sanders_per_sample_plot
    RuntimeAttr? runtime_override_collins_per_sample_plot
    RuntimeAttr? runtime_override_werling_per_sample_plot
    RuntimeAttr? runtime_override_sanitize_outputs
    RuntimeAttr? runtime_override_merge_vcfwide_stat_shards
    RuntimeAttr? runtime_override_merge_vcf_2_bed
    RuntimeAttr? runtime_override_collect_sharded_vcf_stats
    RuntimeAttr? runtime_override_svtk_vcf_2_bed
    RuntimeAttr? runtime_override_split_vcf_to_qc
    RuntimeAttr? runtime_override_merge_subvcf_stat_shards
    RuntimeAttr? runtime_override_merge_svtk_vcf_2_bed
    RuntimeAttr? runtime_override_collect_vids_per_sample
    RuntimeAttr? runtime_override_split_samples_list
    RuntimeAttr? runtime_override_tar_shard_vid_lists
    RuntimeAttr? runtime_override_benchmark_samples
    RuntimeAttr? runtime_override_split_shuffled_list
    RuntimeAttr? runtime_override_merge_and_tar_shard_benchmarks
  }

  call Cluster.Module0506Cluster {
    input:
      cohort_name=cohort_name,
      batches=batches,
      ped_files=ped_files,
      merge_vcfs=merge_cluster_vcfs,
      pesr_vcfs=pesr_vcfs,
      depth_vcfs=depth_vcfs,
      raw_sr_bothside_pass_files=raw_sr_bothside_pass_files,
      raw_sr_background_fail_files=raw_sr_background_fail_files,
      contig_list=contig_list,
      max_shards_per_chrom=max_shards_per_chrom,
      min_variants_per_shard=min_variants_per_shard,
      pe_exclude_list=pe_exclude_list,
      depth_exclude_list=depth_exclude_list,
      min_sr_background_fail_batches=min_sr_background_fail_batches,
      empty_file=empty_file,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_override_update_sr_list=runtime_override_update_sr_list,
      runtime_override_merge_pesr_depth=runtime_override_merge_pesr_depth,
      runtime_override_clean_bothside_pass=runtime_override_clean_bothside_pass,
      runtime_override_clean_background_fail=runtime_override_clean_background_fail,
      runtime_override_merge_fam_file_list=runtime_override_merge_fam_file_list,
      runtime_override_concat=runtime_override_cluster_merge,
      runtime_override_join_vcfs=runtime_override_join_vcfs,
      runtime_override_subset_bothside_pass=runtime_override_subset_bothside_pass,
      runtime_override_subset_background_fail=runtime_override_subset_background_fail,
      runtime_override_subset_sv_type=runtime_override_subset_sv_type,
      runtime_override_concat_sv_types=runtime_override_concat_sv_types,
      runtime_override_shard_vcf_precluster=runtime_override_shard_vcf_precluster,
      runtime_override_svtk_vcf_cluster=runtime_override_svtk_vcf_cluster,
      runtime_override_get_vcf_header_with_members_info_line=runtime_override_get_vcf_header_with_members_info_line,
      runtime_override_concat_shards=runtime_override_concat_shards
  }

  call ComplexResolve.Module0506ComplexResolve {
    input:
      cohort_name=cohort_name,
      merge_vcfs=merge_complex_resolve_vcfs,
      cluster_vcfs=Module0506Cluster.vcfs,
      cluster_bothside_pass_lists=Module0506Cluster.cluster_bothside_pass_lists,
      cluster_background_fail_lists=Module0506Cluster.cluster_background_fail_lists,
      disc_files=disc_files,
      rf_cutoff_files=rf_cutoff_files,
      contig_list=contig_list,
      max_shards_per_chrom=max_shards_per_chrom,
      min_variants_per_shard=min_variants_per_shard,
      cytobands=cytobands,
      mei_bed=mei_bed,
      pe_exclude_list=pe_exclude_list,
      ref_dict=ref_dict,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_override_update_sr_list=runtime_override_update_sr_list,
      runtime_override_breakpoint_overlap_filter=runtime_override_breakpoint_overlap_filter,
      runtime_override_integrate_resolved_vcfs=runtime_override_integrate_resolved_vcfs,
      runtime_override_rename_variants=runtime_override_rename_variants,
      runtime_override_subset_inversions=runtime_override_subset_inversions,
      runtime_override_merge_fam_file_list=runtime_override_merge_fam_file_list,
      runtime_override_concat=runtime_override_complex_resolve_merge,
      runtime_override_get_se_cutoff=runtime_override_get_se_cutoff,
      runtime_override_shard_vcf_cpx=runtime_override_shard_vcf_cpx,
      runtime_override_resolve_prep=runtime_override_resolve_prep,
      runtime_override_resolve_cpx_per_shard=runtime_override_resolve_cpx_per_shard,
      runtime_override_restore_unresolved_cnv_per_shard=runtime_override_restore_unresolved_cnv_per_shard,
      runtime_override_concat_resolved_per_shard=runtime_override_concat_resolved_per_shard
  }

  call ComplexGenotype.Module0506ComplexGenotype {
    input:
      cohort_name=cohort_name,
      batches=batches,
      merge_vcfs=merge_complex_genotype_vcfs,
      complex_resolve_vcfs=Module0506ComplexResolve.complex_resolve_vcfs,
      complex_resolve_vcf_indexes=Module0506ComplexResolve.complex_resolve_vcf_indexes,
      merged_ped_file=Module0506Cluster.merged_ped_file,
      ped_files=ped_files,
      bincov_files=bincov_files,
      depth_gt_rd_sep_files=depth_gt_rd_sep_files,
      median_coverage_files=median_coverage_files,
      bin_exclude=bin_exclude,
      contig_list=contig_list,
      max_shards_per_chrom=max_shards_per_chrom,
      ref_dict=ref_dict,
      linux_docker=linux_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
      runtime_override_ids_from_vcf=runtime_override_ids_from_vcf,
      runtime_override_merge_fam_file_list=runtime_override_merge_fam_file_list,
      runtime_override_concat=runtime_override_complex_genotype_merge,
      runtime_override_ids_from_median=runtime_override_ids_from_median,
      runtime_override_split_vcf_to_genotype=runtime_override_split_vcf_to_genotype,
      runtime_override_concat_cpx_cnv_vcfs=runtime_override_concat_cpx_cnv_vcfs,
      runtime_override_get_cpx_cnv_intervals=runtime_override_get_cpx_cnv_intervals,
      runtime_override_parse_genotypes=runtime_override_parse_genotypes,
      runtime_override_merge_melted_gts=runtime_override_merge_melted_gts,
      runtime_override_split_bed_by_size=runtime_override_split_bed_by_size,
      runtime_override_rd_genotype=runtime_override_rd_genotype,
      runtime_override_concat_melted_genotypes=runtime_override_concat_melted_genotypes
  }

  call Clean.Module0506Clean {
    input:
      cohort_name=cohort_name,
      complex_genotype_vcfs=Module0506ComplexGenotype.complex_genotype_vcfs,
      complex_resolve_bothside_pass_lists=Module0506ComplexResolve.complex_resolve_bothside_pass_lists,
      complex_resolve_background_fail_lists=Module0506ComplexResolve.complex_resolve_background_fail_lists,
      merged_ped_file=Module0506Cluster.merged_ped_file,
      contig_list=contig_list,
      max_shards_per_chrom=max_shards_per_chrom,
      max_shards_per_chrom_clean_vcf_step1=max_shards_per_chrom_clean_vcf_step1,
      min_records_per_shard_clean_vcf_step1=min_records_per_shard_clean_vcf_step1,
      samples_per_clean_vcf_step2_shard=samples_per_clean_vcf_step2_shard,
      outlier_samples_list=outlier_samples_list,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_override_concat_cleaned_vcfs=runtime_override_concat_cleaned_vcfs,
      runtime_override_clean_vcf_1a=runtime_override_clean_vcf_1a,
      runtime_override_clean_vcf_1b=runtime_override_clean_vcf_1b,
      runtime_override_clean_vcf_2=runtime_override_clean_vcf_2,
      runtime_override_clean_vcf_3=runtime_override_clean_vcf_3,
      runtime_override_clean_vcf_4=runtime_override_clean_vcf_4,
      runtime_override_clean_vcf_5=runtime_override_clean_vcf_5,
      runtime_override_drop_redundant_cnvs=runtime_override_drop_redundant_cnvs,
      runtime_override_stitch_fragmented_cnvs=runtime_override_stitch_fragmented_cnvs,
      runtime_override_final_cleanup=runtime_override_final_cleanup,
      runtime_override_split_vcf_to_clean=runtime_override_split_vcf_to_clean,
      runtime_override_combine_step_1_vcfs=runtime_override_combine_step_1_vcfs,
      runtime_override_combine_step_1_sex_chr_revisions=runtime_override_combine_step_1_sex_chr_revisions,
      runtime_override_split_include_list=runtime_override_split_include_list,
      runtime_override_combine_clean_vcf_2=runtime_override_combine_clean_vcf_2,
      runtime_override_combine_revised_4=runtime_override_combine_revised_4,
      runtime_override_combine_multi_ids_4=runtime_override_combine_multi_ids_4
  }

  Array[String] contigs = transpose(read_tsv(contig_list))[0]
  call VcfQc.MasterVcfQc {
    input:
      vcf=Module0506Clean.cleaned_vcf,
      vcf_idx=Module0506Clean.cleaned_vcf_index,
      ped_file=Module0506Cluster.merged_ped_file,
      prefix="~{cohort_name}.cleaned",
      sv_per_shard=10000,
      samples_per_shard=100,
      thousand_genomes_tarballs=thousand_genomes_tarballs,
      hgsv_tarballs=hgsv_tarballs,
      asc_tarballs=asc_tarballs,
      sanders_2015_tarball=sanders_2015_tarball,
      collins_2017_tarball=collins_2017_tarball,
      werling_2018_tarball=werling_2018_tarball,
      contigs=contigs,
      random_seed=random_seed,
      sv_pipeline_qc_docker=sv_pipeline_qc_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_override_collect_vids_per_sample=runtime_override_collect_vids_per_sample
  }

  output {
    File vcf = Module0506Clean.cleaned_vcf
    File vcf_index = Module0506Clean.cleaned_vcf_index
    File vcf_qc = MasterVcfQc.sv_vcf_qc_output

    # If merge_intermediate_vcfs enabled
    File? cluster_vcf = Module0506Cluster.merged_vcf
    File? cluster_vcf_index = Module0506Cluster.merged_vcf_index
    File? complex_resolve_vcf = Module0506ComplexResolve.merged_vcf
    File? complex_resolve_vcf_index = Module0506ComplexResolve.merged_vcf_index
    File? complex_genotype_vcf = Module0506ComplexGenotype.merged_vcf
    File? complex_genotype_vcf_index = Module0506ComplexGenotype.merged_vcf_index
  }
}
