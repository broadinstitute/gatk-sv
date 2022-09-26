version 1.0

import "MainVcfQc.wdl" as vcf_qc
import "Utils.wdl" as util

workflow FilterBatchQc {
  input {
    File? manta_vcf_noOutliers
    File? melt_vcf_noOutliers
    File? wham_vcf_noOutliers
    File? depth_vcf_noOutliers
    File? merged_pesr_vcf

    File? manta_vcf_noOutliers_index
    File? melt_vcf_noOutliers_index
    File? wham_vcf_noOutliers_index
    File? depth_vcf_noOutliers_index
    File? merged_pesr_vcf_index

    String batch
    File ped_file
    Array[Array[String]]? site_level_comparison_datasets    # Array of two-element arrays, one per dataset, each of format [prefix, gs:// path to directory with one BED per population]
    Array[Array[String]]? sample_level_comparison_datasets  # Array of two-element arrays, one per dataset, each of format [prefix, gs:// path to per-sample tarballs]
    File? sample_renaming_tsv # File with mapping to rename sample IDs for compatibility with sample_level_comparison_datasets

    File contig_list
    Int? random_seed
    Int? max_gq

    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_qc_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_plot_qc_vcf_wide
    RuntimeAttr? runtime_override_site_level_benchmark_plot
    RuntimeAttr? runtime_override_plot_qc_per_sample
    RuntimeAttr? runtime_override_plot_qc_per_family
    RuntimeAttr? runtime_override_per_sample_benchmark_plot
    RuntimeAttr? runtime_override_sanitize_outputs
    RuntimeAttr? runtime_attr_ids_from_vcf
    RuntimeAttr? runtime_attr_subset_ped

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_subset_vcf
    RuntimeAttr? runtime_override_merge_vcfwide_stat_shards
    RuntimeAttr? runtime_override_merge_vcf_2_bed

    # overrides for CollectQcVcfWide
    RuntimeAttr? runtime_override_preprocess_vcf
    RuntimeAttr? runtime_override_collect_sharded_vcf_stats
    RuntimeAttr? runtime_override_svtk_vcf_2_bed
    RuntimeAttr? runtime_override_scatter_vcf
    RuntimeAttr? runtime_override_merge_subvcf_stat_shards

    # overrides for CollectSiteLevelBenchmarking
    RuntimeAttr? runtime_override_site_level_benchmark
    RuntimeAttr? runtime_override_merge_site_level_benchmark

    # overrides for CollectQcPerSample
    RuntimeAttr? runtime_override_collect_vids_per_sample
    RuntimeAttr? runtime_override_split_samples_list
    RuntimeAttr? runtime_override_tar_shard_vid_lists
    RuntimeAttr? runtime_override_merge_sharded_per_sample_vid_lists

    # overrides for CollectPerSampleBenchmarking
    RuntimeAttr? runtime_override_benchmark_samples
    RuntimeAttr? runtime_override_split_shuffled_list
    RuntimeAttr? runtime_override_merge_and_tar_shard_benchmarks
  }

  Array[String] algorithms = ["manta", "melt", "wham", "depth", "pesr"]
  Array[File?] vcfs_array = [manta_vcf_noOutliers, melt_vcf_noOutliers, wham_vcf_noOutliers, depth_vcf_noOutliers, merged_pesr_vcf]
  Array[File?] vcf_indexes_array = [manta_vcf_noOutliers_index, melt_vcf_noOutliers_index, wham_vcf_noOutliers_index, depth_vcf_noOutliers_index, merged_pesr_vcf_index]
  Int num_algorithms = length(algorithms)

  call util.GetSampleIdsFromVcf {
    input:
      vcf = select_first(vcfs_array),
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_ids_from_vcf
  }

  call util.SubsetPedFile {
    input:
      ped_file = ped_file,
      sample_list = GetSampleIdsFromVcf.out_file,
      subset_name = batch,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_subset_ped
  }

  Int max_gq_ = select_first([max_gq, 999])

  scatter (i in range(num_algorithms)) {
    if (defined(vcfs_array[i]) && defined(vcf_indexes_array[i])) {
      call vcf_qc.MainVcfQc as VcfQc {
        input:
          vcfs = [select_first([vcfs_array[i]])],
          ped_file=SubsetPedFile.ped_subset_file,
          prefix="${batch}.${algorithms[i]}_FilterBatch_filtered_vcf",
          sv_per_shard=2500,
          samples_per_shard=600,
          site_level_comparison_datasets=site_level_comparison_datasets,
          sample_level_comparison_datasets=sample_level_comparison_datasets,
          sample_renaming_tsv=sample_renaming_tsv,
          primary_contigs_fai=contig_list,
          random_seed=random_seed,
          max_gq=max_gq_,
          sv_base_mini_docker=sv_base_mini_docker,
          sv_pipeline_docker=sv_pipeline_docker,
          sv_pipeline_qc_docker=sv_pipeline_qc_docker,
          runtime_override_subset_vcf=runtime_override_subset_vcf,
          runtime_override_preprocess_vcf=runtime_override_preprocess_vcf,
          runtime_override_plot_qc_vcf_wide=runtime_override_plot_qc_vcf_wide,
          runtime_override_site_level_benchmark_plot=runtime_override_site_level_benchmark_plot,
          runtime_override_per_sample_benchmark_plot=runtime_override_per_sample_benchmark_plot,
          runtime_override_plot_qc_per_sample=runtime_override_plot_qc_per_sample,
          runtime_override_plot_qc_per_family=runtime_override_plot_qc_per_family,
          runtime_override_sanitize_outputs=runtime_override_sanitize_outputs,
          runtime_override_merge_vcfwide_stat_shards=runtime_override_merge_vcfwide_stat_shards,
          runtime_override_merge_vcf_2_bed=runtime_override_merge_vcf_2_bed,
          runtime_override_collect_sharded_vcf_stats=runtime_override_collect_sharded_vcf_stats,
          runtime_override_svtk_vcf_2_bed=runtime_override_svtk_vcf_2_bed,
          runtime_override_scatter_vcf=runtime_override_scatter_vcf,
          runtime_override_merge_subvcf_stat_shards=runtime_override_merge_subvcf_stat_shards,
          runtime_override_site_level_benchmark=runtime_override_site_level_benchmark,
          runtime_override_merge_site_level_benchmark=runtime_override_merge_site_level_benchmark,
          runtime_override_collect_vids_per_sample=runtime_override_collect_vids_per_sample,
          runtime_override_split_samples_list=runtime_override_split_samples_list,
          runtime_override_tar_shard_vid_lists=runtime_override_tar_shard_vid_lists,
          runtime_override_benchmark_samples=runtime_override_benchmark_samples,
          runtime_override_split_shuffled_list=runtime_override_split_shuffled_list,
          runtime_override_merge_sharded_per_sample_vid_lists=runtime_override_merge_sharded_per_sample_vid_lists,
          runtime_override_merge_and_tar_shard_benchmarks=runtime_override_merge_and_tar_shard_benchmarks
      }
    }
  }

  output {
    File? filtered_manta_vcf_qc = VcfQc.sv_vcf_qc_output[0]
    File? filtered_melt_vcf_qc = VcfQc.sv_vcf_qc_output[1]
    File? filtered_wham_vcf_qc = VcfQc.sv_vcf_qc_output[2]
    File? filtered_depth_vcf_qc = VcfQc.sv_vcf_qc_output[3]
    File? filtered_pesr_vcf_qc = VcfQc.sv_vcf_qc_output[4]
  }

}
