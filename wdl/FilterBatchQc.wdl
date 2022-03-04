version 1.0

import "MasterVcfQc.wdl" as vcf_qc
import "Utils.wdl" as util

workflow FilterBatchQc {
  input {
    File? manta_vcf_noOutliers
    File? delly_vcf_noOutliers
    File? melt_vcf_noOutliers
    File? wham_vcf_noOutliers
    File? depth_vcf_noOutliers
    File? merged_pesr_vcf
    String batch
    File ped_file
    Array[File]? thousand_genomes_benchmark_calls
    Array[File]? hgsv_benchmark_calls
    Array[File]? asc_benchmark_calls
    File? sanders_2015_tarball
    File? collins_2017_tarball
    File? werling_2018_tarball

    File contig_list
    Int? random_seed

    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_qc_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_plot_qc_vcf_wide
    RuntimeAttr? runtime_override_thousand_g_benchmark
    RuntimeAttr? runtime_override_thousand_g_plot
    RuntimeAttr? runtime_override_asc_benchmark
    RuntimeAttr? runtime_override_asc_plot
    RuntimeAttr? runtime_override_custom_external
    RuntimeAttr? runtime_override_hgsv_benchmark
    RuntimeAttr? runtime_override_hgsv_plot
    RuntimeAttr? runtime_override_plot_qc_per_sample
    RuntimeAttr? runtime_override_plot_qc_per_family
    RuntimeAttr? runtime_override_sanders_per_sample_plot
    RuntimeAttr? runtime_override_collins_per_sample_plot
    RuntimeAttr? runtime_override_werling_per_sample_plot
    RuntimeAttr? runtime_override_sanitize_outputs
    RuntimeAttr? runtime_attr_ids_from_vcf
    RuntimeAttr? runtime_attr_subset_ped

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_merge_vcfwide_stat_shards
    RuntimeAttr? runtime_override_merge_vcf_2_bed

    # overrides for ShardedQcCollection
    RuntimeAttr? runtime_override_collect_sharded_vcf_stats
    RuntimeAttr? runtime_override_svtk_vcf_2_bed
    RuntimeAttr? runtime_override_split_vcf_to_qc
    RuntimeAttr? runtime_override_merge_subvcf_stat_shards
    RuntimeAttr? runtime_override_merge_svtk_vcf_2_bed

    # overrides for CollectQcPerSample
    RuntimeAttr? runtime_override_collect_vids_per_sample
    RuntimeAttr? runtime_override_split_samples_list
    RuntimeAttr? runtime_override_tar_shard_vid_lists

    # overrides for PerSampleExternalBenchmark
    RuntimeAttr? runtime_override_benchmark_samples
    RuntimeAttr? runtime_override_split_shuffled_list
    RuntimeAttr? runtime_override_merge_and_tar_shard_benchmarks
  }

  Array[String] algorithms = ["manta", "delly", "melt", "wham", "depth", "pesr"]
  Array[File?] vcfs_array = [manta_vcf_noOutliers, delly_vcf_noOutliers, melt_vcf_noOutliers, wham_vcf_noOutliers, depth_vcf_noOutliers, merged_pesr_vcf]
  Int num_algorithms = length(algorithms)

  Array[String] contigs = transpose(read_tsv(contig_list))[0]

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

  scatter (i in range(num_algorithms)) {
    if (defined(vcfs_array[i])) {
      call vcf_qc.MasterVcfQc as VcfQc {
        input:
          vcf = select_first([vcfs_array[i]]),
          ped_file=SubsetPedFile.ped_subset_file,
          prefix="${batch}.${algorithms[i]}_FilterBatch_filtered_vcf",
          sv_per_shard=10000,
          samples_per_shard=100,
          thousand_genomes_benchmark_calls=thousand_genomes_benchmark_calls,
          hgsv_benchmark_calls=hgsv_benchmark_calls,
          asc_benchmark_calls=asc_benchmark_calls,
          sanders_2015_tarball=sanders_2015_tarball,
          collins_2017_tarball=collins_2017_tarball,
          werling_2018_tarball=werling_2018_tarball,
          contigs=contigs,
          random_seed=random_seed,
          sv_base_mini_docker=sv_base_mini_docker,
          sv_pipeline_docker=sv_pipeline_docker,
          sv_pipeline_qc_docker=sv_pipeline_qc_docker,
          runtime_override_plot_qc_vcf_wide=runtime_override_plot_qc_vcf_wide,
          runtime_override_thousand_g_benchmark=runtime_override_thousand_g_benchmark,
          runtime_override_thousand_g_plot=runtime_override_thousand_g_plot,
          runtime_override_asc_benchmark=runtime_override_asc_benchmark,
          runtime_override_asc_plot=runtime_override_asc_plot,
          runtime_override_custom_external=runtime_override_custom_external,
          runtime_override_hgsv_benchmark=runtime_override_hgsv_benchmark,
          runtime_override_hgsv_plot=runtime_override_hgsv_plot,
          runtime_override_plot_qc_per_sample=runtime_override_plot_qc_per_sample,
          runtime_override_plot_qc_per_family=runtime_override_plot_qc_per_family,
          runtime_override_sanders_per_sample_plot=runtime_override_sanders_per_sample_plot,
          runtime_override_collins_per_sample_plot=runtime_override_collins_per_sample_plot,
          runtime_override_werling_per_sample_plot=runtime_override_werling_per_sample_plot,
          runtime_override_sanitize_outputs=runtime_override_sanitize_outputs,
          runtime_override_merge_vcfwide_stat_shards=runtime_override_merge_vcfwide_stat_shards,
          runtime_override_merge_vcf_2_bed=runtime_override_merge_vcf_2_bed,
          runtime_override_collect_sharded_vcf_stats=runtime_override_collect_sharded_vcf_stats,
          runtime_override_svtk_vcf_2_bed=runtime_override_svtk_vcf_2_bed,
          runtime_override_split_vcf_to_qc=runtime_override_split_vcf_to_qc,
          runtime_override_merge_subvcf_stat_shards=runtime_override_merge_subvcf_stat_shards,
          runtime_override_merge_svtk_vcf_2_bed=runtime_override_merge_svtk_vcf_2_bed,
          runtime_override_collect_vids_per_sample=runtime_override_collect_vids_per_sample,
          runtime_override_split_samples_list=runtime_override_split_samples_list,
          runtime_override_tar_shard_vid_lists=runtime_override_tar_shard_vid_lists,
          runtime_override_benchmark_samples=runtime_override_benchmark_samples,
          runtime_override_split_shuffled_list=runtime_override_split_shuffled_list,
          runtime_override_merge_and_tar_shard_benchmarks=runtime_override_merge_and_tar_shard_benchmarks
      }
    }
  }

  output {
    File? filtered_manta_vcf_qc = VcfQc.sv_vcf_qc_output[0]
    File? filtered_delly_vcf_qc = VcfQc.sv_vcf_qc_output[1]
    File? filtered_melt_vcf_qc = VcfQc.sv_vcf_qc_output[2]
    File? filtered_wham_vcf_qc = VcfQc.sv_vcf_qc_output[3]
    File? filtered_depth_vcf_qc = VcfQc.sv_vcf_qc_output[4]
    File? filtered_pesr_vcf_qc = VcfQc.sv_vcf_qc_output[5]
  }

}
