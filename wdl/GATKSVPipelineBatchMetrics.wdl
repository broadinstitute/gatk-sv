version 1.0

import "Module00aBatchMetrics.wdl" as m00a
import "Module00cMetrics.wdl" as m00c
import "Module01Metrics.wdl" as m01
import "Module02Metrics.wdl" as m02
import "Module03Metrics.wdl" as m03
import "Module04Metrics.wdl" as m04
import "Module0506Metrics.wdl" as m0506
import "TestUtils.wdl" as utils

workflow BatchMetrics {
  input {
    String name
    Array[String] samples
    File contig_list
    File contig_index
    String linux_docker
    String sv_pipeline_base_docker

    File? baseline_00a_metrics
    File? baseline_00c_metrics
    File? baseline_01_metrics
    File? baseline_02_metrics
    File? baseline_03_metrics
    File? baseline_04_metrics
    File? baseline_0506_metrics

    # 00a
    Array[File] coverage_counts
    Array[File] pesr_disc
    Array[File] pesr_split
    Array[File]? delly_vcf
    Array[File]? manta_vcf
    Array[File]? melt_vcf
    Array[File]? wham_vcf

    Array[File]? baseline_delly_vcf
    Array[File]? baseline_manta_vcf
    Array[File]? baseline_melt_vcf
    Array[File]? baseline_wham_vcf

    # 00c
    File merged_BAF
    File merged_SR
    File merged_PE
    File merged_bincov
    File merged_dels
    File merged_dups
    File median_cov
    Array[File]? std_delly_vcf
    Array[File]? std_manta_vcf
    Array[File]? std_melt_vcf
    Array[File]? std_wham_vcf

    File? baseline_merged_dels
    File? baseline_merged_dups
    File? baseline_median_cov
    Array[File]? baseline_std_delly_vcf
    Array[File]? baseline_std_manta_vcf
    Array[File]? baseline_std_melt_vcf
    Array[File]? baseline_std_wham_vcf

    # 01
    File merged_depth_vcf
    File? merged_delly_vcf
    File? merged_manta_vcf
    File? merged_wham_vcf
    File? merged_melt_vcf

    File? baseline_merged_depth_vcf
    File? baseline_merged_delly_vcf
    File? baseline_merged_manta_vcf
    File? baseline_merged_wham_vcf
    File? baseline_merged_melt_vcf

    # 02
    File metrics
    File metrics_common

    # 03
    File filtered_pesr_vcf
    File filtered_depth_vcf
    File cutoffs
    File outlier_list
    File filtered_ped_file
    File samples_post_filtering_file

    File? baseline_filtered_pesr_vcf
    File? baseline_filtered_depth_vcf

    # 04
    File genotyped_pesr_vcf
    File genotyped_depth_vcf
    File cutoffs_pesr_pesr
    File cutoffs_pesr_depth
    File cutoffs_depth_pesr
    File cutoffs_depth_depth
    File sr_bothside_pass
    File sr_background_fail

    File? baseline_genotyped_pesr_vcf
    File? baseline_genotyped_depth_vcf

    # 0506
    File? module0506_cluster_vcf
    File? module0506_complex_resolve_vcf
    File? module0506_complex_genotype_vcf
    File module0506_cleaned_vcf

    File? baseline_module0506_cluster_vcf
    File? baseline_module0506_complex_resolve_vcf
    File? baseline_module0506_complex_genotype_vcf
    File? baseline_module0506_cleaned_vcf
  }

  Array[String] samples_post_filter = read_lines(samples_post_filtering_file)

  call m00a.Module00aBatchMetrics {
    input:
      name = name,
      samples = samples,
      coverage_counts = coverage_counts,
      pesr_disc = pesr_disc,
      pesr_split = pesr_split,
      delly_vcf = delly_vcf,
      manta_vcf = manta_vcf,
      melt_vcf = melt_vcf,
      wham_vcf = wham_vcf,
      baseline_delly_vcf = baseline_delly_vcf,
      baseline_manta_vcf = baseline_manta_vcf,
      baseline_melt_vcf = baseline_melt_vcf,
      baseline_wham_vcf = baseline_wham_vcf,
      contig_list = contig_list,
      contig_index = contig_index,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      linux_docker = linux_docker
  }

  call m00c.Module00cMetrics {
    input:
      name = name,
      samples = samples,
      merged_BAF = merged_BAF,
      merged_SR = merged_SR,
      merged_PE = merged_PE,
      merged_bincov = merged_bincov,
      merged_dels = merged_dels,
      merged_dups = merged_dups,
      median_cov = median_cov,
      std_delly_vcf = std_delly_vcf,
      std_manta_vcf = std_manta_vcf,
      std_melt_vcf = std_melt_vcf,
      std_wham_vcf = std_wham_vcf,
      baseline_merged_dels = baseline_merged_dels,
      baseline_merged_dups = baseline_merged_dups,
      baseline_median_cov = baseline_median_cov,
      baseline_std_delly_vcf = baseline_std_delly_vcf,
      baseline_std_manta_vcf = baseline_std_manta_vcf,
      baseline_std_melt_vcf = baseline_std_melt_vcf,
      baseline_std_wham_vcf = baseline_std_wham_vcf,
      contig_list = contig_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      linux_docker = linux_docker
  }

  call m01.Module01Metrics {
    input:
      name = name,
      samples = samples,
      depth_vcf = merged_depth_vcf,
      delly_vcf = merged_delly_vcf,
      manta_vcf = merged_manta_vcf,
      wham_vcf = merged_wham_vcf,
      melt_vcf = merged_melt_vcf,
      baseline_depth_vcf = baseline_merged_depth_vcf,
      baseline_delly_vcf = baseline_merged_delly_vcf,
      baseline_manta_vcf = baseline_merged_manta_vcf,
      baseline_wham_vcf = baseline_merged_wham_vcf,
      baseline_melt_vcf = baseline_merged_melt_vcf,
      contig_list = contig_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      linux_docker = linux_docker
  }

  call m02.Module02Metrics {
    input:
      name = name,
      metrics = metrics,
      metrics_common = metrics_common,
      contig_list = contig_list,
      linux_docker = linux_docker,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call m03.Module03Metrics {
    input:
      name = name,
      samples = samples,
      filtered_pesr_vcf = filtered_pesr_vcf,
      filtered_depth_vcf = filtered_depth_vcf,
      cutoffs = cutoffs,
      outlier_list = outlier_list,
      filtered_ped_file = filtered_ped_file,
      samples_post_filtering_file = samples_post_filtering_file,
      baseline_filtered_pesr_vcf = baseline_filtered_pesr_vcf,
      baseline_filtered_depth_vcf = baseline_filtered_depth_vcf,
      contig_list = contig_list,
      linux_docker = linux_docker,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call m04.Module04Metrics {
    input:
      name = name,
      samples = samples_post_filter,
      genotyped_pesr_vcf = genotyped_pesr_vcf,
      genotyped_depth_vcf = genotyped_depth_vcf,
      cutoffs_pesr_pesr = cutoffs_pesr_pesr,
      cutoffs_pesr_depth = cutoffs_pesr_depth,
      cutoffs_depth_pesr = cutoffs_depth_pesr,
      cutoffs_depth_depth = cutoffs_depth_depth,
      sr_bothside_pass = sr_bothside_pass,
      sr_background_fail = sr_background_fail,
      baseline_genotyped_pesr_vcf = baseline_genotyped_pesr_vcf,
      baseline_genotyped_depth_vcf = baseline_genotyped_depth_vcf,
      contig_list = contig_list,
      linux_docker = linux_docker,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call m0506.Module0506Metrics {
    input:
      name = name,
      samples = samples_post_filter,
      cluster_vcf = module0506_cluster_vcf,
      complex_resolve_vcf = module0506_complex_resolve_vcf,
      complex_genotype_vcf = module0506_complex_genotype_vcf,
      cleaned_vcf = module0506_cleaned_vcf,
      baseline_cluster_vcf = baseline_module0506_cluster_vcf,
      baseline_complex_resolve_vcf = baseline_module0506_complex_resolve_vcf,
      baseline_complex_genotype_vcf = baseline_module0506_complex_genotype_vcf,
      baseline_cleaned_vcf = baseline_module0506_cleaned_vcf,
      contig_list = contig_list,
      linux_docker = linux_docker,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call utils.CatMetrics as CatBatchMetrics {
    input:
      prefix = "batch_sv." + name,
      metric_files = [Module00aBatchMetrics.metrics_file, Module00cMetrics.metrics_file, Module01Metrics.metrics_file, Module02Metrics.metrics_file, Module03Metrics.metrics_file, Module04Metrics.metrics_file, Module0506Metrics.metrics_file],
      linux_docker = linux_docker
  }

  Array[File] defined_baseline_metrics = select_all([baseline_00a_metrics, baseline_00c_metrics, baseline_01_metrics, baseline_02_metrics, baseline_03_metrics, baseline_04_metrics, baseline_0506_metrics])
  if (length(defined_baseline_metrics) > 0) {
    call utils.CatMetrics as CatBaselineMetrics {
      input:
        prefix = "baseline." + name,
        metric_files = defined_baseline_metrics,
        linux_docker = linux_docker
    }
    call utils.PlotMetrics {
      input:
        name = name,
        samples = samples,
        test_metrics = CatBatchMetrics.out,
        base_metrics = CatBaselineMetrics.out,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
  }

  output {
    File metrics_file = CatBatchMetrics.out
    File? metrics_plot_pdf = PlotMetrics.metrics_plot_pdf
    File? metrics_plot_tsv = PlotMetrics.metrics_plot_tsv
  }
}
