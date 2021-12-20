version 1.0

import "FilterBatchSites.wdl" as filter_sites
import "PlotSVCountsPerSample.wdl" as sv_counts
import "FilterBatchSamples.wdl" as filter_outliers
import "Utils.wdl" as util
import "FilterBatchMetrics.wdl" as metrics

workflow FilterBatch {
  input {
    String batch
    File? manta_vcf
    File? delly_vcf
    File? wham_vcf
    File? melt_vcf
    File? depth_vcf
    File evidence_metrics
    File evidence_metrics_common

    Int outlier_cutoff_nIQR
    File? outlier_cutoff_table

    # Module metrics parameters
    # Run module metrics workflow at the end - on by default
    Boolean? run_module_metrics
    String? sv_pipeline_base_docker  # required if run_module_metrics = true
    File? primary_contigs_list  # required if run_module_metrics = true
    File? ped_file  # required if run_module_metrics = true
    File? baseline_filtered_depth_vcf  # baseline files are optional for metrics workflow
    File? baseline_filtered_pesr_vcf

    File original_scores
    File original_cutoffs
    File outlier_samples_list

    String sv_pipeline_docker
    String sv_base_mini_docker
    String linux_docker

    RuntimeAttr? runtime_attr_adjudicate
    RuntimeAttr? runtime_attr_rewrite_scores
    RuntimeAttr? runtime_attr_filter_annotate_vcf
    RuntimeAttr? runtime_attr_ids_from_vcf
    RuntimeAttr? runtime_attr_merge_pesr_vcfs

    RuntimeAttr? runtime_attr_count_svs
    RuntimeAttr? runtime_attr_plot_svcounts
    RuntimeAttr? runtime_attr_cat_outliers_preview

    RuntimeAttr? runtime_attr_identify_outliers
    RuntimeAttr? runtime_attr_exclude_outliers
    RuntimeAttr? runtime_attr_cat_outliers
    RuntimeAttr? runtime_attr_filter_samples
  }

  call filter_sites.FilterBatchSites {
    input:
      batch = batch,
      manta_vcf = manta_vcf,
      melt_vcf = melt_vcf,
      delly_vcf = delly_vcf,
      depth_vcf = depth_vcf,
      wham_vcf = wham_vcf,
      evidence_metrics = evidence_metrics,
      evidence_metrics_common = evidence_metrics_common,
      original_cutoffs = original_cutoffs,
      original_scores = original_scores,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_adjudicate = runtime_attr_adjudicate,
      runtime_attr_rewrite_scores = runtime_attr_rewrite_scores,
      runtime_attr_filter_annotate_vcf = runtime_attr_filter_annotate_vcf
  }

  call filter_outliers.FilterBatchSamples {
    input:
      batch = batch,
      outlier_cutoff_table = outlier_cutoff_table,
      N_IQR_cutoff = outlier_cutoff_nIQR,
      manta_vcf = FilterBatchSites.sites_filtered_manta_vcf,
      delly_vcf = FilterBatchSites.sites_filtered_delly_vcf,
      wham_vcf = FilterBatchSites.sites_filtered_wham_vcf,
      melt_vcf = FilterBatchSites.sites_filtered_melt_vcf,
      depth_vcf = FilterBatchSites.sites_filtered_depth_vcf,
      outlier_samples_list = outlier_samples_list,
      linux_docker = linux_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_identify_outliers = runtime_attr_identify_outliers,
      runtime_attr_exclude_outliers= runtime_attr_exclude_outliers,
      runtime_attr_cat_outliers = runtime_attr_cat_outliers,
      runtime_attr_filter_samples = runtime_attr_filter_samples,
      runtime_attr_ids_from_vcf = runtime_attr_ids_from_vcf,
      runtime_attr_merge_pesr_vcfs = runtime_attr_merge_pesr_vcfs,
      runtime_attr_count_svs = runtime_attr_count_svs
  }

  call sv_counts.PlotSVCountsPerSample {
    input:
      prefix = batch,
      vcfs = [FilterBatchSamples.outlier_filtered_manta_vcf, FilterBatchSamples.outlier_filtered_delly_vcf, FilterBatchSamples.outlier_filtered_wham_vcf, FilterBatchSamples.outlier_filtered_melt_vcf, FilterBatchSamples.outlier_filtered_depth_vcf],
      N_IQR_cutoff = outlier_cutoff_nIQR,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_count_svs = runtime_attr_count_svs,
      runtime_attr_plot_svcounts = runtime_attr_plot_svcounts,
      runtime_attr_cat_outliers_preview = runtime_attr_cat_outliers_preview
  }

  output {
    File? filtered_manta_vcf = FilterBatchSamples.outlier_filtered_manta_vcf
    File? filtered_delly_vcf = FilterBatchSamples.outlier_filtered_delly_vcf
    File? filtered_wham_vcf = FilterBatchSamples.outlier_filtered_wham_vcf
    File? filtered_melt_vcf = FilterBatchSamples.outlier_filtered_melt_vcf
    File? filtered_depth_vcf = FilterBatchSamples.outlier_filtered_depth_vcf
    File? filtered_pesr_vcf = FilterBatchSamples.outlier_filtered_pesr_vcf
    Array[File] sv_counts = PlotSVCountsPerSample.sv_counts
    Array[File] sv_count_plots = PlotSVCountsPerSample.sv_count_plots
    Array[String] batch_samples_postOutlierExclusion = FilterBatchSamples.filtered_batch_samples_list
    File batch_samples_postOutlierExclusion_file = FilterBatchSamples.filtered_batch_samples_file
    File cutoffs = FilterBatchSites.cutoffs
    File scores = FilterBatchSites.scores
    File RF_intermediate_files = FilterBatchSites.RF_intermediate_files
  }
}

