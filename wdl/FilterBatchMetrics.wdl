version 1.0

import "TestUtils.wdl" as tu
import "Utils.wdl" as util
import "Structs.wdl"

workflow FilterBatchMetrics {
  input {
    Array[String] samples
    String name

    File filtered_pesr_vcf
    File filtered_depth_vcf
    File cutoffs
    File outlier_list
    File ped_file
    File samples_post_filtering_file

    File? baseline_filtered_pesr_vcf
    File? baseline_filtered_depth_vcf

    File contig_list
    String linux_docker
    String sv_pipeline_base_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_subset_ped
    RuntimeAttr? runtime_attr_pesr_vcf_metrics
    RuntimeAttr? runtime_attr_depth_vcf_metrics
    RuntimeAttr? runtime_attr_cutoff_outlier_metrics
    RuntimeAttr? runtime_attr_cat_metrics
  }

  Array[String] samples_post_filtering = read_lines(samples_post_filtering_file)

  call tu.VCFMetrics as PESR_VCF_Metrics {
    input:
      vcf = filtered_pesr_vcf,
      baseline_vcf = baseline_filtered_pesr_vcf,
      samples = samples_post_filtering,
      prefix = "filtered_pesr",
      types = "DEL,DUP,INS,INV,BND",
      contig_list = contig_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_pesr_vcf_metrics
  }

  call tu.VCFMetrics as Depth_VCF_Metrics {
    input:
      vcf = filtered_depth_vcf,
      baseline_vcf = baseline_filtered_depth_vcf,
      samples = samples_post_filtering,
      prefix = "filtered_depth",
      types = "DEL,DUP",
      contig_list = contig_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_depth_vcf_metrics
  }

  call util.SubsetPedFile {
    input:
      ped_file = ped_file,
      sample_list = samples_post_filtering_file,
      subset_name = name,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_subset_ped
  }

  call tu.CutoffAndOutlierMetrics {
    input:
      cutoffs = cutoffs,
      outlier_list = outlier_list,
      filtered_ped_file = SubsetPedFile.ped_subset_file,
      samples = samples,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_cutoff_outlier_metrics
  }

  call tu.CatMetrics {
    input:
      prefix = "FilterBatch." + name,
      metric_files = [PESR_VCF_Metrics.out, Depth_VCF_Metrics.out, CutoffAndOutlierMetrics.out],
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_cat_metrics
  }

  output {
    File metrics_file = CatMetrics.out
  }
}
