version 1.0

import "TestUtils.wdl" as tu

workflow Module03Metrics {
  input {
    Array[String] samples
    String name

    File filtered_pesr_vcf
    File filtered_depth_vcf
    File cutoffs
    File outlier_list
    File filtered_ped_file
    File samples_post_filtering_file

    File? baseline_filtered_pesr_vcf
    File? baseline_filtered_depth_vcf

    File contig_list
    String linux_docker
    String sv_pipeline_base_docker
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
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call tu.VCFMetrics as Depth_VCF_Metrics {
    input:
      vcf = filtered_depth_vcf,
      baseline_vcf = baseline_filtered_depth_vcf,
      samples = samples_post_filtering,
      prefix = "filtered_depth",
      types = "DEL,DUP",
      contig_list = contig_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call tu.CutoffAndOutlierMetrics {
    input:
      cutoffs = cutoffs,
      outlier_list = outlier_list,
      filtered_ped_file = filtered_ped_file,
      samples = samples,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call tu.CatMetrics {
    input:
      prefix = "module03." + name,
      metric_files = [PESR_VCF_Metrics.out, Depth_VCF_Metrics.out, CutoffAndOutlierMetrics.out],
      linux_docker = linux_docker
  }

  output {
    File metrics_file = CatMetrics.out
  }
}
