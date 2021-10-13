version 1.0

import "TestUtils.wdl" as tu

workflow GenotypeBatchMetrics {
  input {
    Array[String] samples
    String name

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

    File contig_list
    String linux_docker
    String sv_pipeline_base_docker

    RuntimeAttr? runtime_attr_pesr_metrics
    RuntimeAttr? runtime_attr_depth_metrics
    RuntimeAttr? runtime_attr_cutoff_metrics
    RuntimeAttr? runtime_attr_id_list_metrics
    RuntimeAttr? runtime_attr_cat_metrics
  }

  call tu.VCFMetrics as PESR_VCF_Metrics {
    input:
      vcf = genotyped_pesr_vcf,
      baseline_vcf = baseline_genotyped_pesr_vcf,
      samples = samples,
      prefix = "genotyped_pesr",
      types = "DEL,DUP,INS,INV,BND",
      contig_list = contig_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_pesr_metrics
  }

  call tu.VCFMetrics as Depth_VCF_Metrics {
    input:
      vcf = genotyped_depth_vcf,
      baseline_vcf = baseline_genotyped_depth_vcf,
      samples = samples,
      prefix = "genotyped_depth",
      types = "DEL,DUP",
      contig_list = contig_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_depth_metrics
  }

  call tu.GenotypingCutoffMetrics as Cutoff_PESR_PESR {
    input:
      cutoffs = cutoffs_pesr_pesr,
      name = "pesr_pesr",
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_cutoff_metrics
  }

  call tu.GenotypingCutoffMetrics as Cutoff_PESR_Depth {
    input:
      cutoffs = cutoffs_pesr_depth,
      name = "pesr_depth",
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_cutoff_metrics
  }

  call tu.GenotypingCutoffMetrics as Cutoff_Depth_PESR {
    input:
      cutoffs = cutoffs_depth_pesr,
      name = "depth_pesr",
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_cutoff_metrics
  }

  call tu.GenotypingCutoffMetrics as Cutoff_Depth_Depth {
    input:
      cutoffs = cutoffs_depth_depth,
      name = "depth_depth",
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_cutoff_metrics
  }

  call tu.IdListMetrics as Background_Fail {
    input:
      id_list = sr_background_fail,
      name = "sr_background_fail",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_id_list_metrics
  }

  call tu.IdListMetrics as Bothside_Pass {
    input:
      id_list = sr_bothside_pass,
      name = "sr_bothside_pass",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_id_list_metrics
  }

  call tu.CatMetrics {
    input:
      prefix = "GenotypeBatch." + name,
      metric_files = [PESR_VCF_Metrics.out, Depth_VCF_Metrics.out, Cutoff_PESR_PESR.out, Cutoff_PESR_Depth.out, Cutoff_Depth_PESR.out, Cutoff_Depth_Depth.out, Background_Fail.out, Bothside_Pass.out],
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_cat_metrics
  }

  output {
    File metrics_file = CatMetrics.out
  }
}
