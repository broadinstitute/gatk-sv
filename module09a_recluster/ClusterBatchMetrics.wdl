version 1.0

import "TestUtils.wdl" as tu
import "Utils.wdl" as util

workflow ClusterBatchMetrics {
  input {
    Array[String]? samples
    String name

    File depth_vcf
    File? manta_vcf
    File? wham_vcf
    File? melt_vcf
    File? scramble_vcf

    File? baseline_depth_vcf
    File? baseline_manta_vcf
    File? baseline_wham_vcf
    File? baseline_melt_vcf
    File? baseline_scramble_vcf

    File contig_list
    String? sv_base_mini_docker  # required if not providing samples array
    String sv_pipeline_base_docker
    String linux_docker

    RuntimeAttr? runtime_attr_sample_ids_from_vcf
    RuntimeAttr? runtime_attr_depth_metrics
    RuntimeAttr? runtime_attr_manta_metrics
    RuntimeAttr? runtime_attr_melt_metrics
    RuntimeAttr? runtime_attr_scramble_metrics
    RuntimeAttr? runtime_attr_wham_metrics
    RuntimeAttr? runtime_attr_cat_metrics
  }

  if (!defined(samples)) {
    call util.GetSampleIdsFromVcf {
      input:
        vcf = depth_vcf,
        sv_base_mini_docker = select_first([sv_base_mini_docker]),
        runtime_attr_override = runtime_attr_sample_ids_from_vcf
    }
  }

  call tu.VCFMetrics as depth_metrics {
    input:
      vcf = depth_vcf,
      baseline_vcf = baseline_depth_vcf,
      samples = select_first([samples, GetSampleIdsFromVcf.out_array]),
      prefix = "depth_clustered",
      types = "DEL,DUP",
      contig_list = contig_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_depth_metrics
  }
  if (defined(manta_vcf)) {
    call tu.VCFMetrics as manta_metrics {
      input:
        vcf = select_first([manta_vcf]),
        baseline_vcf = baseline_manta_vcf,
        samples = select_first([samples, GetSampleIdsFromVcf.out_array]),
        prefix = "manta_clustered",
        types = "DEL,DUP,INS,INV,BND",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_manta_metrics
    }
  }
  if (defined(melt_vcf)) {
    call tu.VCFMetrics as melt_metrics {
      input:
        vcf = select_first([melt_vcf]),
        baseline_vcf = baseline_melt_vcf,
        samples = select_first([samples, GetSampleIdsFromVcf.out_array]),
        prefix = "melt_clustered",
        types = "INS",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_melt_metrics
    }
  }
  if (defined(scramble_vcf)) {
    call tu.VCFMetrics as scramble_metrics {
      input:
        vcf = select_first([scramble_vcf]),
        baseline_vcf = baseline_scramble_vcf,
        samples = select_first([samples, GetSampleIdsFromVcf.out_array]),
        prefix = "scramble_clustered",
        types = "INS",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_scramble_metrics
    }
  }
  if (defined(wham_vcf)) {
    call tu.VCFMetrics as wham_metrics {
      input:
        vcf = select_first([wham_vcf]),
        baseline_vcf = baseline_wham_vcf,
        samples = select_first([samples, GetSampleIdsFromVcf.out_array]),
        prefix = "wham_clustered",
        types = "DEL,DUP",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_wham_metrics
    }
  }

  call tu.CatMetrics {
    input:
      prefix = "ClusterBatch." + name,
      metric_files = select_all([depth_metrics.out, manta_metrics.out, melt_metrics.out, scramble_metrics.out, wham_metrics.out]),
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_cat_metrics
  }

  output {
    File metrics_file = CatMetrics.out
  }
}
