version 1.0

import "TestUtils.wdl" as tu
import "Utils.wdl" as util

workflow GatherBatchEvidenceMetrics {
  input {
    Array[String] samples
    String name
    File merged_BAF
    File merged_SR
    File merged_PE
    File merged_bincov
    File merged_dels
    File merged_dups
    File median_cov

    File? baseline_merged_dels
    File? baseline_merged_dups
    File? baseline_median_cov

    File contig_list
    String sv_pipeline_base_docker
    String linux_docker

    RuntimeAttr? runtime_attr_baf_metrics
    RuntimeAttr? runtime_attr_sr_metrics
    RuntimeAttr? runtime_attr_pe_metrics
    RuntimeAttr? runtime_attr_bincov_metrics
    RuntimeAttr? runtime_attr_medcov_metrics
    RuntimeAttr? runtime_attr_merged_del
    RuntimeAttr? runtime_attr_merged_dup
    RuntimeAttr? runtime_attr_cat_metrics

    # Do not use
    File? FILE_NONE_
  }

  call tu.BAFMetrics {
    input:
      baf_file = merged_BAF,
      samples = samples,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = select_first([runtime_attr_baf_metrics, {"mem_gb": 30, "disk_gb": 100}])
  }
  call tu.SRMetrics {
    input:
      sr_file = merged_SR,
      samples = samples,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = select_first([runtime_attr_sr_metrics, {"mem_gb": 30, "disk_gb": 100, "preemptible_tries": 0}])
  }
  call tu.PEMetrics {
    input:
      pe_file = merged_PE,
      samples = samples,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = select_first([runtime_attr_pe_metrics, {"mem_gb": 15, "disk_gb": 100, "preemptible_tries": 0}])
  }
  call tu.BincovMetrics {
    input:
      bincov_matrix = merged_bincov,
      samples = samples,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_bincov_metrics
  }
  call tu.MedcovMetrics {
    input:
      medcov_file = median_cov,
      samples = samples,
      baseline_medcov_file = baseline_median_cov,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_medcov_metrics
  }
  if (defined(baseline_merged_dels)) {
    call tu.MergedDepthMetricsWithBaseline as MergedDelWithBaseline {
      input:
        bed = merged_dels,
        baseline_bed = select_first([baseline_merged_dels]),
        type = "DEL",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_merged_del
    }
  }
  if (!defined(baseline_merged_dels)) { # else
    call tu.MergedDepthMetricsWithoutBaseline as MergedDel {
      input:
        bed = merged_dels,
        type = "DEL",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_merged_del
    }
  }
  if (defined(baseline_merged_dups)) {
    call tu.MergedDepthMetricsWithBaseline as MergedDupWithBaseline {
      input:
        bed = merged_dups,
        baseline_bed = select_first([baseline_merged_dups]),
        type = "DUP",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_merged_dup
    }
  }
  if (!defined(baseline_merged_dups)) { # else
    call tu.MergedDepthMetricsWithoutBaseline as MergedDup {
      input:
        bed = merged_dups,
        type = "DUP",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_merged_dup
    }
  }

  File del_metrics = select_first([MergedDelWithBaseline.out, MergedDel.out])
  File dup_metrics = select_first([MergedDupWithBaseline.out, MergedDup.out])

  call tu.CatMetrics {
    input:
      prefix = "GatherBatchEvidence." + name,
      metric_files = [BAFMetrics.out, SRMetrics.out, PEMetrics.out, BincovMetrics.out, MedcovMetrics.out, del_metrics, dup_metrics],
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_cat_metrics
  }

  output {
    File metrics_file = CatMetrics.out
  }
}
