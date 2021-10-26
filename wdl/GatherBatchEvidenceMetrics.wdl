version 1.0

import "TestUtils.wdl" as tu

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
    Array[File]? std_delly_vcf
    Array[File]? std_manta_vcf
    Array[File]? std_melt_vcf
    Array[File]? std_scramble_vcf
    Array[File]? std_wham_vcf

    File? baseline_merged_dels
    File? baseline_merged_dups
    File? baseline_median_cov
    Array[File]? baseline_std_delly_vcf
    Array[File]? baseline_std_manta_vcf
    Array[File]? baseline_std_melt_vcf
    Array[File]? baseline_std_scramble_vcf
    Array[File]? baseline_std_wham_vcf

    File contig_list
    String sv_pipeline_base_docker
    String linux_docker

    RuntimeAttr? runtime_attr_delly_metrics
    RuntimeAttr? runtime_attr_manta_metrics
    RuntimeAttr? runtime_attr_melt_metrics
    RuntimeAttr? runtime_attr_scramble_metrics
    RuntimeAttr? runtime_attr_wham_metrics
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

  scatter (i in range(length(samples))) {
    File? baseline_delly_vcf_i = if defined(baseline_std_delly_vcf) then select_first([baseline_std_delly_vcf])[i] else FILE_NONE_
    File? baseline_manta_vcf_i = if defined(baseline_std_manta_vcf) then select_first([baseline_std_manta_vcf])[i] else FILE_NONE_
    File? baseline_melt_vcf_i = if defined(baseline_std_melt_vcf) then select_first([baseline_std_melt_vcf])[i] else FILE_NONE_
    File? baseline_scramble_vcf_i = if defined(baseline_std_scramble_vcf) then select_first([baseline_std_scramble_vcf])[i] else FILE_NONE_
    File? baseline_wham_vcf_i = if defined(baseline_std_wham_vcf) then select_first([baseline_std_wham_vcf])[i] else FILE_NONE_
    if (defined(std_delly_vcf)) {
      call tu.VCFMetrics as Delly_Std_Metrics {
        input:
          vcf = select_first([std_delly_vcf])[i],
          baseline_vcf = baseline_delly_vcf_i,
          samples = [samples[i]],
          prefix = "delly_std_" + samples[i],
          types = "DEL,DUP,INS,INV,BND",
          contig_list = contig_list,
          sv_pipeline_base_docker = sv_pipeline_base_docker,
          runtime_attr_override = runtime_attr_delly_metrics
      }
    }
    if (defined(std_manta_vcf)) {
      call tu.VCFMetrics as Manta_Std_Metrics {
        input:
          vcf = select_first([std_manta_vcf])[i],
          baseline_vcf = baseline_manta_vcf_i,
          samples = [samples[i]],
          prefix = "manta_std_" + samples[i],
          types = "DEL,DUP,INS,INV,BND",
          contig_list = contig_list,
          sv_pipeline_base_docker = sv_pipeline_base_docker,
          runtime_attr_override = runtime_attr_manta_metrics
      }
    }
    if (defined(std_melt_vcf)) {
      call tu.VCFMetrics as Melt_Std_Metrics {
        input:
          vcf = select_first([std_melt_vcf])[i],
          baseline_vcf = baseline_melt_vcf_i,
          samples = [samples[i]],
          prefix = "melt_std_" + samples[i],
          types = "INS",
          contig_list = contig_list,
          sv_pipeline_base_docker = sv_pipeline_base_docker,
          runtime_attr_override = runtime_attr_melt_metrics
      }
    }
    if (defined(std_scramble_vcf)) {
      call tu.VCFMetrics as Scramble_Std_Metrics {
        input:
          vcf = select_first([std_scramble_vcf])[i],
          baseline_vcf = baseline_scramble_vcf_i,
          samples = [samples[i]],
          prefix = "scramble_std_" + samples[i],
          types = "INS",
          contig_list = contig_list,
          sv_pipeline_base_docker = sv_pipeline_base_docker,
          runtime_attr_override = runtime_attr_scramble_metrics
      }
    }
    if (defined(std_wham_vcf)) {
      call tu.VCFMetrics as Wham_Std_Metrics {
        input:
          vcf = select_first([std_wham_vcf])[i],
          baseline_vcf = baseline_wham_vcf_i,
          samples = [samples[i]],
          prefix = "wham_std_" + samples[i],
          types = "DEL,DUP",
          contig_list = contig_list,
          sv_pipeline_base_docker = sv_pipeline_base_docker,
          runtime_attr_override = runtime_attr_wham_metrics
      }
    }
    Array[File] sample_metric_files_ = select_all([Delly_Std_Metrics.out, Manta_Std_Metrics.out, Melt_Std_Metrics.out, Scramble_Std_Metrics.out, Wham_Std_Metrics.out])
  }
  Array[File] sample_metric_files = flatten(sample_metric_files_)

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
      metric_files = flatten([sample_metric_files, [BAFMetrics.out, SRMetrics.out, PEMetrics.out, BincovMetrics.out, MedcovMetrics.out, del_metrics, dup_metrics]]),
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_cat_metrics
  }

  output {
    File metrics_file = CatMetrics.out
  }
}
