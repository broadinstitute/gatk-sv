version 1.0

import "TestUtils.wdl" as tu
import "Structs.wdl"

workflow SingleSampleMetrics {
  input {
    String name
    Array[String] ref_samples
    String case_sample
    File? wgd_scores
    File? sample_pe
    File? sample_sr
    File? sample_counts

    File? cleaned_vcf
    File? final_vcf
    File? genotyped_pesr_vcf
    File? genotyped_depth_vcf
    File? non_genotyped_unique_depth_calls_vcf

    File? baseline_cleaned_vcf
    File? baseline_final_vcf
    File? baseline_genotyped_pesr_vcf
    File? baseline_genotyped_depth_vcf
    File? baseline_non_genotyped_unique_depth_calls_vcf

    File contig_list
    String linux_docker
    String sv_pipeline_base_docker

    RuntimeAttr? runtime_attr_sr_metrics
    RuntimeAttr? runtime_attr_pe_metrics
    RuntimeAttr? runtime_attr_counts_metrics
    RuntimeAttr? runtime_attr_vcf_metrics
    RuntimeAttr? runtime_attr_pesr_metrics
    RuntimeAttr? runtime_attr_depth_metrics
    RuntimeAttr? runtime_attr_unique_depth_metrics
    RuntimeAttr? runtime_attr_wgd_metrics
    RuntimeAttr? runtime_attr_cat_metrics
  }

  Array[String] samples = flatten([[case_sample], ref_samples])

  if (defined(sample_sr)) {
    call tu.SRMetrics {
      input:
        sr_file = select_first([sample_sr]),
        samples = [case_sample],
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_sr_metrics
    }
  }

  if (defined(sample_pe)) {
    call tu.PEMetrics {
      input:
        pe_file = select_first([sample_pe]),
        samples = [case_sample],
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_pe_metrics
    }
  }

  if (defined(sample_counts)) {
    call tu.CountsMetrics {
      input:
        counts_file = select_first([sample_counts]),
        sample_id = case_sample,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_counts_metrics
    }
  }

  if (defined(cleaned_vcf)) {
    call tu.VCFMetrics as Cleaned_VCF_Metrics {
      input:
        vcf = select_first([cleaned_vcf]),
        baseline_vcf = baseline_cleaned_vcf,
        samples = samples,
        prefix = "cleaned",
        types = "DEL,DUP,INS,INV,CTX,CNV,CPX,BND",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_vcf_metrics
    }
  }

  if (defined(final_vcf)) {
    call tu.VCFMetrics as Final_VCF_Metrics {
      input:
        vcf = select_first([final_vcf]),
        baseline_vcf = baseline_final_vcf,
        samples = samples,
        prefix = "final",
        types = "DEL,DUP,INS,INV,CTX,CNV,CPX,BND",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_vcf_metrics
    }
  }

  if (defined(genotyped_pesr_vcf)) {
    call tu.VCFMetrics as Genotyped_PESR_VCF_Metrics {
      input:
        vcf = select_first([genotyped_pesr_vcf]),
        baseline_vcf = baseline_genotyped_pesr_vcf,
        samples = samples,
        prefix = "genotyped_pesr",
        types = "DEL,DUP,INS,INV,BND",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_pesr_metrics
    }
  }

  if (defined(genotyped_depth_vcf)) {
    call tu.VCFMetrics as Genotyped_Depth_VCF_Metrics {
      input:
        vcf = select_first([genotyped_depth_vcf]),
        baseline_vcf = baseline_genotyped_depth_vcf,
        samples = samples,
        prefix = "genotyped_depth",
        types = "DEL,DUP",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_depth_metrics
    }
  }

  if (defined(non_genotyped_unique_depth_calls_vcf)) {
    call tu.VCFMetrics as NonGenotypedUniqueDepthCallsVCFMetrics {
      input:
        vcf = select_first([non_genotyped_unique_depth_calls_vcf]),
        baseline_vcf = baseline_non_genotyped_unique_depth_calls_vcf,
        samples = samples,
        prefix = "non_genotyped_uniq_depth",
        types = "DEL,DUP",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_unique_depth_metrics
    }
  }

  if (defined(wgd_scores)) {
    call SingleSampleWGDMetrics {
      input:
        wgd_scores = select_first([wgd_scores]),
        linux_docker = linux_docker,
        runtime_attr_override = runtime_attr_wgd_metrics
    }
  }

  call tu.CatMetrics {
    input:
      prefix = "single_sample." + name,
      metric_files = select_all([SingleSampleWGDMetrics.out, SRMetrics.out, PEMetrics.out, CountsMetrics.out, Genotyped_PESR_VCF_Metrics.out, Genotyped_Depth_VCF_Metrics.out, Cleaned_VCF_Metrics.out, Final_VCF_Metrics.out, NonGenotypedUniqueDepthCallsVCFMetrics.out]),
      search_string = case_sample,
      replace_string = "sample",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_cat_metrics
  }

  output {
    File metrics_file = CatMetrics.out
  }
}

task SingleSampleWGDMetrics {
  input {
    File wgd_scores
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr runtime_default = object {
    mem_gb: 1.0,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  output {
    File out = "wgd.tsv"
  }
  command <<<

    SCORE=$(zcat ~{wgd_scores} | tail -n1 | cut -f2)
    echo "wgd_score_sample	${SCORE}" > wgd.tsv

  >>>
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: linux_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }
}
