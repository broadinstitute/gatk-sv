version 1.0

import "TestUtils.wdl" as tu

workflow SingleSampleMetrics {
  input {
    String name
    Array[String] ref_samples
    String case_sample
    File wgd_scores
    File sample_pe
    File sample_sr
    File sample_counts

    File cleaned_vcf
    File final_vcf
    File genotyped_pesr_vcf
    File genotyped_depth_vcf
    File non_genotyped_unique_depth_calls_vcf

    File? baseline_cleaned_vcf
    File? baseline_final_vcf
    File? baseline_genotyped_pesr_vcf
    File? baseline_genotyped_depth_vcf
    File? baseline_non_genotyped_unique_depth_calls_vcf

    File contig_list
    String linux_docker
    String sv_pipeline_base_docker
  }

  Array[String] samples = flatten([[case_sample], ref_samples])

  call tu.SRMetrics {
    input:
      sr_file = sample_sr,
      samples = [case_sample],
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call tu.PEMetrics {
    input:
      pe_file = sample_pe,
      samples = [case_sample],
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call tu.CountsMetrics {
    input:
      counts_file = sample_counts,
      sample_id = case_sample,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call tu.VCFMetrics as Cleaned_VCF_Metrics {
    input:
      vcf = cleaned_vcf,
      baseline_vcf = baseline_cleaned_vcf,
      samples = samples,
      prefix = "cleaned",
      types = "DEL,DUP,INS,INV,CTX,CNV,CPX,BND",
      contig_list = contig_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call tu.VCFMetrics as Final_VCF_Metrics {
    input:
      vcf = final_vcf,
      baseline_vcf = baseline_final_vcf,
      samples = samples,
      prefix = "final",
      types = "DEL,DUP,INS,INV,CTX,CNV,CPX,BND",
      contig_list = contig_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call tu.VCFMetrics as Genotyped_PESR_VCF_Metrics {
    input:
      vcf = genotyped_pesr_vcf,
      baseline_vcf = baseline_genotyped_pesr_vcf,
      samples = samples,
      prefix = "genotyped_pesr",
      types = "DEL,DUP,INS,INV,BND",
      contig_list = contig_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call tu.VCFMetrics as Genotyped_Depth_VCF_Metrics {
    input:
      vcf = genotyped_depth_vcf,
      baseline_vcf = baseline_genotyped_depth_vcf,
      samples = samples,
      prefix = "genotyped_depth",
      types = "DEL,DUP",
      contig_list = contig_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call tu.VCFMetrics as NonGenotypedUniqueDepthCallsVCFMetrics {
    input:
      vcf = non_genotyped_unique_depth_calls_vcf,
      baseline_vcf = baseline_non_genotyped_unique_depth_calls_vcf,
      samples = samples,
      prefix = "non_genotyped_uniq_depth",
      types = "DEL,DUP",
      contig_list = contig_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call SingleSampleWGDMetrics {
    input:
      wgd_scores = wgd_scores,
      linux_docker = linux_docker
  }

  call tu.CatMetrics {
    input:
      prefix = "single_sample." + name,
      metric_files = [SingleSampleWGDMetrics.out, SRMetrics.out, PEMetrics.out, CountsMetrics.out, Genotyped_PESR_VCF_Metrics.out, Genotyped_Depth_VCF_Metrics.out, Cleaned_VCF_Metrics.out, Final_VCF_Metrics.out, NonGenotypedUniqueDepthCallsVCFMetrics.out],
      search_string = case_sample,
      replace_string = "sample",
      linux_docker = linux_docker
  }

  output {
    File metrics_file = CatMetrics.out
  }
}

task SingleSampleWGDMetrics {
  input {
    File wgd_scores
    String linux_docker
    Float mem_gib = 1
    Int disk_gb = 10
    Int preemptible_attempts = 3
  }

  output {
    File out = "wgd.tsv"
  }
  command <<<

    SCORE=$(zcat ~{wgd_scores} | tail -n1 | cut -f2)
    echo "wgd_score_sample	${SCORE}" > wgd.tsv

  >>>
  runtime {
    cpu: 1
    memory: "~{mem_gib} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    bootDiskSizeGb: 10
    docker: linux_docker
    preemptible: preemptible_attempts
    maxRetries: 1
  }

}
