version 1.0

import "TestUtils.wdl" as tu

workflow Module05_06Metrics {
  input {
    Array[String] samples
    String name

    File final_vcf
    File cleaned_vcf

    File? baseline_final_vcf
    File? baseline_cleaned_vcf

    File contig_list
    String linux_docker
    String sv_pipeline_base_docker
  }

  call tu.VCFMetrics as Final_VCF_Metrics {
    input:
      vcf = final_vcf,
      baseline_vcf = baseline_final_vcf,
      samples = samples,
      prefix = "regenotyped",
      types = "DEL,DUP,INS,INV,CTX,CPX,BND",
      contig_list = contig_list,
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

  call tu.CatMetrics {
    input:
      prefix = "module05_06." + name,
      metric_files = [Final_VCF_Metrics.out, Cleaned_VCF_Metrics.out],
      linux_docker = linux_docker
  }

  output {
    File metrics_file = CatMetrics.out
  }
}
