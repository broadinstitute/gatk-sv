version 1.0

import "TestUtils.wdl" as tu

workflow Module0506Metrics {
  input {
    Array[String] samples
    String name

    File? cluster_vcf
    File? complex_resolve_vcf
    File? complex_genotype_vcf
    File cleaned_vcf

    File? baseline_cluster_vcf
    File? baseline_complex_resolve_vcf
    File? baseline_complex_genotype_vcf
    File? baseline_cleaned_vcf

    File contig_list
    String linux_docker
    String sv_pipeline_base_docker

    # Do not use
    Array[File]? NONE_FILE_ARRAY_
  }

  if (defined(cluster_vcf)) {
    call tu.VCFMetrics as ClusterMetrics {
      input:
        vcf = select_first([cluster_vcf]),
        baseline_vcf = baseline_cluster_vcf,
        samples = samples,
        prefix = "cluster",
        types = "DEL,DUP,INS,INV,CTX,CPX,BND",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
  }

  if (defined(complex_resolve_vcf)) {
    call tu.VCFMetrics as ComplexResolveMetrics {
      input:
        vcf = select_first([complex_resolve_vcf]),
        baseline_vcf = baseline_complex_resolve_vcf,
        samples = samples,
        prefix = "cpx_resolve",
        types = "DEL,DUP,INS,INV,CTX,CPX,BND",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
  }

  if (defined(complex_genotype_vcf)) {
    call tu.VCFMetrics as ComplesGenotypeMetrics {
      input:
        vcf = select_first([complex_genotype_vcf]),
        baseline_vcf = baseline_complex_genotype_vcf,
        samples = samples,
        prefix = "cpx_genotype",
        types = "DEL,DUP,INS,INV,CTX,CPX,BND",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
  }

  call tu.VCFMetrics as CleanedMetrics {
    input:
      vcf = cleaned_vcf,
      baseline_vcf = baseline_cleaned_vcf,
      samples = samples,
      prefix = "cleaned",
      types = "DEL,DUP,INS,INV,CTX,CPX,BND",
      contig_list = contig_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call tu.CatMetrics {
    input:
      prefix = "module0506." + name,
      metric_files = select_all([ClusterMetrics.out, ComplexResolveMetrics.out, ComplesGenotypeMetrics.out, CleanedMetrics.out]),
      linux_docker = linux_docker
  }

  output {
    File metrics_file = CatMetrics.out
  }
}
