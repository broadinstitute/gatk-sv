version 1.0

import "TestUtils.wdl" as tu
import "Utils.wdl" as util

workflow MakeCohortVcfMetrics {
  input {
    File? samples_list
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
    String sv_pipeline_docker
    String? sv_base_mini_docker  # required if not providing samples array

    RuntimeAttr? runtime_attr_sample_ids_from_vcf
    RuntimeAttr? runtime_attr_vcf_metrics
    RuntimeAttr? runtime_attr_cat_metrics

    # Do not use
    Array[File]? NONE_FILE_ARRAY_
  }

  if (!defined(samples_list)) {
    call util.GetSampleIdsFromVcf {
      input:
        vcf = cleaned_vcf,
        sv_base_mini_docker = select_first([sv_base_mini_docker]),
        runtime_attr_override = runtime_attr_sample_ids_from_vcf
    }
  }
  Array[String] samples = read_lines(select_first([samples_list, GetSampleIdsFromVcf.out_file]))

  if (defined(cluster_vcf)) {
    call tu.VCFMetrics as ClusterMetrics {
      input:
        vcf = select_first([cluster_vcf]),
        baseline_vcf = baseline_cluster_vcf,
        samples = samples,
        prefix = "cluster",
        types = "DEL,DUP,INS,INV,CTX,CPX,BND",
        contig_list = contig_list,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_vcf_metrics
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
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_vcf_metrics
    }
  }

  if (defined(complex_genotype_vcf)) {
    call tu.VCFMetrics as ComplexGenotypeMetrics {
      input:
        vcf = select_first([complex_genotype_vcf]),
        baseline_vcf = baseline_complex_genotype_vcf,
        samples = samples,
        prefix = "cpx_genotype",
        types = "DEL,DUP,INS,INV,CTX,CPX,BND",
        contig_list = contig_list,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_vcf_metrics
    }
  }

  call tu.VCFMetrics as CleanedMetrics {
    input:
      vcf = cleaned_vcf,
      baseline_vcf = baseline_cleaned_vcf,
      samples = samples,
      prefix = "cleaned",
      types = "DEL,DUP,INS,INV,CTX,CPX,BND,CNV",
      contig_list = contig_list,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_vcf_metrics
  }

  call tu.CatMetrics {
    input:
      prefix = "MakeCohortVcf." + name,
      metric_files = select_all([ClusterMetrics.out, ComplexResolveMetrics.out, ComplexGenotypeMetrics.out, CleanedMetrics.out]),
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_cat_metrics
  }

  output {
    File metrics_file = CatMetrics.out
  }
}
