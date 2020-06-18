version 1.0

import "TestUtils.wdl" as tu

workflow Module01Metrics {
  input {
    Array[String] samples
    String name

    File depth_vcf
    File? delly_vcf
    File? manta_vcf
    File? wham_vcf
    File? melt_vcf

    File? baseline_depth_vcf
    File? baseline_delly_vcf
    File? baseline_manta_vcf
    File? baseline_wham_vcf
    File? baseline_melt_vcf

    File contig_list
    String sv_pipeline_base_docker
    String linux_docker
  }

  call tu.VCFMetrics as depth_metrics {
    input:
      vcf = depth_vcf,
      baseline_vcf = baseline_depth_vcf,
      samples = samples,
      prefix = "depth_clustered",
      types = "DEL,DUP",
      contig_list = contig_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }
  if (defined(delly_vcf)) {
    call tu.VCFMetrics as delly_metrics {
      input:
        vcf = select_first([delly_vcf]),
        baseline_vcf = baseline_delly_vcf,
        samples = samples,
        prefix = "delly_clustered",
        types = "DEL,DUP,INS,INV,BND",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
  }
  if (defined(manta_vcf)) {
    call tu.VCFMetrics as manta_metrics {
      input:
        vcf = select_first([manta_vcf]),
        baseline_vcf = baseline_manta_vcf,
        samples = samples,
        prefix = "manta_clustered",
        types = "DEL,DUP,INS,INV,BND",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
  }
  if (defined(melt_vcf)) {
    call tu.VCFMetrics as melt_metrics {
      input:
        vcf = select_first([melt_vcf]),
        baseline_vcf = baseline_melt_vcf,
        samples = samples,
        prefix = "melt_clustered",
        types = "INS",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
  }
  if (defined(wham_vcf)) {
    call tu.VCFMetrics as wham_metrics {
      input:
        vcf = select_first([wham_vcf]),
        baseline_vcf = baseline_wham_vcf,
        samples = samples,
        prefix = "wham_clustered",
        types = "DEL,DUP",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
  }

  call tu.CatMetrics {
    input:
      prefix = "module01." + name,
      metric_files = select_all([depth_metrics.out, delly_metrics.out, manta_metrics.out, melt_metrics.out, wham_metrics.out]),
      linux_docker = linux_docker
  }

  output {
    File metrics_file = CatMetrics.out
  }
}
