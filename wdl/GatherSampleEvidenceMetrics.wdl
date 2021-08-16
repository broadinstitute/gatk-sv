version 1.0

import "TestUtils.wdl" as tu

workflow GatherSampleEvidenceMetrics {
  input {
    String sample
    File? coverage_counts
    File? pesr_disc
    File? pesr_split
    File? delly_vcf
    File? manta_vcf
    File? melt_vcf
    File? wham_vcf

    File? baseline_delly_vcf
    File? baseline_manta_vcf
    File? baseline_melt_vcf
    File? baseline_wham_vcf

    File contig_list
    File contig_index
    Int min_size = 50
    String sv_pipeline_base_docker
  }

  if (defined(delly_vcf)) {
    call tu.StandardizeVCF as Delly_Std {
      input:
        vcf = select_first([delly_vcf]),
        sample_id = sample,
        caller = "delly",
        contig_index = contig_index,
        min_size = min_size,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
    if (defined(baseline_delly_vcf)) {
      call tu.StandardizeVCF as Delly_Std_Base {
        input:
          vcf = select_first([baseline_delly_vcf]),
          sample_id = sample,
          caller = "delly",
          contig_index = contig_index,
          min_size = min_size,
          sv_pipeline_base_docker = sv_pipeline_base_docker
      }
    }
    call tu.VCFMetrics as Delly_Metrics {
      input:
        vcf = Delly_Std.out,
        baseline_vcf = Delly_Std_Base.out,
        samples = [sample],
        prefix = "delly_" + sample,
        types = "DEL,DUP,INS,INV,BND",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
  }
  if (defined(manta_vcf)) {
    call tu.StandardizeVCF as Manta_Std {
      input:
        vcf = select_first([manta_vcf]),
        sample_id = sample,
        caller = "manta",
        contig_index = contig_index,
        min_size = min_size,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
    if (defined(baseline_manta_vcf)) {
      call tu.StandardizeVCF as Manta_Std_Base {
        input:
          vcf = select_first([baseline_manta_vcf]),
          sample_id = sample,
          caller = "manta",
          contig_index = contig_index,
          min_size = min_size,
          sv_pipeline_base_docker = sv_pipeline_base_docker
      }
    }
    call tu.VCFMetrics as Manta_Metrics {
      input:
        vcf = Manta_Std.out,
        baseline_vcf = Manta_Std_Base.out,
        samples = [sample],
        prefix = "manta_" + sample,
        types = "DEL,DUP,INS,INV,BND",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
  }
  if (defined(melt_vcf)) {
    call tu.StandardizeVCF as Melt_Std {
      input:
        vcf = select_first([melt_vcf]),
        sample_id = sample,
        caller = "melt",
        contig_index = contig_index,
        min_size = min_size,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
    if (defined(baseline_melt_vcf)) {
      call tu.StandardizeVCF as Melt_Std_Base {
        input:
          vcf = select_first([baseline_melt_vcf]),
          sample_id = sample,
          caller = "melt",
          contig_index = contig_index,
          min_size = min_size,
          sv_pipeline_base_docker = sv_pipeline_base_docker
      }
    }
    call tu.VCFMetrics as Melt_Metrics {
      input:
        vcf = Melt_Std.out,
        baseline_vcf = Melt_Std_Base.out,
        samples = [sample],
        prefix = "melt_" + sample,
        types = "DEL,DUP,INS,INV,BND",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
  }
  if (defined(wham_vcf)) {
    call tu.StandardizeVCF as Wham_Std {
      input:
        vcf = select_first([wham_vcf]),
        sample_id = sample,
        caller = "wham",
        contig_index = contig_index,
        min_size = min_size,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
    if (defined(baseline_wham_vcf)) {
      call tu.StandardizeVCF as Wham_Std_Base {
        input:
          vcf = select_first([baseline_wham_vcf]),
          sample_id = sample,
          caller = "wham",
          contig_index = contig_index,
          min_size = min_size,
          sv_pipeline_base_docker = sv_pipeline_base_docker
      }
    }
    call tu.VCFMetrics as Wham_Metrics {
      input:
        vcf = Wham_Std.out,
        baseline_vcf = Wham_Std_Base.out,
        samples = [sample],
        prefix = "wham_" + sample,
        types = "DEL,DUP,INS,INV,BND",
        contig_list = contig_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
  }
  if (defined(pesr_split)) {
    call tu.SRMetrics {
      input:
        sr_file = select_first([pesr_split]),
        samples = [sample],
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
  }
  if (defined(pesr_disc)) {
    call tu.PEMetrics {
      input:
        pe_file = select_first([pesr_disc]),
        samples = [sample],
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
  }
  if (defined(coverage_counts)) {
    call tu.CountsMetrics {
      input:
        counts_file = select_first([coverage_counts]),
        sample_id = sample,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
  }

  output {
    Array[File] sample_metrics_files = select_all([Delly_Metrics.out, Manta_Metrics.out, Melt_Metrics.out, Wham_Metrics.out, SRMetrics.out, PEMetrics.out, CountsMetrics.out])
  }
}
