version 1.0

import "TestUtils.wdl" as tu

workflow Module00aBatchMetrics {
  input {
    Array[String] samples
    String name
    Array[File]? coverage_counts
    Array[File]? pesr_disc
    Array[File]? pesr_split
    Array[File]? delly_vcf
    Array[File]? manta_vcf
    Array[File]? melt_vcf
    Array[File]? wham_vcf

    Array[File]? baseline_delly_vcf
    Array[File]? baseline_manta_vcf
    Array[File]? baseline_melt_vcf
    Array[File]? baseline_wham_vcf

    File contig_list
    File contig_index
    Int min_size = 50
    String sv_pipeline_base_docker
    String linux_docker
  }

  scatter (i in range(length(samples))) {
    if (defined(delly_vcf)) {
      call tu.StandardizeVCF as Delly_Std {
        input:
          vcf = select_first([delly_vcf])[i],
          sample_id = samples[i],
          caller = "delly",
          contig_index = contig_index,
          min_size = min_size,
          sv_pipeline_base_docker = sv_pipeline_base_docker
      }
      if (defined(baseline_delly_vcf)) {
        call tu.StandardizeVCF as Delly_Std_Base {
          input:
            vcf = select_first([baseline_delly_vcf])[i],
            sample_id = samples[i],
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
          samples = [samples[i]],
          prefix = "delly_" + samples[i],
          types = "DEL,DUP,INS,INV,BND",
          contig_list = contig_list,
          sv_pipeline_base_docker = sv_pipeline_base_docker
      }
    }
    if (defined(manta_vcf)) {
      call tu.StandardizeVCF as Manta_Std {
        input:
          vcf = select_first([manta_vcf])[i],
          sample_id = samples[i],
          caller = "manta",
          contig_index = contig_index,
          min_size = min_size,
          sv_pipeline_base_docker = sv_pipeline_base_docker
      }
      if (defined(baseline_manta_vcf)) {
        call tu.StandardizeVCF as Manta_Std_Base {
          input:
            vcf = select_first([baseline_manta_vcf])[i],
            sample_id = samples[i],
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
          samples = [samples[i]],
          prefix = "manta_" + samples[i],
          types = "DEL,DUP,INS,INV,BND",
          contig_list = contig_list,
          sv_pipeline_base_docker = sv_pipeline_base_docker
      }
    }
    if (defined(melt_vcf)) {
      call tu.StandardizeVCF as Melt_Std {
        input:
          vcf = select_first([melt_vcf])[i],
          sample_id = samples[i],
          caller = "melt",
          contig_index = contig_index,
          min_size = min_size,
          sv_pipeline_base_docker = sv_pipeline_base_docker
      }
      if (defined(baseline_melt_vcf)) {
        call tu.StandardizeVCF as Melt_Std_Base {
          input:
            vcf = select_first([baseline_melt_vcf])[i],
            sample_id = samples[i],
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
          samples = [samples[i]],
          prefix = "melt_" + samples[i],
          types = "DEL,DUP,INS,INV,BND",
          contig_list = contig_list,
          sv_pipeline_base_docker = sv_pipeline_base_docker
      }
    }
    if (defined(wham_vcf)) {
      call tu.StandardizeVCF as Wham_Std {
        input:
          vcf = select_first([wham_vcf])[i],
          sample_id = samples[i],
          caller = "wham",
          contig_index = contig_index,
          min_size = min_size,
          sv_pipeline_base_docker = sv_pipeline_base_docker
      }
      if (defined(baseline_wham_vcf)) {
        call tu.StandardizeVCF as Wham_Std_Base {
          input:
            vcf = select_first([baseline_wham_vcf])[i],
            sample_id = samples[i],
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
          samples = [samples[i]],
          prefix = "wham_" + samples[i],
          types = "DEL,DUP,INS,INV,BND",
          contig_list = contig_list,
          sv_pipeline_base_docker = sv_pipeline_base_docker
      }
    }
    if (defined(pesr_split)) {
      call tu.SRMetrics {
        input:
          sr_file = select_first([pesr_split])[i],
          samples = [samples[i]],
          sv_pipeline_base_docker = sv_pipeline_base_docker
      }
    }
    if (defined(pesr_disc)) {
      call tu.PEMetrics {
        input:
          pe_file = select_first([pesr_disc])[i],
          samples = [samples[i]],
          sv_pipeline_base_docker = sv_pipeline_base_docker
      }
    }
    if (defined(coverage_counts)) {
      call tu.CountsMetrics {
        input:
          counts_file = select_first([coverage_counts])[i],
          sample_id = samples[i],
          sv_pipeline_base_docker = sv_pipeline_base_docker
      }
    }
    Array[File] sample_metric_files = select_all([Delly_Metrics.out, Manta_Metrics.out, Melt_Metrics.out, Wham_Metrics.out, SRMetrics.out, PEMetrics.out, CountsMetrics.out])
  }

  call tu.CatMetrics {
    input:
      prefix = "module00a." + name,
      metric_files = flatten(sample_metric_files),
      linux_docker = linux_docker
  }

  output {
    File metrics_file = CatMetrics.out
  }
}
