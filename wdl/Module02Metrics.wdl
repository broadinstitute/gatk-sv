version 1.0

import "TestUtils.wdl" as tu

workflow Module02Metrics {
  input {
    String name
    File metrics
    File metrics_common
    File contig_list
    String sv_pipeline_base_docker
    String linux_docker
  }

  call tu.MetricsFileMetrics {
    input:
      metrics_file = metrics,
      contig_list = contig_list,
      common = false,
      prefix = "module02.non_common." + name,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call tu.MetricsFileMetrics as CommonMetricsFileMetrics {
    input:
      metrics_file = metrics_common,
      contig_list = contig_list,
      common = true,
      prefix = "module02.common." + name,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  call tu.CatMetrics {
    input:
      prefix = "module02." + name,
      metric_files = [MetricsFileMetrics.out, CommonMetricsFileMetrics.out],
      linux_docker = linux_docker
  }

  output {
    File metrics_file = CatMetrics.out
  }
}
