version 1.0

import "Module02.wdl" as module
import "Module02Metrics.wdl" as metrics
import "TestUtils.wdl" as utils

workflow Module02Test {
  input {
    String test_name
    Array[String] samples
    String base_metrics
  }

  call module.Module02

  call metrics.Module02Metrics {
    input:
      name = test_name,
      metrics = Module02.metrics,
      metrics_common = Module02.metrics_common
  }

  call utils.PlotMetrics {
    input:
      name = test_name,
      samples = samples,
      test_metrics = Module02Metrics.metrics_file,
      base_metrics = base_metrics
  }

  output {
    File metrics = Module02Metrics.metrics_file
    File metrics_plot_pdf = PlotMetrics.metrics_plot_pdf
    File metrics_plot_tsv = PlotMetrics.metrics_plot_tsv
  }
}
