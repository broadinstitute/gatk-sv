version 1.0

import "Module05_06.wdl" as module
import "Module05_06Metrics.wdl" as metrics
import "TestUtils.wdl" as utils

workflow Module05_06Test {
  input {
    String test_name
    Array[String] samples
    String base_metrics
  }

  call module.Module05_06

  call metrics.Module05_06Metrics {
    input:
      name = test_name,
      samples = samples,
      final_vcf = Module05_06.vcf_cpx,
      cleaned_vcf = Module05_06.vcf
  }

  call utils.PlotMetrics {
    input:
      name = test_name,
      samples = samples,
      test_metrics = Module05_06Metrics.metrics_file,
      base_metrics = base_metrics
  }

  output {
    File metrics = Module05_06Metrics.metrics_file
    File metrics_plot_pdf = PlotMetrics.metrics_plot_pdf
    File metrics_plot_tsv = PlotMetrics.metrics_plot_tsv
  }
}
