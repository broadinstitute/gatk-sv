version 1.0

import "Module01.wdl" as module
import "Module01Metrics.wdl" as metrics
import "TestUtils.wdl" as utils

workflow Module01Test {
  input {
    String test_name
    Array[String] samples
    String base_metrics
  }

  call module.Module01

  call metrics.Module01Metrics {
    input:
      name = test_name,
      samples = samples,
      depth_vcf = Module01.depth_vcf,
      delly_vcf = Module01.delly_vcf,
      manta_vcf = Module01.manta_vcf,
      wham_vcf = Module01.wham_vcf,
      melt_vcf = Module01.melt_vcf
  }

  call utils.PlotMetrics {
    input:
      name = test_name,
      samples = samples,
      test_metrics = Module01Metrics.metrics_file,
      base_metrics = base_metrics
  }

  output {
    File metrics = Module01Metrics.metrics_file
    File metrics_plot_pdf = PlotMetrics.metrics_plot_pdf
    File metrics_plot_tsv = PlotMetrics.metrics_plot_tsv
  }
}
