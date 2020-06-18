version 1.0

import "Module00c.wdl" as module
import "Module00cMetrics.wdl" as metrics
import "TestUtils.wdl" as utils

workflow Module00cTest {
  input {
    String test_name
    Array[String] samples
    String base_metrics
  }

  call module.Module00c {
    input:
      samples = samples
  }

  call metrics.Module00cMetrics {
    input:
      name = test_name,
      samples = samples,
      merged_BAF = select_first([Module00c.merged_BAF]),
      merged_SR = Module00c.merged_SR,
      merged_PE = Module00c.merged_PE,
      merged_bincov = Module00c.merged_bincov,
      merged_dels = Module00c.merged_dels,
      merged_dups = Module00c.merged_dups,
      median_cov = Module00c.median_cov,
      std_delly_vcf = Module00c.std_delly_vcf,
      std_manta_vcf = Module00c.std_manta_vcf,
      std_melt_vcf = Module00c.std_melt_vcf,
      std_wham_vcf = Module00c.std_wham_vcf
  }

  call utils.PlotMetrics {
    input:
      name = test_name,
      samples = samples,
      test_metrics = Module00cMetrics.metrics_file,
      base_metrics = base_metrics
  }

  output {
    File metrics = Module00cMetrics.metrics_file
    File metrics_plot_pdf = PlotMetrics.metrics_plot_pdf
    File metrics_plot_tsv = PlotMetrics.metrics_plot_tsv
  }
}
