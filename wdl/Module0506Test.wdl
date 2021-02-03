version 1.0

import "Module0506.wdl" as module
import "Module0506Metrics.wdl" as metrics
import "TestUtils.wdl" as utils

workflow Module0506Test {
  meta {
    allowNestedInputs: true
  }

  input {
    String test_name
    Array[String] samples
    String base_metrics
  }

  call module.Module0506 {
    input:
      merge_cluster_vcfs = true,
      merge_complex_resolve_vcfs = true,
      merge_complex_genotype_vcfs = true
  }

  call metrics.Module0506Metrics {
    input:
      name = test_name,
      samples = samples,
      cluster_vcf = Module0506.cluster_vcf,
      complex_resolve_vcf = Module0506.complex_resolve_vcf,
      complex_genotype_vcf = Module0506.complex_genotype_vcf,
      cleaned_vcf = Module0506.vcf
  }

  call utils.PlotMetrics {
    input:
      name = test_name,
      samples = samples,
      test_metrics = Module0506Metrics.metrics_file,
      base_metrics = base_metrics
  }

  output {
    File metrics = Module0506Metrics.metrics_file
    File metrics_plot_pdf = PlotMetrics.metrics_plot_pdf
    File metrics_plot_tsv = PlotMetrics.metrics_plot_tsv
  }
}
