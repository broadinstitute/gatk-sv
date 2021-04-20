version 1.0

import "Module03.wdl" as module
import "Module03Metrics.wdl" as metrics
import "TestUtils.wdl" as utils

workflow Module03Test {
  meta {
    allowNestedInputs: true
  }

  input {
    String test_name
    Array[String] samples
    String base_metrics
  }

  call module.Module03

  call metrics.Module03Metrics {
    input:
      name = test_name,
      samples = samples,
      filtered_pesr_vcf = Module03.filtered_pesr_vcf,
      filtered_depth_vcf = select_first([Module03.filtered_depth_vcf]),
      cutoffs = Module03.cutoffs,
      outlier_list = Module03.outlier_samples_excluded_file,
      samples_post_filtering_file = Module03.batch_samples_postOutlierExclusion_file
  }

  call utils.PlotMetrics {
    input:
      name = test_name,
      samples = samples,
      test_metrics = Module03Metrics.metrics_file,
      base_metrics = base_metrics
  }

  output {
    File metrics = Module03Metrics.metrics_file
    File metrics_plot_pdf = PlotMetrics.metrics_plot_pdf
    File metrics_plot_tsv = PlotMetrics.metrics_plot_tsv
  }
}
