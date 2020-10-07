version 1.0

import "GATKSVPipelineSingleSample.wdl" as module
import "TestUtils.wdl" as utils

workflow GATKSVPipelineSingleSampleTest {
  meta {
    allowNestedInputs: true
  }

  input {
    String test_name
    String case_sample
    Array[String] ref_samples
    String base_metrics
  }

  call module.GATKSVPipelineSingleSample {
    input:
      sample_id = case_sample,
      ref_samples = ref_samples
  }

  Array[String] samples = flatten([[case_sample], ref_samples])

  call utils.PlotMetrics {
    input:
      name = test_name,
      samples = samples,
      test_metrics = GATKSVPipelineSingleSample.metrics_file,
      base_metrics = base_metrics
  }

  output {
    File metrics = GATKSVPipelineSingleSample.metrics_file
    File metrics_plot_pdf = PlotMetrics.metrics_plot_pdf
    File metrics_plot_tsv = PlotMetrics.metrics_plot_tsv
  }
}
