version 1.0

import "Module00aBatch.wdl" as module
import "Module00aBatchMetrics.wdl" as metrics
import "TestUtils.wdl" as utils

workflow Module00aBatchTest {
  input {
    String test_name
    String base_metrics
    Array[String] samples

    Array[File]? NONE_ARR_
  }

  call module.Module00aBatch {
    input:
      sample_ids = samples
  }

  Array[File]? coverage_counts = if length(select_all(Module00aBatch.coverage_counts)) > 0 then select_all(Module00aBatch.coverage_counts) else NONE_ARR_
  Array[File]? pesr_disc = if length(select_all(Module00aBatch.pesr_disc)) > 0 then select_all(Module00aBatch.pesr_disc) else NONE_ARR_
  Array[File]? pesr_split = if length(select_all(Module00aBatch.pesr_split)) > 0 then select_all(Module00aBatch.pesr_split) else NONE_ARR_
  Array[File]? delly_vcf = if length(select_all(Module00aBatch.delly_vcf)) > 0 then select_all(Module00aBatch.delly_vcf) else NONE_ARR_
  Array[File]? manta_vcf = if length(select_all(Module00aBatch.manta_vcf)) > 0 then select_all(Module00aBatch.manta_vcf) else NONE_ARR_
  Array[File]? melt_vcf = if length(select_all(Module00aBatch.melt_vcf)) > 0 then select_all(Module00aBatch.melt_vcf) else NONE_ARR_
  Array[File]? wham_vcf = if length(select_all(Module00aBatch.wham_vcf)) > 0 then select_all(Module00aBatch.wham_vcf) else NONE_ARR_

  call metrics.Module00aBatchMetrics {
    input:
      name = test_name,
      samples = samples,
      coverage_counts = coverage_counts,
      pesr_disc = pesr_disc,
      pesr_split = pesr_split,
      delly_vcf = delly_vcf,
      manta_vcf = manta_vcf,
      melt_vcf = melt_vcf,
      wham_vcf = wham_vcf
  }

  call utils.PlotMetrics {
    input:
      name = test_name,
      samples = samples,
      test_metrics = Module00aBatchMetrics.metrics_file,
      base_metrics = base_metrics
  }

  output {
    File metrics = Module00aBatchMetrics.metrics_file
    File metrics_plot_pdf = PlotMetrics.metrics_plot_pdf
    File metrics_plot_tsv = PlotMetrics.metrics_plot_tsv
  }
}
