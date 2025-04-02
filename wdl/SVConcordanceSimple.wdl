version 1.0

import "SVConcordance.wdl" as conc
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow SVConcordanceSimple {
  input {
    File eval_vcf
    File eval_vcf_idx
    File truth_vcf
    File truth_vcf_idx

    String output_prefix

    Float? pesr_interval_overlap = 0.5
    Float? pesr_size_similarity = 0.0
    Int? pesr_breakend_window = 5000

    Float? depth_interval_overlap = 0.5
    Float? depth_size_similarity = 0.0
    Int? depth_breakend_window = 10000000

    Float? mixed_interval_overlap = 0.5
    Float? mixed_size_similarity = 0.0
    Int? mixed_breakend_window = 10000000

    File? clustering_config
    File? stratification_config
    Array[String]? track_names
    Array[File]? track_intervals
    Int? stratification_num_breakpoint_overlaps = 0
    Float? stratification_overlap_fraction = 0.5

    File reference_dict

    String sv_base_mini_docker
    String gatk_docker

    RuntimeAttr? runtime_attr_sv_concordance
    RuntimeAttr? runtime_attr_sort_vcf
  }

  call conc.SVConcordanceTask {
      input:
        truth_vcf=truth_vcf,
        eval_vcf=eval_vcf,
        output_prefix="~{output_prefix}.unsorted",
        additional_args="--pesr-interval-overlap ~{pesr_interval_overlap} --pesr-size-similarity ~{pesr_size_similarity} --pesr-breakend-window ~{pesr_breakend_window} --depth-interval-overlap ~{depth_interval_overlap} --depth-size-similarity ~{depth_size_similarity} --depth-breakend-window ~{depth_breakend_window} --mixed-interval-overlap ~{mixed_interval_overlap} --mixed-size-similarity ~{mixed_size_similarity} --mixed-breakend-window ~{mixed_breakend_window} --stratification-num-breakpoint-overlaps ~{stratification_num_breakpoint_overlaps} --stratification-overlap-fraction ~{stratification_overlap_fraction}",
        clustering_config=clustering_config,
        stratification_config=stratification_config,
        track_names=track_names,
        track_intervals=track_intervals,
        reference_dict=reference_dict,
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_sv_concordance
  }

  call tasks_cohort.SortVcf {
      input:
        vcf=SVConcordanceTask.out_unsorted,
        outfile_prefix="~{output_prefix}.sorted",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_sort_vcf
    }

  output {
    File conc_vcf = SortVcf.out
    File conc_vcf_idx = SortVcf.out_index
  }
}