version 1.0

import "SVConcordance.wdl" as conc
import "TasksMakeCohortVcf.wdl" as tasks_cohort
import "Utils.wdl" as utils

workflow SVConcordanceSimplePerSample {
  input {
    File eval_vcf
    File eval_vcf_idx
    File truth_vcf
    File truth_vcf_idx

    Float? pesr_interval_overlap = 0.5
    Float? pesr_size_similarity = 0.0
    Int? pesr_breakend_window = 500

    Float? depth_interval_overlap = 0.8
    Float? depth_size_similarity = 0.0
    Int? depth_breakend_window = 10000000

    Float? mixed_interval_overlap = 0.8
    Float? mixed_size_similarity = 0.0
    Int? mixed_breakend_window = 1000

    File? clustering_config
    File? stratification_config
    Array[String]? track_names
    Array[File]? track_intervals
    Int? stratify_num_breakpoint_overlaps = 0
    Float? stratify_overlap_fraction = 0.5

    File reference_dict

    Array[String] sample_ids
    String output_prefix

    String sv_base_mini_docker
    String gatk_docker

    RuntimeAttr? runtime_attr_sv_concordance
    RuntimeAttr? runtime_attr_sort_vcf
  }

  scatter (sample in sample_ids) {
    call utils.SubsetVcfToSample as SubsetEval {
      input:
        vcf = eval_vcf,
        vcf_idx = eval_vcf_idx,
        sample = sample,
        remove_sample = false,
        sv_base_mini_docker = sv_base_mini_docker
    }

    call utils.SubsetVcfToSample as SubsetTruth {
      input:
        vcf = truth_vcf,
        vcf_idx = truth_vcf_idx,
        sample = sample,
        remove_sample = false,
        sv_base_mini_docker = sv_base_mini_docker
    }

    call conc.SVConcordanceTask as Concordance {
      input:
        truth_vcf = SubsetTruth.vcf_subset,
        eval_vcf = SubsetEval.vcf_subset,
        output_prefix = "~{output_prefix}.unsorted",
        additional_args = "--pesr-interval-overlap ~{pesr_interval_overlap} --pesr-size-similarity ~{pesr_size_similarity} --pesr-breakend-window ~{pesr_breakend_window} --depth-interval-overlap ~{depth_interval_overlap} --depth-size-similarity ~{depth_size_similarity} --depth-breakend-window ~{depth_breakend_window} --mixed-interval-overlap ~{mixed_interval_overlap} --mixed-size-similarity ~{mixed_size_similarity} --mixed-breakend-window ~{mixed_breakend_window} --stratify-num-breakpoint-overlaps ~{stratify_num_breakpoint_overlaps} --stratify-overlap-fraction ~{stratify_overlap_fraction}",
        clustering_config = clustering_config,
        stratification_config = stratification_config,
        track_names = track_names,
        track_intervals = track_intervals,
        reference_dict = reference_dict,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_sv_concordance
    }

    call tasks_cohort.SortVcf as Sort {
      input:
        vcf = Concordance.out_unsorted,
        outfile_prefix = "~{output_prefix}.sorted",
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_sort_vcf
    }
  }

  output {
    Array[File] conc_vcfs = Sort.out
    Array[File] conc_vcf_idxs = Sort.out_index
  }
}
