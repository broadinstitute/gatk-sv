version 1.0

import "SVConcordance.wdl" as conc
import "Utils.wdl" as utils

workflow SVConcordanceSimple {
  input {
    File eval_vcf
    File eval_vcf_idx
    File truth_vcf
    File truth_vcf_idx

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
    Int? stratify_num_breakpoint_overlaps = 0
    Float? stratify_overlap_fraction = 0.5

    File reference_dict

    String output_prefix
    File? sample_ids

    String sv_base_mini_docker
    String gatk_docker

    RuntimeAttr? runtime_attr_sv_concordance
    RuntimeAttr? runtime_attr_sort_vcf
    RuntimeAttr? runtime_attr_subset_vcf
  }

  # If sample_ids is provided, subset both VCFs to only include those samples
  if (defined(sample_ids)) {
    call utils.SubsetVcfBySamplesList as SubsetEvalVcf {
      input:
        vcf = eval_vcf,
        vcf_idx = eval_vcf_idx,
        list_of_samples = select_first([sample_ids]),
        outfile_name = "~{output_prefix}.eval.subset.vcf.gz",
        remove_samples = false,
        remove_private_sites = true,
        keep_af = true,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_subset_vcf
    }

    call utils.SubsetVcfBySamplesList as SubsetTruthVcf {
      input:
        vcf = truth_vcf,
        vcf_idx = truth_vcf_idx,
        list_of_samples = select_first([sample_ids]),
        outfile_name = "~{output_prefix}.truth.subset.vcf.gz",
        remove_samples = false,
        remove_private_sites = true,
        keep_af = true,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_subset_vcf
    }
  }

  call conc.SVConcordanceTask {
      input:
        truth_vcf=select_first([SubsetTruthVcf.vcf_subset, truth_vcf]),
        eval_vcf=select_first([SubsetEvalVcf.vcf_subset, eval_vcf]),
        output_prefix="~{output_prefix}.unsorted",
        additional_args="--pesr-interval-overlap ~{pesr_interval_overlap} --pesr-size-similarity ~{pesr_size_similarity} --pesr-breakend-window ~{pesr_breakend_window} --depth-interval-overlap ~{depth_interval_overlap} --depth-size-similarity ~{depth_size_similarity} --depth-breakend-window ~{depth_breakend_window} --mixed-interval-overlap ~{mixed_interval_overlap} --mixed-size-similarity ~{mixed_size_similarity} --mixed-breakend-window ~{mixed_breakend_window} --stratify-num-breakpoint-overlaps ~{stratify_num_breakpoint_overlaps} --stratify-overlap-fraction ~{stratify_overlap_fraction}",
        clustering_config=clustering_config,
        stratification_config=stratification_config,
        track_names=track_names,
        track_intervals=track_intervals,
        reference_dict=reference_dict,
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_sv_concordance
  }

  output {
    File conc_vcf = SoSVConcordanceTaskrtVcf.concordance_vcf
    File conc_vcf_idx = SVConcordanceTask.concordance_vcf_index
  }
}