version 1.0

import "Structs.wdl"
import "Utils.wdl" as utils
import "SVConcordance.wdl" as concordance

workflow SVConcordanceCNVBySample {
  input {
    # Input VCFs - must be formatted using FormatVcfForGatk (if unsure, check for ECN FORMAT field)
    File vcf_eval
    File vcf_eval_index
    File vcf_truth
    File vcf_truth_index
    
    # Sample IDs to process
    Array[String] sample_ids
    String output_prefix
    File reference_dict

    # Optional clustering parameters
    Float? pesr_interval_overlap = 0.5
    Float? pesr_size_similarity = 0.0
    Int? pesr_breakend_window = 5000

    Float? depth_interval_overlap = 0.5
    Float? depth_size_similarity = 0.0
    Int? depth_breakend_window = 10000000

    Float? mixed_interval_overlap = 0.5
    Float? mixed_size_similarity = 0.0
    Int? mixed_breakend_window = 10000000

    # Optional stratification parameters
    File? clustering_config
    File? stratification_config
    Array[String]? track_names
    Array[File]? track_intervals
    Int? stratify_num_breakpoint_overlaps = 0
    Float? stratify_overlap_fraction = 0.5
    
    # Docker images
    String sv_base_mini_docker
    String linux_docker
    String gatk_docker
    
    # Runtime parameters
    RuntimeAttr? runtime_attr_subset_vcf
    RuntimeAttr? runtime_attr_sv_concordance
  }
  
  # Convert SVTYPE=DUP to INS in both vcf_eval and vcf_truth
  call utils.ConvertDupToIns as ConvertEvalVcf {
    input:
      vcf=vcf_eval,
      outfile_name="~{output_prefix}.eval.DUPtoINS.vcf.gz",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_subset_vcf
  }
  call utils.ConvertDupToIns as ConvertTruthVcf {
    input:
      vcf=vcf_truth,
      outfile_name="~{output_prefix}.truth.DUPtoINS.vcf.gz",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_subset_vcf
  }

  scatter (sample_id in sample_ids) {
    call utils.SubsetVcfToSample as SubsetEvalVcf {
      input:
        vcf=ConvertEvalVcf.out,
        sample=sample_id,
        outfile_name="~{output_prefix}.~{sample_id}.eval.vcf.gz",
        remove_sample = false,
        remove_private_sites=true,
        keep_af = true,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_subset_vcf
    }
    
    # Subset truth VCF to current sample
    call utils.SubsetVcfToSample as SubsetTruthVcf {
      input:
        vcf=ConvertTruthVcf.out,
        sample=sample_id,
        outfile_name="~{output_prefix}.~{sample_id}.truth.vcf.gz",
        remove_private_sites=true,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_subset_vcf
    }
    
    call concordance.SVConcordanceTask as SVConcordanceSample {
      input:
        eval_vcf=SubsetEvalVcf.vcf_subset,
        truth_vcf=SubsetTruthVcf.vcf_subset,
        output_prefix="~{output_prefix}.~{sample_id}",
        additional_args="--pesr-interval-overlap ~{pesr_interval_overlap} --pesr-size-similarity ~{pesr_size_similarity} --pesr-breakend-window ~{pesr_breakend_window} --depth-interval-overlap ~{depth_interval_overlap} --depth-size-similarity ~{depth_size_similarity} --depth-breakend-window ~{depth_breakend_window} --mixed-interval-overlap ~{mixed_interval_overlap} --mixed-size-similarity ~{mixed_size_similarity} --mixed-breakend-window ~{mixed_breakend_window} --stratify-num-breakpoint-overlaps ~{stratify_num_breakpoint_overlaps} --stratify-overlap-fraction ~{stratify_overlap_fraction}",
        clustering_config=clustering_config,
        stratification_config=stratification_config,
        track_names=track_names,
        track_intervals=track_intervals,
        reference_dict=reference_dict,
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_sv_concordance,
    }
  }
  
  call utils.TarFiles as TarConcordanceFiles {
    input:
      files=flatten([SVConcordanceSample.out, SVConcordanceSample.out_index]),
      prefix=output_prefix + ".concordance_results",
      linux_docker=linux_docker,
      runtime_attr_override=runtime_attr_subset_vcf
  }
  
  output {
    File concordance_results_tar = TarConcordanceFiles.out
  }
}