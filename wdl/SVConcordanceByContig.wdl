version 1.0

import "Structs.wdl"
import "SVConcordance.wdl" as svconcordance

workflow SVConcordanceByContig {
  input {
    # Vcfs must be formatted using FormatVcfForGatk (if unsure, check for ECN FORMAT field)
    # Per-contig vcfs, ordered to match contig_list. eval_vcfs must be indexed (.tbi)
    Array[File] eval_vcfs
    Array[File] truth_vcfs
    String output_prefix

    File contig_list
    File reference_dict

    Int n_per_split = 10000

    # Stratification parameters
    File? clustering_config
    File? stratification_config
    Array[String]? track_names
    Array[File]? track_intervals

    String gatk_docker
    String sv_base_mini_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_scatter_vcf
    RuntimeAttr? runtime_attr_sv_concordance
    RuntimeAttr? runtime_attr_concat
  }

  Array[String] contigs = read_lines(contig_list)

  scatter (i in range(length(contigs))) {
    call svconcordance.SVConcordance {
      input:
        eval_vcf=eval_vcfs[i],
        truth_vcf=truth_vcfs[i],
        output_prefix="~{output_prefix}.~{contigs[i]}",
        contig=contigs[i],
        reference_dict=reference_dict,
        n_per_split=n_per_split,
        clustering_config=clustering_config,
        stratification_config=stratification_config,
        track_names=track_names,
        track_intervals=track_intervals,
        java_mem_fraction=java_mem_fraction,
        gatk_docker=gatk_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_scatter_vcf=runtime_attr_scatter_vcf,
        runtime_attr_sv_concordance=runtime_attr_sv_concordance,
        runtime_attr_concat=runtime_attr_concat
    }
  }

  output {
    Array[File] concordance_vcfs = SVConcordance.concordance_vcf
    Array[File] concordance_vcf_indexes = SVConcordance.concordance_vcf_index
  }
}
