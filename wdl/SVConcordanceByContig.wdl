version 1.0

import "Structs.wdl"
import "SVConcordance.wdl" as svconcordance

workflow SVConcordanceByContig {
  input {
    # Vcfs must be formatted using FormatVcfForGatk (if unsure, check for ECN FORMAT field)
    File eval_vcf
    File truth_vcf
    String output_prefix

    String contig
    File reference_dict

    # Stratification parameters
    File? clustering_config
    File? stratification_config
    Array[String]? track_names
    Array[File]? track_intervals

    String gatk_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_sv_concordance
  }

  call svconcordance.SVConcordanceTask as SVConcordance {
    input:
      eval_vcf=eval_vcf,
      truth_vcf=truth_vcf,
      output_prefix="~{output_prefix}.concordance.~{contig}",
      contig=contig,
      clustering_config=clustering_config,
      stratification_config=stratification_config,
      track_names=track_names,
      track_intervals=track_intervals,
      reference_dict=reference_dict,
      java_mem_fraction=java_mem_fraction,
      gatk_docker=gatk_docker,
      runtime_attr_override=runtime_attr_sv_concordance
  }

  output {
    File concordance_vcf = SVConcordance.out
    File concordance_vcf_index = SVConcordance.out_index
  }
}
