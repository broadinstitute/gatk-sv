version 1.0

import "Structs.wdl"
import "SVConcordance.wdl" as svc

workflow SVConcordanceScattered {
  input {
    Array[File] eval_vcfs
    Array[File] truth_vcf
    String output_prefix

    File contig_list
    File reference_dict

    String gatk_docker
    String sv_base_mini_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_sv_concordance
    RuntimeAttr? runtime_attr_sort_vcf
    RuntimeAttr? runtime_override_concat_shards
  }

  scatter (i in range(length(eval_vcfs))) {
    call svc.SVConcordance {
      input:
        eval_vcf=eval_vcfs[i],
        truth_vcf=truth_vcf[i],
        output_prefix=output_prefix,
        contig_list=contig_list,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        sv_base_mini_docker=sv_base_mini_docker,
        gatk_docker=gatk_docker,
        runtime_attr_sv_concordance=runtime_attr_sv_concordance,
        runtime_attr_sort_vcf=runtime_attr_sort_vcf,
        runtime_override_concat_shards=runtime_override_concat_shards
    }
  }

  output {
    Array[File] concordance_vcf = SVConcordance.concordance_vcf
    Array[File] concordance_vcf_index = SVConcordance.concordance_vcf_index
  }
}
