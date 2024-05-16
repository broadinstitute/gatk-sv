version 1.0

import "Structs.wdl"
import "SVConcordance.wdl" as svc
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow SVConcordanceScattered {
  input {
    # Either a single whole-genome vcf or a list of contig-sharded vcfs
    Array[File] eval_vcfs
    # Either a single whole-genome vcf or a list of contig-sharded vcfs
    Array[File] truth_vcfs
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

  Array[String] contigs = read_lines(contig_list)
  scatter (i in range(length(contigs))) {
    call svc.SVConcordanceTask {
      input:
        eval_vcf=if length(eval_vcfs) == 1 then eval_vcfs[0] else eval_vcfs[i],
        truth_vcf=if length(truth_vcfs) == 1 then truth_vcfs[0] else truth_vcfs[i],
        output_prefix="~{output_prefix}.concordance.~{contigs[i]}.unsorted",
        contig=contigs[i],
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_sv_concordance
    }

    call tasks_cohort.SortVcf {
      input:
        vcf=SVConcordanceTask.out_unsorted,
        outfile_prefix="~{output_prefix}.concordance.~{contigs[i]}",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_sort_vcf
    }
  }

  output {
    Array[File] concordance_vcf = SortVcf.out
    Array[File] concordance_vcf_index = SortVcf.out_index
  }
}
