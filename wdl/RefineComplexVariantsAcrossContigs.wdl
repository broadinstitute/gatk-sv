version 1.0

import "Structs.wdl"
import "RefineComplexVariants2.wdl" as per_contig

workflow RefineComplexVariantsAcrossContigs {
  input {
    Array[File] vcfs
    File primary_contigs_list
    String prefix

    Array[File] batch_sample_lists
    Array[File] pe_files
    Array[File] depth_del_beds
    Array[File] depth_dup_beds

    Int n_per_split = 10000
    Int min_pe_cpx = 3
    Int min_pe_ctx = 3

    String gatk_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_scatter_vcf
    RuntimeAttr? runtime_attr_refine
    RuntimeAttr? runtime_attr_concat
  }

  Array[String] contigs = read_lines(primary_contigs_list)

  scatter (i in range(length(vcfs))) {
    call per_contig.RefineComplexVariants {
      input:
        vcf = vcfs[i],
        prefix = "~{prefix}.~{contigs[i]}",
        batch_sample_lists = batch_sample_lists,
        pe_files = pe_files,
        depth_del_beds = depth_del_beds,
        depth_dup_beds = depth_dup_beds,
        n_per_split = n_per_split,
        min_pe_cpx = min_pe_cpx,
        min_pe_ctx = min_pe_ctx,
        gatk_docker = gatk_docker,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_scatter_vcf = runtime_attr_scatter_vcf,
        runtime_attr_refine = runtime_attr_refine,
        runtime_attr_concat = runtime_attr_concat
    }
  }

  output {
    Array[File] cpx_refined_vcfs = RefineComplexVariants.cpx_refined_vcf
    Array[File] cpx_refined_vcf_indexes = RefineComplexVariants.cpx_refined_vcf_index
  }
}
