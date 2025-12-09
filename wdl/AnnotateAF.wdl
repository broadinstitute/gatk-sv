version 1.0

import "Structs.wdl"
import "ShardedAnnotateVcf.wdl" as sharded_annotate_vcf
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow AnnotateAF {

  input {
    File vcf
    File vcf_index
    Array[String] contigs
    String prefix

    File? protein_coding_gtf
    File? noncoding_bed
    Int? promoter_window
    Int? max_breakend_as_cnv_length
    String? svannotate_additional_args

    File? sample_pop_assignments
    File? sample_keep_list
    File? ped_file
    File? par_bed
    File? allosomes_list
    Int   sv_per_shard

    Boolean annotate_external_af = true
    Boolean annotate_internal_af = true
    Boolean annotate_functional_consequences = true

    File? external_af_ref_bed
    String? external_af_ref_prefix
    Array[String]? external_af_population

    String annotate_af_docker
    String sv_base_mini_docker
    String gatk_docker

    RuntimeAttr? runtime_attr_svannotate
    RuntimeAttr? runtime_attr_scatter_vcf
    RuntimeAttr? runtime_attr_subset_vcf_by_samples_list
    RuntimeAttr? runtime_attr_compute_AFs
    RuntimeAttr? runtime_attr_modify_vcf
    RuntimeAttr? runtime_attr_split_ref_bed
    RuntimeAttr? runtime_attr_split_query_vcf
    RuntimeAttr? runtime_attr_bedtools_closest
    RuntimeAttr? runtime_attr_select_matched_svs
    RuntimeAttr? runtime_attr_concat
    RuntimeAttr? runtime_attr_preconcat
    RuntimeAttr? runtime_attr_fix_header
  }

  scatter (contig in contigs) {
    call sharded_annotate_vcf.ShardedAnnotateVcf {
      input:
        vcf = vcf,
        vcf_idx = vcf_index,
        contig = contig,
        prefix = prefix,
        protein_coding_gtf = protein_coding_gtf,
        noncoding_bed = noncoding_bed,
        promoter_window = promoter_window,
        svannotate_additional_args = svannotate_additional_args,
        max_breakend_as_cnv_length = max_breakend_as_cnv_length,

        annotate_external_af = annotate_external_af,
        annotate_internal_af = annotate_internal_af,
        annotate_functional_consequences = annotate_functional_consequences,

        sample_pop_assignments = sample_pop_assignments,
        sample_keep_list = sample_keep_list,
        ped_file = ped_file,
        par_bed = par_bed,
        sv_per_shard = sv_per_shard,
        allosomes_list = allosomes_list,

        ref_bed = external_af_ref_bed,
        ref_prefix = external_af_ref_prefix,
        population = external_af_population,

        gatk_docker = gatk_docker,
        sv_pipeline_docker = annotate_af_docker,
        sv_base_mini_docker = sv_base_mini_docker,

        runtime_attr_svannotate = runtime_attr_svannotate,
        runtime_attr_scatter_vcf = runtime_attr_scatter_vcf,
        runtime_attr_subset_vcf_by_samples_list = runtime_attr_subset_vcf_by_samples_list,
        runtime_attr_compute_AFs  = runtime_attr_compute_AFs,
        runtime_attr_modify_vcf = runtime_attr_modify_vcf,
        runtime_attr_split_ref_bed  = runtime_attr_split_ref_bed,
        runtime_attr_split_query_vcf  = runtime_attr_split_query_vcf,
        runtime_attr_bedtools_closest = runtime_attr_bedtools_closest,
        runtime_attr_select_matched_svs = runtime_attr_select_matched_svs
    }
  }

  # ShardedAnnotateVcf.sharded_annotated_vcf is is an Array[Array[File]] with one inner Array[File] of shards per contig
  Array[File] vcfs_for_concatenation = flatten(ShardedAnnotateVcf.sharded_annotated_vcf)
  Array[File] vcf_idxs_for_concatenation = flatten(ShardedAnnotateVcf.sharded_annotated_vcf_idx)

  call MiniTasks.ConcatVcfs {
    input:
      vcfs=vcfs_for_concatenation,
      vcfs_idx=vcf_idxs_for_concatenation,
      allow_overlaps=true,
      outfile_prefix="~{prefix}.annotated",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_concat
  }

  output {
    File af_annotated_vcf = ConcatVcfs.concat_vcf
    File af_annotated_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}
