version 1.0

import "Structs.wdl"
import "ShardedAnnotateVcf.wdl" as sharded_annotate_vcf
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "HailMerge.wdl" as HailMerge

workflow AnnotateVcf {

  input {
    File vcf  # GATK-SV VCF for annotation. Index .tbi must be located at the same path
    File contig_list  # Ordered list of contigs to annotate that are present in the input VCF
    String prefix

    File protein_coding_gtf
    File? noncoding_bed
    Int? promoter_window
    Int? max_breakend_as_cnv_length
    String? svannotate_additional_args

    File? sample_pop_assignments  # Two-column file with sample ID & pop assignment. "." for pop will ignore sample
    File? sample_keep_list              # List of samples to be retained from the output vcf
    File? ped_file                # Used for M/F AF calculations
    File? par_bed
    File? allosomes_list
    Int   sv_per_shard

    File? ref_bed              # File with external allele frequencies
    String? ref_prefix         # prefix name for external AF call set (required if ref_bed set)
    Array[String]? population  # populations to annotate external AF for (required if ref_bed set)

    Boolean use_hail
    String? gcs_project

    String sv_pipeline_docker
    String? sv_pipeline_hail_docker
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
    RuntimeAttr? runtime_attr_hail_merge
    RuntimeAttr? runtime_attr_fix_header
  }

  Array[String] contigs = read_lines(contig_list)

  scatter (contig in contigs) {
    call sharded_annotate_vcf.ShardedAnnotateVcf {
      input:
        vcf = vcf,
        vcf_idx = vcf + ".tbi",
        contig = contig,
        prefix = prefix,
        protein_coding_gtf = protein_coding_gtf,
        noncoding_bed = noncoding_bed,
        promoter_window = promoter_window,
        svannotate_additional_args = svannotate_additional_args,
        max_breakend_as_cnv_length = max_breakend_as_cnv_length,

        sample_pop_assignments = sample_pop_assignments,
        sample_keep_list = sample_keep_list,
        ped_file = ped_file,
        par_bed = par_bed,
        sv_per_shard = sv_per_shard,
        allosomes_list = allosomes_list,

        ref_bed = ref_bed,
        ref_prefix = ref_prefix,
        population = population,

        use_hail = use_hail,
        gcs_project = gcs_project,

        gatk_docker = gatk_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_hail_docker = sv_pipeline_hail_docker,

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

  # Concat VCF shards with or without hail
  # ShardedAnnotateVcf.sharded_annotated_vcf is is an Array[Array[File]] with one inner Array[File] of shards per contig
  Array[File] vcfs_for_concatenation = flatten(ShardedAnnotateVcf.sharded_annotated_vcf)
  Array[File] vcf_idxs_for_concatenation = flatten(ShardedAnnotateVcf.sharded_annotated_vcf_idx)
  if (use_hail) {
    call HailMerge.HailMerge {
      input:
        vcfs=vcfs_for_concatenation,
        prefix="~{prefix}.annotated",
        gcs_project=gcs_project,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_hail_docker=select_first([sv_pipeline_hail_docker]),
        runtime_override_preconcat=runtime_attr_preconcat,
        runtime_override_hail_merge=runtime_attr_hail_merge,
        runtime_override_fix_header=runtime_attr_fix_header
    }
  }

  if (!use_hail) {
    call MiniTasks.ConcatVcfs {
      input:
        vcfs=vcfs_for_concatenation,
        vcfs_idx=vcf_idxs_for_concatenation,
        allow_overlaps=true,
        outfile_prefix="~{prefix}.annotated",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_concat
    }
  }

  output {
    File annotated_vcf = select_first([ConcatVcfs.concat_vcf, HailMerge.merged_vcf])
    File annotated_vcf_index = select_first([ConcatVcfs.concat_vcf_idx, HailMerge.merged_vcf_index])
  }
}
