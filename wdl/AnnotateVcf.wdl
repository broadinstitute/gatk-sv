version 1.0

import "Structs.wdl"
import "ShardedAnnotateVcf.wdl" as sharded_annotate_vcf
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "HailMerge.wdl" as HailMerge

workflow AnnotateVcf {

  input {
    Array[File] vcf_list  # Must be either single full VCF (array of length 1) or array of VCFs sharded by contig. Outputs will match
    Array[File] vcf_idx_list
    File contig_list
    Array[String] prefix_list
    Boolean sharded_by_contig  # True if providing a vcf_list sharded by contig. False if providing a single full VCF

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
    RuntimeAttr? runtime_attr_concat_vcfs
    RuntimeAttr? runtime_attr_shard_vcf
    RuntimeAttr? runtime_attr_subset_vcf_by_samples_list
    RuntimeAttr? runtime_attr_compute_AFs
    RuntimeAttr? runtime_attr_modify_vcf
    RuntimeAttr? runtime_attr_split_ref_bed
    RuntimeAttr? runtime_attr_split_query_vcf
    RuntimeAttr? runtime_attr_bedtools_closest
    RuntimeAttr? runtime_attr_select_matched_svs
    RuntimeAttr? runtime_attr_concat_sharded_cluster
    RuntimeAttr? runtime_attr_preconcat_sharded_cluster
    RuntimeAttr? runtime_attr_hail_merge_sharded_cluster
    RuntimeAttr? runtime_attr_fix_header_sharded_cluster
  }

  Array[String] contigs = read_lines(contig_list)

  scatter (i in range(length(contigs))) {
    Int array_index = if (sharded_by_contig) then i else 0
    call sharded_annotate_vcf.ShardedAnnotateVcf {
      input:
        vcf = vcf_list[array_index],
        vcf_idx = vcf_idx_list[array_index],
        contig = contigs[i],
        prefix = prefix_list[array_index],
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
        runtime_attr_subset_vcf_by_samples_list = runtime_attr_subset_vcf_by_samples_list,
        runtime_attr_compute_AFs  = runtime_attr_compute_AFs,
        runtime_attr_modify_vcf = runtime_attr_modify_vcf,
        runtime_attr_split_ref_bed  = runtime_attr_split_ref_bed,
        runtime_attr_split_query_vcf  = runtime_attr_split_query_vcf,
        runtime_attr_bedtools_closest = runtime_attr_bedtools_closest,
        runtime_attr_select_matched_svs = runtime_attr_select_matched_svs,
        runtime_attr_concat_sharded_cluster = runtime_attr_concat_sharded_cluster,
        runtime_attr_preconcat_sharded_cluster  = runtime_attr_preconcat_sharded_cluster,
        runtime_attr_hail_merge_sharded_cluster = runtime_attr_hail_merge_sharded_cluster,
        runtime_attr_fix_header_sharded_cluster = runtime_attr_fix_header_sharded_cluster
    }
  }

  # Concat VCFs to the contig level or fully depending on format of input
  # ShardedAnnotateVcf.sharded_annotated_vcf is is an Array[Array[File]] with one inner Array[File] of shards per contig
  Array[Array[File]] vcfs_for_concatenation = if sharded_by_contig then ShardedAnnotateVcf.sharded_annotated_vcf else [flatten(ShardedAnnotateVcf.sharded_annotated_vcf)]
  Array[Array[File]] vcf_idxs_for_concatenation = if sharded_by_contig then ShardedAnnotateVcf.sharded_annotated_vcf_idx else [flatten(ShardedAnnotateVcf.sharded_annotated_vcf_idx)]
  if (use_hail) {
    scatter (i in range(length(vcfs_for_concatenation))) {
      call HailMerge.HailMerge {
        input:
          vcfs=vcfs_for_concatenation[i],
          prefix="~{prefix_list[i]}.annotated",
          gcs_project=gcs_project,
          sv_base_mini_docker=sv_base_mini_docker,
          sv_pipeline_docker=sv_pipeline_docker,
          sv_pipeline_hail_docker=select_first([sv_pipeline_hail_docker]),
          runtime_override_preconcat=runtime_attr_preconcat_sharded_cluster,
          runtime_override_hail_merge=runtime_attr_hail_merge_sharded_cluster,
          runtime_override_fix_header=runtime_attr_fix_header_sharded_cluster
      }
    }
  }

  if (!use_hail) {
    scatter (i in range(length(vcfs_for_concatenation))) {
      call MiniTasks.ConcatVcfs {
        input:
          vcfs=vcfs_for_concatenation[i],
          vcfs_idx=vcf_idxs_for_concatenation[i],
          allow_overlaps=true,
          outfile_prefix="~{prefix_list[i]}.annotated",
          sv_base_mini_docker=sv_base_mini_docker,
          runtime_attr_override=runtime_attr_concat_sharded_cluster
      }
    }
  }

  output {
    Array[File] output_vcf_list = select_first([ConcatVcfs.concat_vcf, HailMerge.merged_vcf])
    Array[File] output_vcf_idx_list = select_first([ConcatVcfs.concat_vcf_idx, HailMerge.merged_vcf_index])
  }
}