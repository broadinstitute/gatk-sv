version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "HailMerge.wdl" as HailMerge
import "AnnotateFunctionalConsequences.wdl" as func
import "PruneAndAddVafs.wdl" as pav
import "AnnotateExternalAF.wdl" as eaf

workflow ShardedAnnotateVcf {

  input {
    File vcf
    File vcf_idx
    String prefix
    String contig

    File protein_coding_gtf
    File? noncoding_bed
    Int? promoter_window
    Int? max_breakend_as_cnv_length
    String? svannotate_additional_args

    Int max_shards_per_chrom_step1
    Int min_records_per_shard_step1

    File? sample_pop_assignments  # Two-column file with sample ID & pop assignment. "." for pop will ignore sample
    File sample_list
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
    RuntimeAttr? runtime_attr_compute_AFs
    RuntimeAttr? runtime_attr_combine_vcfs
    RuntimeAttr? runtime_attr_modify_vcf
    RuntimeAttr? runtime_attr_combine_vcfs
    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_split_ref_bed
    RuntimeAttr? runtime_attr_split_query_vcf
    RuntimeAttr? runtime_attr_bedtools_closest
    RuntimeAttr? runtime_attr_select_matched_svs
    RuntimeAttr? runtime_attr_scatter_vcf
    RuntimeAttr? runtime_attr_fix_ends_rescale_GQ
    RuntimeAttr? runtime_attr_concat_sharded_cluster
    RuntimeAttr? runtime_attr_preconcat_sharded_cluster
    RuntimeAttr? runtime_attr_hail_merge_sharded_cluster
    RuntimeAttr? runtime_attr_fix_header_sharded_cluster
    RuntimeAttr? runtime_attr_get_vcf_header_with_members_info_line
  }

  call MiniTasks.ScatterVcf {
    input:
      vcf = vcf,
      prefix = prefix,
      records_per_shard = sv_per_shard,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_scatter_vcf
  }

  scatter (i in range(length(ScatterVcf.shards))) {

    call func.AnnotateFunctionalConsequences {
      input:
        vcf = ScatterVcf.shards[i],
        vcf_index = ScatterVcf.shards_idx[i],
        prefix = "~{prefix}.~{i}",
        protein_coding_gtf = protein_coding_gtf,
        noncoding_bed = noncoding_bed,
        promoter_window = promoter_window,
        max_breakend_as_cnv_length = max_breakend_as_cnv_length,
        additional_args = svannotate_additional_args,
        gatk_docker = gatk_docker,
        runtime_attr_svannotate = runtime_attr_svannotate
    }

    call pav.PruneAndAddVafs as PruneAndAddVafs {
      input:
        vcf                    = AnnotateFunctionalConsequences.annotated_vcf,
        vcf_idx                = AnnotateFunctionalConsequences.annotated_vcf_index,
        prefix                 = prefix,
        contig                 = contig,
        ped_file               = ped_file,
        par_bed                = par_bed,
        sample_list            = sample_list,
        allosomes_list         = allosomes_list,
        sample_pop_assignments = sample_pop_assignments,

        sv_base_mini_docker     = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_shard_vcf    = runtime_attr_shard_vcf,
        runtime_attr_compute_AFs  = runtime_attr_compute_AFs,
        runtime_attr_combine_vcfs = runtime_attr_combine_vcfs,
        runtime_attr_concat_vcfs  = runtime_attr_concat_vcfs
    }

    if (defined(ref_bed)) {
      call eaf.AnnotateExternalAF as AnnotateExternalAF {
        input:
          vcf     = PruneAndAddVafs.output_vcf,
          vcf_idx = PruneAndAddVafs.output_vcf_idx,
          ref_bed = select_first([ref_bed]),
          population = select_first([population]),
          ref_prefix = select_first([ref_prefix]),
          prefix = prefix,
          contigs = [contig],
          max_shards_per_chrom_step1 = max_shards_per_chrom_step1,
          min_records_per_shard_step1 = min_records_per_shard_step1,
          sv_base_mini_docker = sv_base_mini_docker,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_modify_vcf = runtime_attr_modify_vcf,
          runtime_attr_split_vcf = runtime_attr_split_vcf,
          runtime_attr_combine_vcfs = runtime_attr_combine_vcfs,
          runtime_attr_split_ref_bed = runtime_attr_split_ref_bed,
          runtime_attr_split_query_vcf = runtime_attr_split_query_vcf,
          runtime_attr_bedtools_closest = runtime_attr_bedtools_closest,
          runtime_attr_select_matched_svs = runtime_attr_select_matched_svs
      }
    }

  }

  #Array[File?] sharded_annotated_vcf = select_first([AnnotateExternalAF.annotated_vcf, PruneAndAddVafs.output_vcf])
  #Array[File?] sharded_annotated_vcf_idx = select_first([AnnotateExternalAF.annotated_vcf_tbi, PruneAndAddVafs.output_vcf_idx])
  Array[File] sharded_annotated_vcf = PruneAndAddVafs.output_vcf
  Array[File] sharded_annotated_vcf_idx = PruneAndAddVafs.output_vcf_idx


  if (length(sharded_annotated_vcf) == 0) {
    call MiniTasks.GetVcfHeaderWithMembersInfoLine as GetVcfHeader_annotated {
      input:
        vcf_gz=vcf,
        prefix="~{prefix}.annotated",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_get_vcf_header_with_members_info_line
    }
  }

  if (length(sharded_annotated_vcf) > 0) {
    if (use_hail) {
      call HailMerge.HailMerge as ConcatVcfsHail_annotated {
        input:
          vcfs=sharded_annotated_vcf,
          prefix="~{prefix}.annotated",
          gcs_project=gcs_project,
          sv_base_mini_docker=sv_base_mini_docker,
          sv_pipeline_docker=sv_pipeline_docker,
          sv_pipeline_hail_docker=select_first([sv_pipeline_hail_docker]),
          runtime_attr_preconcat=runtime_attr_preconcat_sharded_cluster,
          runtime_attr_hail_merge=runtime_attr_hail_merge_sharded_cluster,
          runtime_attr_fix_header=runtime_attr_fix_header_sharded_cluster
      }
    }

    if (!use_hail) {
      call MiniTasks.ConcatVcfs as ConcatVcfs_annotated {
        input:
          vcfs=sharded_annotated_vcf,
          vcfs_idx=sharded_annotated_vcf_idx,
          allow_overlaps=true,
          outfile_prefix="~{prefix}.annotatedd",
          sv_base_mini_docker=sv_base_mini_docker,
          runtime_attr_override=runtime_attr_concat_sharded_cluster
      }
    }

  }


  output {
    File output_vcf = select_first([GetVcfHeader_annotated.out, ConcatVcfs_annotated.concat_vcf, ConcatVcfsHail_annotated.merged_vcf])
    File output_vcf_idx = select_first([GetVcfHeader_annotated.out_idx, ConcatVcfs_annotated.concat_vcf_idx, ConcatVcfsHail_annotated.merged_vcf_index])
  }
}

