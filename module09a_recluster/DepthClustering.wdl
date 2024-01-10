version 1.0

import "Structs.wdl"
import "TasksClusterBatch.wdl" as tasks_cluster
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow ClusterDepth {
  input {
    File del_bed
    File dup_bed
    String batch
    File ploidy_table

    Int records_per_bed_shard

    File contig_list
    File sample_list
    File exclude_intervals
    Float exclude_overlap_fraction

    String? clustering_algorithm
    Float depth_interval_overlap

    File? contig_subset_list
    File? gatk_to_svtk_script
    File? cnv_bed_to_gatk_vcf_script

    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_scatter_bed
    RuntimeAttr? runtime_attr_cnv_bed_to_gatk_vcf
    RuntimeAttr? runtime_attr_exclude_intervals_depth
    RuntimeAttr? runtime_attr_svcluster
    RuntimeAttr? runtime_attr_gatk_to_svtk_vcf
    RuntimeAttr? runtime_override_concat_vcfs
  }

  call tasks_cluster.ScatterCompressedBedOmitHeaders as ScatterDel {
    input:
      bed=del_bed,
      records_per_shard=records_per_bed_shard,
      prefix="~{batch}.cluster_batch.depth.del.shard_",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_scatter_bed
  }

  call tasks_cluster.ScatterCompressedBedOmitHeaders as ScatterDup {
    input:
      bed=dup_bed,
      records_per_shard=records_per_bed_shard,
      prefix="~{batch}.cluster_batch.depth.dup.shard_",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_scatter_bed
  }

  scatter (i in range(length(ScatterDel.out))) {
    call tasks_cluster.CNVBedToGatkVcf as DelBedToVcf {
      input:
        bed=ScatterDel.out[i],
        script=cnv_bed_to_gatk_vcf_script,
        contig_list=contig_list,
        sample_list=sample_list,
        ploidy_table=ploidy_table,
        reference_fasta_fai=reference_fasta_fai,
        output_prefix="~{batch}.cluster_batch.depth.gatk_formatted.del.shard_~{i}",
        vid_prefix="~{batch}_raw_depth_DEL_shard_~{i}_",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_cnv_bed_to_gatk_vcf
    }
  }

  scatter (i in range(length(ScatterDup.out))) {
    call tasks_cluster.CNVBedToGatkVcf as DupBedToVcf {
      input:
        bed=ScatterDup.out[i],
        script=cnv_bed_to_gatk_vcf_script,
        contig_list=contig_list,
        sample_list=sample_list,
        ploidy_table=ploidy_table,
        reference_fasta_fai=reference_fasta_fai,
        output_prefix="~{batch}.cluster_batch.depth.gatk_formatted.dup.shard_~{i}",
        vid_prefix="~{batch}_raw_depth_DUP_shard_~{i}_",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_cnv_bed_to_gatk_vcf
    }
  }

  Array[String] contigs = transpose(read_tsv(select_first([contig_subset_list, contig_list])))[0]
  scatter (contig in contigs) {
    call tasks_cluster.SVCluster {
      input:
        vcfs=flatten([DelBedToVcf.out, DupBedToVcf.out]),
        ploidy_table=ploidy_table,
        output_prefix="~{batch}.cluster_batch.depth.~{contig}.clustered",
        contig=contig,
        fast_mode=true,
        algorithm=clustering_algorithm,
        depth_sample_overlap=0,
        depth_interval_overlap=depth_interval_overlap,
        depth_breakend_window=0,
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        variant_prefix="~{batch}_depth_~{contig}_",
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_svcluster
    }
    call tasks_cluster.ExcludeIntervalsByIntervalOverlap {
      input:
        vcf=SVCluster.out,
        overlap_fraction=exclude_overlap_fraction,
        reference_fasta_fai=reference_fasta_fai,
        output_prefix="~{batch}.cluster_batch.depth.~{contig}.exclude_intervals",
        intervals=exclude_intervals,
        intervals_index=exclude_intervals + ".tbi",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_exclude_intervals_depth
    }
    call tasks_cluster.GatkToSvtkVcf {
      input:
        vcf=ExcludeIntervalsByIntervalOverlap.out,
        output_prefix="~{batch}.cluster_batch.depth.~{contig}.svtk_formatted",
        script=gatk_to_svtk_script,
        source="depth",
        contig_list=contig_list,
        remove_formats="CN",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_gatk_to_svtk_vcf
    }
  }

  call tasks_cohort.ConcatVcfs {
    input:
      vcfs=GatkToSvtkVcf.out,
      vcfs_idx=GatkToSvtkVcf.out_index,
      naive=true,
      outfile_prefix="~{batch}.cluster_batch.depth",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_vcfs
  }

  output {
    File clustered_vcf = ConcatVcfs.concat_vcf
    File clustered_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}

