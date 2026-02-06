version 1.0

import "Structs.wdl"
import "FormatVcfForGatk.wdl" as format_vcf
import "TasksClusterBatch.wdl" as tasks_cluster
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow MergeBatchSites {
  input {
    String cohort             # Cohort name or project prefix for all cohort-level outputs
    Array[File] ploidy_tables
    Array[File] depth_vcfs    # Filtered depth VCFs across batches
    Array[File] pesr_vcfs     # Filtered PESR VCFs across batches

    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    String sv_base_mini_docker
    String sv_pipeline_docker
    String gatk_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_svtk_to_gatk_vcf
    RuntimeAttr? runtime_attr_merge_pesr
    RuntimeAttr? runtime_attr_merge_depth
    RuntimeAttr? runtime_override_concat
  }

  scatter (i in range(length(depth_vcfs))) {
    call format_vcf.FormatVcf as FormatDepth {
      input:
        vcf=depth_vcfs[i],
        ploidy_table=ploidy_tables[i],
        output_prefix=basename(depth_vcfs[i]) + ".reformatted",
        args="--add-sr-pos",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_svtk_to_gatk_vcf
    }
  }

  scatter (i in range(length(pesr_vcfs))) {
    call format_vcf.FormatVcf as FormatPesr {
      input:
        vcf=pesr_vcfs[i],
        ploidy_table=ploidy_tables[i],
        output_prefix=basename(pesr_vcfs[i]) + ".reformatted",
        args="--add-sr-pos",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_svtk_to_gatk_vcf
    }
  }

  call tasks_cluster.SVCluster as MergePesr {
    input:
      vcfs=FormatPesr.out,
      output_prefix="~{cohort}.merge_batch_sites.pesr",
      sites_only=true,
      omit_members=true,
      algorithm="SINGLE_LINKAGE",
      depth_sample_overlap=0,
      depth_interval_overlap=1,
      depth_size_similarity=1,
      depth_breakend_window=0,
      mixed_sample_overlap=0,
      mixed_interval_overlap=1,
      mixed_size_similarity=1,
      mixed_breakend_window=0,
      pesr_sample_overlap=0,
      pesr_interval_overlap=1,
      pesr_size_similarity=1,
      pesr_breakend_window=0,
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_fasta_fai,
      reference_dict=reference_dict,
      java_mem_fraction=java_mem_fraction,
      gatk_docker=gatk_docker,
      runtime_attr_override=runtime_attr_merge_pesr
  }

  call tasks_cluster.SVCluster as MergeDepth {
    input:
      vcfs=FormatDepth.out,
      output_prefix="~{cohort}.merge_batch_sites.depth",
      sites_only=true,
      omit_members=true,
      algorithm="SINGLE_LINKAGE",
      depth_sample_overlap=0,
      depth_interval_overlap=1,
      depth_size_similarity=1,
      depth_breakend_window=0,
      mixed_sample_overlap=0,
      mixed_interval_overlap=1,
      mixed_size_similarity=1,
      depth_breakend_window=0,
      pesr_sample_overlap=0,
      pesr_interval_overlap=1,
      pesr_size_similarity=1,
      pesr_breakend_window=0,
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_fasta_fai,
      reference_dict=reference_dict,
      java_mem_fraction=java_mem_fraction,
      gatk_docker=gatk_docker,
      runtime_attr_override=runtime_attr_merge_depth
  }

  call tasks_cohort.ConcatVcfs {
    input:
      vcfs = [MergePesr.out, MergeDepth.out],
      vcfs_idx = [MergePesr.out_index, MergeDepth.out_index],
      allow_overlaps = true,
      outfile_prefix = cohort + ".merge_batch_sites",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_concat
  }

  output {
    File merge_batch_sites_vcf = ConcatVcfs.concat_vcf
    File merge_batch_sites_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}
