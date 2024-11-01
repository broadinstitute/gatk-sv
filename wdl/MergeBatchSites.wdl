version 1.0

import "Structs.wdl"
import "TasksClusterBatch.wdl" as tasks_cluster

workflow MergeBatchSites {
  input {
    String cohort             # Cohort name or project prefix for all cohort-level outputs
    Array[File] depth_vcfs    # Filtered depth VCFs across batches
    Array[File] pesr_vcfs     # Filtered PESR VCFs across batches
    Array[File]? ref_pop_depth_vcfs   # Reference population depth VCFs
    Array[File]? ref_pop_pesr_vcfs   # Reference population PESR VCFs

    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    String gatk_docker
    Float? java_mem_fraction
    RuntimeAttr? runtime_attr_merge_pesr
    RuntimeAttr? runtime_attr_merge_depth
  }

  call tasks_cluster.SVCluster as MergePesr {
    input:
      vcfs=flatten(select_all([pesr_vcfs, ref_pop_pesr_vcfs])),
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
      depth_breakend_window=0,
      pesr_sample_overlap=0,
      pesr_interval_overlap=1,
      pesr_size_similarity=1,
      pesr_breakend_window=0,
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_fasta_fai,
      reference_dict=reference_dict,
      java_mem_fraction=java_mem_fraction,
      variant_prefix="~{cohort}_mbs_pesr",
      gatk_docker=gatk_docker,
      runtime_attr_override=runtime_attr_merge_pesr
  }

  call tasks_cluster.SVCluster as MergeDepth {
    input:
      vcfs=flatten(select_all([depth_vcfs, ref_pop_depth_vcfs])),
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
      variant_prefix="~{cohort}_mbs_depth",
      gatk_docker=gatk_docker,
      runtime_attr_override=runtime_attr_merge_depth
  }

  output {
    File cohort_pesr_vcf = MergePesr.out
    File cohort_pesr_vcf_index = MergePesr.out_index
    File cohort_depth_vcf = MergeDepth.out
    File cohort_depth_vcf_index = MergeDepth.out_index
  }
}
