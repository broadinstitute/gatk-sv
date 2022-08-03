version 1.0

import "Structs.wdl"
import "TasksClusterBatch.wdl" as tasks_cluster
import "TasksMakeCohortVcf.wdl" as tasks_cohort
import "SVConcordance.wdl" as svc

# Clusters raw call VCFs across batches - to be used for preparing raw calls for SV concordance analysis

workflow JoinRawCalls {
  input {

    String cohort

    # ClusterBatch outputs
    Array[File]? clustered_manta_vcfs
    Array[File]? clustered_melt_vcfs
    Array[File]? clustered_scramble_vcfs
    Array[File]? clustered_wham_vcfs
    Array[File]? clustered_depth_vcfs

    File ploidy_table

    String? preprocess_args

    File contig_list
    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_prepare_truth
    RuntimeAttr? runtime_attr_svcluster
    RuntimeAttr? runtime_override_concat_vcfs_pesr
  }

  Array[File] vcfs_ = flatten(select_all([clustered_manta_vcfs, clustered_melt_vcfs, clustered_scramble_vcfs, clustered_wham_vcfs, clustered_depth_vcfs]))
  scatter (i in range(length(vcfs_))) {
    call svc.PreprocessVcf {
      input:
        vcf=vcfs_[i],
        ploidy_table=ploidy_table,
        args=preprocess_args,
        output_prefix="~{cohort}.join_raw_calls.preprocess_~{i}",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_prepare_truth
    }
  }

  Array[String] contigs = transpose(read_tsv(contig_list))[0]
  scatter (contig in contigs) {
    call tasks_cluster.SVCluster {
      input:
        vcfs=PreprocessVcf.out,
        ploidy_table=ploidy_table,
        output_prefix="~{cohort}.join_raw_calls.~{contig}",
        contig=contig,
        fast_mode=true,
        algorithm="SINGLE_LINKAGE",
        pesr_sample_overlap=0,
        mixed_sample_overlap=0,
        depth_sample_overlap=0,
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        variant_prefix="~{cohort}_~{contig}_",
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_svcluster
    }
  }

  call tasks_cohort.ConcatVcfs {
    input:
      vcfs=SVCluster.out,
      vcfs_idx=SVCluster.out_index,
      naive=true,
      outfile_prefix="~{cohort}.join_raw_calls",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_vcfs_pesr
  }

  output {
    File joined_raw_calls_vcf = ConcatVcfs.concat_vcf
    File joined_raw_calls_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}
