version 1.0

import "Structs.wdl"
import "FormatVcfForGatk.wdl" as format
import "TasksClusterBatch.wdl" as tasks_cluster
import "TasksMakeCohortVcf.wdl" as tasks_cohort

# Clusters raw call VCFs across batches - to be used for preparing raw calls for SV concordance analysis

workflow JoinRawCalls {
  input {

    String prefix

    # ClusterBatch outputs
    Array[File]? clustered_manta_vcfs
    Array[File]? clustered_manta_vcf_indexes
    Array[File]? clustered_melt_vcfs
    Array[File]? clustered_melt_vcf_indexes
    Array[File]? clustered_scramble_vcfs
    Array[File]? clustered_scramble_vcf_indexes
    Array[File]? clustered_wham_vcfs
    Array[File]? clustered_wham_vcf_indexes
    Array[File]? clustered_depth_vcfs
    Array[File]? clustered_depth_vcf_indexes

    File ped_file

    File contig_list
    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    String? chr_x
    String? chr_y

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_create_ploidy
    RuntimeAttr? runtime_override_concat_input_vcfs
    RuntimeAttr? runtime_attr_prepare_truth
    RuntimeAttr? runtime_attr_svcluster
    RuntimeAttr? runtime_override_concat_vcfs_pesr
  }

  call tasks_cluster.CreatePloidyTableFromPed {
    input:
      ped_file=ped_file,
      contig_list=contig_list,
      retain_female_chr_y=false,
      chr_x=chr_x,
      chr_y=chr_y,
      output_prefix="~{prefix}.ploidy",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_create_ploidy
  }

  Array[Array[File]] vcf_matrix = transpose(select_all([clustered_manta_vcfs, clustered_melt_vcfs, clustered_scramble_vcfs, clustered_wham_vcfs, clustered_depth_vcfs]))
  Array[Array[File]] vcf_index_matrix = transpose(select_all([clustered_manta_vcf_indexes, clustered_melt_vcf_indexes, clustered_scramble_vcf_indexes, clustered_wham_vcf_indexes, clustered_depth_vcf_indexes]))
  scatter (i in range(length(vcf_matrix))) {
    call tasks_cohort.ConcatVcfs as ConcatInputVcfs {
      input:
        vcfs=vcf_matrix[i],
        vcfs_idx=vcf_index_matrix[i],
        allow_overlaps=true,
        outfile_prefix="~{prefix}.join_raw_calls.concat_batch_~{i}",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_concat_input_vcfs
    }
  }

  scatter (i in range(length(ConcatInputVcfs.concat_vcf))) {
    call format.FormatVcfForGatk {
      input:
        vcf=ConcatInputVcfs.concat_vcf[i],
        ped_file=ped_file,
        contig_list=contig_list,
        prefix="~{prefix}.join_raw_calls.format_~{i}",
        sv_pipeline_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker
    }
  }

  scatter (contig in read_lines(contig_list)) {
    call tasks_cluster.SVCluster {
      input:
        vcfs=FormatVcfForGatk.gatk_formatted_vcf,
        ploidy_table=CreatePloidyTableFromPed.out,
        output_prefix="~{prefix}.join_raw_calls.~{contig}",
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
        variant_prefix="~{prefix}_~{contig}_",
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_svcluster
    }
  }

  call tasks_cohort.ConcatVcfs {
    input:
      vcfs=SVCluster.out,
      vcfs_idx=SVCluster.out_index,
      naive=true,
      outfile_prefix="~{prefix}.join_raw_calls",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_vcfs_pesr
  }

  output {
    File joined_raw_calls_vcf = ConcatVcfs.concat_vcf
    File joined_raw_calls_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}
