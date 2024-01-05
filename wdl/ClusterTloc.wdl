version 1.0

import "PESRClustering.wdl" as pesr
import "TasksClusterBatch.wdl" as tasks
import "Utils.wdl" as util

workflow ClusterTloc {
  input {
    String batch

    Array[File] manta_tloc_vcfs

    File ped_file

    # Reference
    File contig_list
    File reference_fasta
    File reference_fasta_fai
    File reference_dict
    String? chr_x
    String? chr_y

    # PESR-based variant clustering
    Int? pesr_min_size
    File pesr_exclude_intervals
    Float pesr_interval_overlap
    Int pesr_breakend_window
    String? pesr_clustering_algorithm

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    String linux_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_tar_files
    RuntimeAttr? runtime_attr_create_ploidy
    RuntimeAttr? runtime_attr_prepare_pesr_vcfs
    RuntimeAttr? runtime_attr_svcluster_manta_tloc
    RuntimeAttr? runtime_override_concat_vcfs_pesr
    RuntimeAttr? runtime_attr_gatk_to_svtk_vcf_pesr
    RuntimeAttr? runtime_attr_exclude_intervals_pesr
  }

  call util.TarFiles {
    input:
      files=manta_tloc_vcfs,
      prefix="{batch}.clustered.manta_tloc",
      linux_docker=linux_docker,
      runtime_attr_override=runtime_attr_tar_files
  }

  # TODO : properly set allosome ploidy, which creates problems in RDTest for allosomes at the moment
  call tasks.CreatePloidyTableFromPed {
    input:
      ped_file=ped_file,
      contig_list=contig_list,
      retain_female_chr_y=true,
      chr_x=chr_x,
      chr_y=chr_y,
      output_prefix="~{batch}.ploidy",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_create_ploidy
  }

  call pesr.ClusterPESR {
    input:
      vcf_tar=TarFiles.out,
      ploidy_table=CreatePloidyTableFromPed.out,
      batch=batch,
      caller="manta_tloc",
      min_size=select_first([pesr_min_size, 50]),
      exclude_intervals=pesr_exclude_intervals,
      contig_list=contig_list,
      pesr_interval_overlap=pesr_interval_overlap,
      pesr_breakend_window=pesr_breakend_window,
      clustering_algorithm=pesr_clustering_algorithm,
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_fasta_fai,
      reference_dict=reference_dict,
      java_mem_fraction=java_mem_fraction,
      gatk_docker=gatk_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_prepare_pesr_vcfs=runtime_attr_prepare_pesr_vcfs,
      runtime_attr_svcluster=runtime_attr_svcluster_manta_tloc,
      runtime_override_concat_vcfs_pesr=runtime_override_concat_vcfs_pesr,
      runtime_attr_gatk_to_svtk_vcf=runtime_attr_gatk_to_svtk_vcf_pesr,
      runtime_attr_exclude_intervals_pesr=runtime_attr_exclude_intervals_pesr
  }

  output {
    File? clustered_manta_tloc_vcf = ClusterPESR.clustered_vcf
    File? clustered_manta_tloc_vcf_index = ClusterPESR.clustered_vcf_index
  }
}
