version 1.0

import "Structs.wdl"
import "FormatVcfForGatk.wdl" as fvcf
import "TasksClusterBatch.wdl" as tasks_cluster

workflow MergeBatchSites {
  input {
    String cohort             # Cohort name or project prefix for all cohort-level outputs
    File ped_file
    Array[File] depth_vcfs    # Filtered depth VCFs across batches
    Array[File] pesr_vcfs     # Filtered PESR VCFs across batches

    Array[File]? ref_pop_depth_vcfs   # Reference population depth VCFs
    Array[File]? ref_pop_pesr_vcfs   # Reference population PESR VCFs

    File contig_list
    File reference_fasta
    File reference_fasta_fai
    File reference_dict
    String? chr_x
    String? chr_y

    String sv_pipeline_docker
    String gatk_docker
    Float? java_mem_fraction
    RuntimeAttr? runtime_attr_create_ploidy
    RuntimeAttr? runtime_attr_svtk_to_gatk_vcf
    RuntimeAttr? runtime_attr_merge_pesr
    RuntimeAttr? runtime_attr_merge_depth
    RuntimeAttr? runtime_attr_gatk_to_svtk_vcf
  }

  call tasks_cluster.CreatePloidyTableFromPed {
    input:
      ped_file=ped_file,
      contig_list=contig_list,
      retain_female_chr_y=true,
      chr_x=chr_x,
      chr_y=chr_y,
      output_prefix="~{cohort}.ploidy.retain_female_chr_y",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_create_ploidy
  }

  scatter (depth_vcf in depth_vcfs) {
    call fvcf.FormatVcf as FormatDepth {
      input:
        vcf=depth_vcf,
        ploidy_table=CreatePloidyTableFromPed.out,
        args="--fix-end",
        output_prefix=basename(depth_vcf) + ".reformatted",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_svtk_to_gatk_vcf
    }
  }

  scatter (pesr_vcf in pesr_vcfs) {
    call fvcf.FormatVcf as FormatPesr {
      input:
        vcf=pesr_vcf,
        ploidy_table=CreatePloidyTableFromPed.out,
        args="--fix-end",
        output_prefix=basename(pesr_vcf) + ".reformatted",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_svtk_to_gatk_vcf
    }
  }

  call tasks_cluster.SVCluster as MergePesr {
    input:
      vcfs=flatten(select_all([FormatPesr.out, ref_pop_pesr_vcfs])),
      output_prefix="~{cohort}.merge_batch_sites.pesr_unformatted",
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
      runtime_attr_override=runtime_attr_merge_pesr
  }

  call tasks_cluster.SVCluster as MergeDepth {
    input:
      vcfs=flatten(select_all([FormatDepth.out, ref_pop_depth_vcfs])),
      output_prefix="~{cohort}.merge_batch_sites.depth_unformatted",
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

  call tasks_cluster.GatkToSvtkVcf as GatkToSvtkVcfDepth {
    input:
      vcf=MergeDepth.out,
      output_prefix="~{cohort}.merge_batch_sites.depth",
      source="depth",
      contig_list=contig_list,
      remove_formats="CN",
      remove_infos="END2",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_gatk_to_svtk_vcf
  }

  call tasks_cluster.GatkToSvtkVcf as GatkToSvtkVcfPesr {
    input:
      vcf=MergePesr.out,
      output_prefix="~{cohort}.merge_batch_sites.pesr",
      source="manta",
      contig_list=contig_list,
      remove_formats="CN",
      remove_infos="END2",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_gatk_to_svtk_vcf
  }

  output {
    File cohort_pesr_vcf = GatkToSvtkVcfPesr.out
    File cohort_pesr_vcf_index = GatkToSvtkVcfPesr.out_index
    File cohort_depth_vcf = GatkToSvtkVcfDepth.out
    File cohort_depth_vcf_index = GatkToSvtkVcfDepth.out_index
  }
}
