version 1.0

import "FilterGenotypes.wdl" as filter
import "MainVcfQc.wdl" as qc

# Wrapper for FilterGenotypes for using sharded VCFs

workflow FilterGenotypesSharded {
  input {
    # VCFs sharded by position
    Array[File] vcfs
    String qc_output_prefix
    File primary_contigs_fai

    # For MainVcfQc
    Boolean run_qc = true
    String qc_bcftools_preprocessing_options = "-e 'FILTER~\"UNRESOLVED\" || FILTER~\"HIGH_NCR\"'"
    Int qc_sv_per_shard = 2500
    Int qc_samples_per_shard = 600

    String linux_docker
    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_recalibrate_scatter
    RuntimeAttr? runtime_attr_recalibrate_gq
    RuntimeAttr? runtime_attr_recalibrate_concat

    RuntimeAttr? runtime_attr_scatter_for_optim
    RuntimeAttr? runtime_attr_make_vcf_table
    RuntimeAttr? runtime_attr_merge_tables
    RuntimeAttr? runtime_attr_optim_cutoffs

    RuntimeAttr? runtime_attr_scatter_for_filter
    RuntimeAttr? runtime_attr_filter
    RuntimeAttr? runtime_attr_filter_concat
  }

  scatter (vcf in vcfs) {
    call filter.FilterGenotypes {
      input:
        vcf = vcf,
        run_qc = false,
        primary_contigs_fai = primary_contigs_fai,
        linux_docker = linux_docker,
        gatk_docker = gatk_docker,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_recalibrate_scatter=runtime_attr_recalibrate_scatter,
        runtime_attr_recalibrate_gq=runtime_attr_recalibrate_gq,
        runtime_attr_recalibrate_concat=runtime_attr_recalibrate_concat,
        runtime_attr_scatter_for_optim=runtime_attr_scatter_for_optim,
        runtime_attr_make_vcf_table=runtime_attr_make_vcf_table,
        runtime_attr_merge_tables=runtime_attr_merge_tables,
        runtime_attr_optim_cutoffs=runtime_attr_optim_cutoffs,
        runtime_attr_scatter_for_filter=runtime_attr_scatter_for_filter,
        runtime_attr_filter=runtime_attr_filter,
        runtime_attr_filter_concat=runtime_attr_filter_concat
    }
  }

  if (run_qc) {
    call qc.MainVcfQc {
      input:
        vcfs=FilterGenotypes.filtered_vcf,
        prefix="~{qc_output_prefix}.filter_genotypes",
        bcftools_preprocessing_options=qc_bcftools_preprocessing_options,
        primary_contigs_fai=primary_contigs_fai,
        sv_per_shard=qc_sv_per_shard,
        samples_per_shard=qc_samples_per_shard,
        sv_pipeline_qc_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  output {
    Array[File] filtered_vcf = FilterGenotypes.filtered_vcf
    Array[File] filtered_vcf_index = FilterGenotypes.filtered_vcf_index
    File? filter_genotypes_main_vcf_qc_tarball = MainVcfQc.sv_vcf_qc_output

    # For optional analysis
    Array[File?] vcf_optimization_table = FilterGenotypes.vcf_optimization_table
    Array[File?] sl_cutoff_qc_tarball = FilterGenotypes.sl_cutoff_qc_tarball
    Array[File] unfiltered_recalibrated_vcf = FilterGenotypes.unfiltered_recalibrated_vcf
    Array[File] unfiltered_recalibrated_vcf_index = FilterGenotypes.unfiltered_recalibrated_vcf_index
  }
}