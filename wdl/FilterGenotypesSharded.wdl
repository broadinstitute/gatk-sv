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
        sv_pipeline_docker = sv_pipeline_docker
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
    File? main_vcf_qc_tarball = MainVcfQc.sv_vcf_qc_output

    # For optional analysis
    Array[File?] vcf_optimization_table = FilterGenotypes.vcf_optimization_table
    Array[File?] sl_cutoff_qc_tarball = FilterGenotypes.sl_cutoff_qc_tarball
    Array[File] unfiltered_recalibrated_vcf = FilterGenotypes.unfiltered_recalibrated_vcf
    Array[File] unfiltered_recalibrated_vcf_index = FilterGenotypes.unfiltered_recalibrated_vcf_index
  }
}