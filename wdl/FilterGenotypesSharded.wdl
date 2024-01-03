version 1.0

import "FilterGenotypes.wdl" as filter
import "MainVcfQc.wdl" as qc

# Wrapper for FilterGenotypes for using sharded VCFs

workflow FilterGenotypesSharded {
  input {
    # Chromosome-sharded cleaned vcfs
    Array[File] vcfs
    String qc_output_prefix
    Boolean run_qc = true
  }

  scatter (vcf in vcfs) {
    call filter.FilterGenotypes {
      input:
        vcf = vcf,
        run_qc = false
    }
  }

  if (run_qc) {
    call qc.MainVcfQc {
      input:
        vcfs=FilterGenotypes.filtered_vcf,
        prefix="~{qc_output_prefix}.filter_genotypes"
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