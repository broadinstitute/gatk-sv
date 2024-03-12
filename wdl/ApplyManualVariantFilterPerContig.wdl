version 1.0

import "Structs.wdl"
import "ApplyManualVariantFilter.wdl" as apply

workflow ApplyManualVariantFilterPerContig {
  input {
    String prefix
    Array[File] vcfs
    Array[File] vcf_indexes
    String filter_name
    String bcftools_filter  # supplied to bcftools view -e "<filter>"
    File primary_contigs_list

    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_hard_filter_vcf
  }

  Array[String] contigs = read_lines(primary_contigs_list)

  scatter (i in range(length(vcfs))) {
    call apply.HardFilterVcf {
      input:
        prefix = "~{prefix}.~{contigs[i]}",
        vcf = vcfs[i],
        vcf_index = vcf_indexes[i],
        filter_name = filter_name,
        bcftools_filter = bcftools_filter,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_hard_filter_vcf
    }
  }

  output {
    Array[File] manual_filtered_vcf = HardFilterVcf.hard_filtered_vcf
    Array[File] manual_filtered_vcf_index = HardFilterVcf.hard_filtered_vcf_index
  }
}
