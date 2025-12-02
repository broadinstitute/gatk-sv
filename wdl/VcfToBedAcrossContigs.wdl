version 1.0

import "VcfToBed.wdl" as bed

workflow VcfToBedAcrossContigs {
  input {
    Array[File] vcfs
    String prefix
    File primary_contigs_list
    String sv_pipeline_docker
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_svtk_vcf2bed
  }

  Array[String] contigs = read_lines(primary_contigs_list)

  scatter (i in range(length(vcfs))) {
    call bed.VcfToBed {
      input:
        vcf = vcfs[i],
        prefix = "~{prefix}.~{contigs[i]}",
        sv_pipeline_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_svtk_vcf2bed=runtime_attr_svtk_vcf2bed
    }
  }

  output {
    Array[File] beds = VcfToBed.bed
  }
}
