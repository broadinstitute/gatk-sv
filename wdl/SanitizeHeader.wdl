version 1.0

import "FilterGenotypes.wdl" as filter

workflow SanitizeHeaderAcrossContigs {
  input {
    Array[File] vcfs
    File contig_list

    String prefix
    String drop_fields

    String sv_pipeline_docker
  }

  Array[String] contigs = read_lines(contig_list)

  scatter (i in range(length(vcfs))) {
    call filter.SanitizeHeader {
      input:
        vcf = vcfs[i],
        vcf_index = vcfs[i] + ".tbi",
        prefix = "~{prefix}.~{contigs[i]}.sanitized",
        drop_fields=drop_fields,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  output {
    Array[File] sanitized_vcfs = SanitizeHeader.out
    Array[File] sanitized_vcf_idxs = SanitizeHeader.out_index
  }
}
