version 1.0

import "ApplyNCRAndRefArtifactFiltersPerContig.wdl" as per_contig

workflow ApplyNCRAndRefArtifactFilters {
  input {
    Array[File] vcfs
    String label = "ncr_and_refartifact"

    File? apply_filters_script

    String sv_pipeline_docker
    String sv_base_mini_docker
  }

  scatter (vcf in vcfs) {
    String base = basename(vcf, ".vcf.gz")
    call per_contig.ApplyNCRAndRefArtifactFiltersPerContig {
      input:
        vcf = vcf,
        prefix = "~{base}.~{label}",
        apply_filters_script = apply_filters_script,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_base_mini_docker = sv_base_mini_docker
    }
  }


  output {
    Array[File] filtered_vcfs = ApplyNCRAndRefArtifactFiltersPerContig.filtered_vcf
    Array[File] filtered_vcf_indexes = ApplyNCRAndRefArtifactFiltersPerContig.filtered_vcf_index
  }
}