version 1.0

import "ApplyNCRAndRefArtifactFiltersPerContig.wdl" as per_contig

workflow ApplyNCRAndRefArtifactFilters {
  input {
    Array[File] vcfs
    File primary_contigs_list
    String cohort_id
    String label = "ncr_and_refartifact"
    File ploidy_table

    Float? no_call_rate_cutoff

    File? apply_filters_script

    String sv_pipeline_docker
    String sv_base_mini_docker
  }

  Array[String] contigs = read_lines(primary_contigs_list)

  scatter (i in range(length(vcfs))) {
    call per_contig.ApplyNCRAndRefArtifactFiltersPerContig {
      input:
        vcf = vcfs[i],
        prefix = "~{cohort_id}.~{label}.~{contigs[i]}",
        cohort_id = cohort_id,
        ploidy_table = ploidy_table,
        no_call_rate_cutoff=no_call_rate_cutoff,
        apply_filters_script = apply_filters_script,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_base_mini_docker = sv_base_mini_docker
    }
  }


  output {
    Array[File] filtered_vcfs = ApplyNCRAndRefArtifactFiltersPerContig.filtered_vcf
    Array[File] filtered_vcf_indexes = ApplyNCRAndRefArtifactFiltersPerContig.filtered_vcf_index
    Array[File] id_rename_maps = ApplyNCRAndRefArtifactFiltersPerContig.id_rename_map
  }
}