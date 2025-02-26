version 1.0

import "TasksClusterBatch.wdl" as tasks_cluster

workflow RemoveExcludedRegions {
  input {
    File vcf
    File vcf_index

    Float overlap_fraction
    File intervals
    File intervals_index
    File reference_fasta_fai

    String output_prefix

    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_exclude
  }

  call tasks_cluster.ExcludeIntervalsByIntervalOverlap as ExcludeIntervals {
      input:
        vcf=vcf,
        overlap_fraction=overlap_fraction,
        intervals=intervals,
        intervals_index=intervals_index,
        reference_fasta_fai=reference_fasta_fai,
        output_prefix=output_prefix,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_exclude
  }

  output {
    File excluded_vcf = ExcludeIntervals.out
    File excluded_vcf_idx = ExcludeIntervals.out_index
  }
}
