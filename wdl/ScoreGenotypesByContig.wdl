version 1.0

import "Structs.wdl"
import "RecalibrateGq.wdl" as recalibrate_gq

# Multi-contig version of ScoreGenotypes. Runs GQ recalibration independently on each per-contig vcf.
# SL cutoff optimization (truth_json) is not supported; run ScoreGenotypes for that.
workflow ScoreGenotypesByContig {

  input {
    # Per-contig vcfs, ordered to match contig_list. Each vcf must have an index at <vcf>.tbi
    Array[File] vcfs
    File contig_list

    File gq_recalibrator_model_file

    Array[String] recalibrate_gq_args = []
    Array[File] genome_tracks = []

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_recalibrate_gq
  }

  Array[String] contigs = read_lines(contig_list)

  scatter (i in range(length(contigs))) {
    call recalibrate_gq.RecalibrateGq {
      input:
        vcf = vcfs[i],
        vcf_index = vcfs[i] + ".tbi",
        gq_recalibrator_model_file=gq_recalibrator_model_file,
        recalibrate_gq_args=recalibrate_gq_args,
        genome_tracks=genome_tracks,
        gatk_docker=gatk_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_recalibrate_gq=runtime_attr_recalibrate_gq
    }
  }

  output {
    Array[File] unfiltered_recalibrated_vcfs = RecalibrateGq.filtered_vcf
    Array[File] unfiltered_recalibrated_vcf_indexes = RecalibrateGq.filtered_vcf_index
  }
}
