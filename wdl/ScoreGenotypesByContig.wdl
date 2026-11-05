version 1.0

import "Structs.wdl"
import "RecalibrateGq.wdl" as recalibrate_gq
import "TasksMakeCohortVcf.wdl" as tasks_cohort
import "ScoreGenotypes.wdl" as score_genotypes

# Multi-contig version of ScoreGenotypes.
# GQ recalibration is run independently on each per-contig vcf. If truth_json is provided, the
# SL optimization tables from all contigs are pooled, and a single genome-wide set of SL cutoffs
# is optimized, matching the behavior of ScoreGenotypes on a whole-genome vcf.
workflow ScoreGenotypesByContig {

  input {
    # Per-contig vcfs, ordered to match contig_list. Each vcf must have an index at <vcf>.tbi
    Array[File] vcfs
    String output_prefix
    File contig_list

    File? truth_json
    File gq_recalibrator_model_file

    Array[String] recalibrate_gq_args = []
    Array[File] genome_tracks = []
    Float fmax_beta = 0.4
    Int optimize_vcf_records_per_shard = 50000

    String linux_docker
    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
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
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  if (defined(truth_json)) {
    scatter (i in range(length(contigs))) {
      call tasks_cohort.ScatterVcf as ScatterForOptimization {
        input:
          vcf=RecalibrateGq.filtered_vcf[i],
          records_per_shard=optimize_vcf_records_per_shard,
          prefix="~{output_prefix}.~{contigs[i]}.filter_genotypes_scatter",
          sv_pipeline_docker=sv_pipeline_docker
      }

      scatter (j in range(length(ScatterForOptimization.shards))) {
        call score_genotypes.MakeVcfTable {
          input:
            vcf=ScatterForOptimization.shards[j],
            truth_json=select_first([truth_json]),
            output_prefix="~{output_prefix}.~{contigs[i]}.vcf_table_shard_~{j}",
            sv_pipeline_docker=sv_pipeline_docker
        }
      }
    }

    call tasks_cohort.ConcatHeaderedTextFiles {
      input:
        text_files=flatten(MakeVcfTable.out),
        gzipped=true,
        output_filename="~{output_prefix}.vcf_table.tsv.gz",
        linux_docker=linux_docker
    }

    call score_genotypes.OptimizeCutoffs {
      input:
        table=ConcatHeaderedTextFiles.out,
        fmax_beta=fmax_beta,
        output_prefix="~{output_prefix}.sl_optimization",
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  output {
    Array[File] unfiltered_recalibrated_vcfs = RecalibrateGq.filtered_vcf
    Array[File] unfiltered_recalibrated_vcf_indexes = RecalibrateGq.filtered_vcf_index

    File? vcf_optimization_table = ConcatHeaderedTextFiles.out
    File? sl_cutoff_table = OptimizeCutoffs.sl_cutoff_table
    File? sl_cutoff_qc_tarball = OptimizeCutoffs.out
  }
}
