version 1.0

import "CallGenomicDisorderCNVs.wdl" as gd
import "Structs.wdl"

# Runs CallGenomicDisorderCNVs on several batches: one scatter shard per batch, each with its own
# BAF matrix and high-resolution RD matrix (both tabix-indexed; the .tbi is taken from the same path).
# batches, baf_matrices and high_res_rd_matrices are parallel arrays (same length and order).
workflow CallGenomicDisorderCNVsMultiBatch {
  input {
    Array[String] batches
    Array[File] baf_matrices
    Array[File] high_res_rd_matrices

    File reference_dict
    File reference_fasta

    File gd_table
    File segdup_bed
    File centromere_bed
    File acrocentric_arm_bed
    File custom_mask_bed
    File hard_inclusion_bed
    Array[File]? flank_exclusion_intervals
    File par_bed
    File gaps_bed
    File gtf
    File? truth_table
    File ploidy_table

    Int rebinned_interval_size = 10000
    String sv_pipeline_docker
    String gatk_docker
    String? preprocess_args
    String? infer_args
    String? call_args
    String? eval_args
    String? plot_args
    RuntimeAttr? runtime_attr_condense
    RuntimeAttr? runtime_attr_gd
  }

  scatter (i in range(length(batches))) {
    call gd.CallGenomicDisorderCNVs {
      input:
        batch = batches[i],
        baf_matrix = baf_matrices[i],
        high_res_rd_matrix = high_res_rd_matrices[i],
        reference_dict = reference_dict,
        reference_fasta = reference_fasta,
        gd_table = gd_table,
        segdup_bed = segdup_bed,
        centromere_bed = centromere_bed,
        acrocentric_arm_bed = acrocentric_arm_bed,
        custom_mask_bed = custom_mask_bed,
        hard_inclusion_bed = hard_inclusion_bed,
        flank_exclusion_intervals = flank_exclusion_intervals,
        par_bed = par_bed,
        gaps_bed = gaps_bed,
        gtf = gtf,
        truth_table = truth_table,
        ploidy_table = ploidy_table,
        rebinned_interval_size = rebinned_interval_size,
        sv_pipeline_docker = sv_pipeline_docker,
        gatk_docker = gatk_docker,
        preprocess_args = preprocess_args,
        infer_args = infer_args,
        call_args = call_args,
        eval_args = eval_args,
        plot_args = plot_args,
        runtime_attr_condense = runtime_attr_condense,
        runtime_attr_gd = runtime_attr_gd
    }
  }

  output {
    Array[File] gd_output_tarballs = CallGenomicDisorderCNVs.gd_output_tarball
  }
}
