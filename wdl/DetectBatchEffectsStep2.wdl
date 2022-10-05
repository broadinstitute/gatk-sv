##########################
## EXPERIMENTAL WORKFLOW
##########################

# This experimental workflow is the second in a two-part process to identify and label
# variants discovered by GATK-SV that exhibit evidence of allele frequency distortions
# specific to an individual batch (or subset of batches) of samples


version 1.0

import "TasksMakeCohortVcf.wdl" as MiniTasks
import "DetectBatchEffectsStep2Helper.wdl" as Step2Helper

workflow DetectBatchEffectsStep2 {
  input {
    Array[File] plus_vs_minus_pvalues
    Array[File] pre_vs_post_pvalues
    Array[File] pairwise_minus_pvalues
    Array[File] pairwise_plus_pvalues
    Array[File] one_vs_all_minus_pvalues
    Array[File] one_vs_all_plus_pvalues
    Array[File] vcfs
    Int sv_per_shard_annotation
    File sample_batch_assignments
    String prefix

    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_override_shard_vcfs
    RuntimeAttr? runtime_attr_override_annotate_batch_effect
    RuntimeAttr? runtime_attr_override_combine_sharded_vcfs
   }

  # Concatenate batch effect test statistics for all input shards
  call MiniTasks.ConcatStats as concat_plus_vs_minus_pv {
       input:
         shard_bed_files = plus_vs_minus_pvalues,
         prefix = "~{prefix}.plus_vs_minus",
         sv_base_mini_docker = sv_base_mini_docker
  }
  call MiniTasks.ConcatStats as concat_pre_vs_post_pv {
       input:
         shard_bed_files = pre_vs_post_pvalues,
         prefix = "~{prefix}.pre_vs_post",
         sv_base_mini_docker = sv_base_mini_docker
  }
  call MiniTasks.ConcatStats as concat_pairwise_minus_pv {
       input:
         shard_bed_files = pairwise_minus_pvalues,
         prefix = "~{prefix}.pairwise_minus",
         sv_base_mini_docker = sv_base_mini_docker
  }
  call MiniTasks.ConcatStats as concat_pairwise_plus_pv {
       input:
         shard_bed_files = pairwise_plus_pvalues,
         prefix = "~{prefix}.pairwise_plus",
         sv_base_mini_docker = sv_base_mini_docker
  }
  call MiniTasks.ConcatStats as concat_one_vs_all_minus_pv {
       input:
         shard_bed_files = one_vs_all_minus_pvalues,
         prefix = "~{prefix}.one_vs_all_minus",
         sv_base_mini_docker = sv_base_mini_docker
  }
  call MiniTasks.ConcatStats as concat_one_vs_all_plus_pv {
       input:
         shard_bed_files = one_vs_all_plus_pvalues,
         prefix = "~{prefix}.one_vs_all_plus",
         sv_base_mini_docker = sv_base_mini_docker
  }

  # Integrate results from different batch effect tests
  call IntegrateBatchEffectPvalues {
    input:
      plus_vs_minus_pv = concat_plus_vs_minus_pv.merged_stat,
      pre_vs_post_pv = concat_pre_vs_post_pv.merged_stat,
      pairwise_minus_pv = concat_pairwise_minus_pv.merged_stat,
      pairwise_plus_pv = concat_pairwise_plus_pv.merged_stat,
      one_vs_all_minus_pv = concat_one_vs_all_minus_pv.merged_stat,
      one_vs_all_plus_pv = concat_one_vs_all_plus_pv.merged_stat,
      prefix = prefix,
      sv_pipeline_docker = sv_pipeline_docker
  }

  scatter (vcf in vcfs) {
    call Step2Helper.AnnotateBatchEffects as ScatteredAnnotation {
      input:
        vcf = vcf,
        vcf_idx = vcf + ".tbi",
        sv_per_shard = sv_per_shard_annotation,
        sample_batch_assignments = sample_batch_assignments,
        SVID_info = IntegrateBatchEffectPvalues.SVID_info,
        SVID_filter = IntegrateBatchEffectPvalues.SVID_filter,
        SVID_format = IntegrateBatchEffectPvalues.SVID_format,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override_shard_vcfs = runtime_attr_override_shard_vcfs,
        runtime_attr_override_annotate_batch_effect = runtime_attr_override_annotate_batch_effect,
        runtime_attr_override_combine_sharded_vcfs = runtime_attr_override_combine_sharded_vcfs
    }
  }

  output{
    Array[File] annotated_vcfs = ScatteredAnnotation.annotated_vcf
    Array[File] annotated_vcf_idxs = ScatteredAnnotation.annotated_vcf_idx
  }
}


# Integrate batch effect statistics from different comparisons
task IntegrateBatchEffectPvalues {
  input{
    File plus_vs_minus_pv
    File pre_vs_post_pv
    File pairwise_minus_pv
    File pairwise_plus_pv
    File one_vs_all_minus_pv
    File one_vs_all_plus_pv
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 5,
    disk_gb: 10 + (5 * ceil(size([plus_vs_minus_pv, pre_vs_post_pv, pairwise_minus_pv, pairwise_plus_pv, one_vs_all_minus_pv, one_vs_all_plus_pv], "GB"))),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  
  command <<<
    set -euo pipefail

    Rscript /opt/sv-pipeline/scripts/benchmark/identify_SVs_failed_BatchEffect.R \
      --stat_plus_vs_minus  ~{plus_vs_minus_pv} \
      --stat_pre_vs_post    ~{pre_vs_post_pv} \
      --stat_pairwise_minus   ~{pairwise_minus_pv} \
      --stat_pairwise_plus  ~{pairwise_plus_pv} \
      --stat_one_vs_all_minus ~{one_vs_all_minus_pv} \
      --stat_one_vs_all_plus  ~{one_vs_all_plus_pv} \
      --out_fail_info   ~{prefix}.SVID_failed.info_col.tsv \
      --out_fail_filter ~{prefix}.SVID_failed.filter_col.tsv \
      --out_fail_format ~{prefix}.SVID_failed.format_col.tsv
  >>>

  output {
    File SVID_info = "~{prefix}.SVID_failed.info_col.tsv"
    File SVID_filter = "~{prefix}.SVID_failed.filter_col.tsv"
    File SVID_format = "~{prefix}.SVID_failed.format_col.tsv"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
