##########################
## EXPERIMENTAL WORKFLOW
##########################


version 1.0

import "prune_add_af.wdl" as calcAF
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "Module07DetectBatchEffects.wdl" as DetectBatchEffects

workflow DetectBatchEffectsGenomeWide {
    input{
        Array[File] vcf_list
        Array[File] vcf_idx_list
        Array[File] vcf_preMinGQ_list
        Array[File] vcf_preMinGQ_idx_list
        Array[File] contiglists
        File batches_list
        File pcrminus_batches_list
        File sample_pop_assignments
        File sample_batch_assignments
        File excludesamples_list #empty file if need be
        File famfile
        File? par_bed
        Int variants_per_shard
        String prefix

        # Optional inputs if PCR+ samples are in callset    
        File? pcrplus_batches_list
        File? pcrplus_samples_list

        String sv_base_mini_docker
        String sv_benchmark_docker
        String sv_pipeline_docker
        String sv_pipeline_updates_docker
        String sv_pipeline_base_docker

        RuntimeAttr? runtime_attr_override_concat_stat
        RuntimeAttr? runtime_attr_merge_labeled_vcfs
        RuntimeAttr? runtime_attr_override_annotate_batch_effect
        RuntimeAttr? runtime_attr_override_integrate_batch_effect_pvalues
        RuntimeAttr? runtime_attr_override_pairwise_pv_integration_PCRMINUS
        RuntimeAttr? runtime_attr_override_pairwise_pv_integration_PCRPLUS
        RuntimeAttr? runtime_attr_override_plus_minus_pv_integration
        RuntimeAttr? runtime_attr_override_one_vs_all_integration_PCRMINUS
        RuntimeAttr? runtime_attr_override_one_vs_all_integration_PCRPLUS
     }


    scatter (i in range(length(vcf_list))) {
        call DetectBatchEffects.DetectBatchEffects as DetectBatchEffects{
            input:
                vcf = vcf_list[i],
                vcf_idx = vcf_idx_list[i],
                vcf_preMinGQ = vcf_preMinGQ_list[i],
                vcf_preMinGQ_idx = vcf_preMinGQ_idx_list[i],
                batches_list = batches_list,
                pcrminus_batches_list = pcrminus_batches_list,
                sample_pop_assignments = sample_pop_assignments,
                sample_batch_assignments = sample_batch_assignments,
                excludesamples_list = excludesamples_list,
                famfile = famfile,
                contiglist = contiglists[i], 
                par_bed = par_bed,
                variants_per_shard = variants_per_shard,
                prefix = prefix,
                pcrplus_batches_list = pcrplus_batches_list,
                pcrplus_samples_list = pcrplus_samples_list,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_benchmark_docker = sv_base_mini_docker, 
                sv_pipeline_docker = sv_pipeline_docker,
                sv_pipeline_base_docker = sv_pipeline_base_docker,
                sv_pipeline_updates_docker = sv_pipeline_updates_docker,

                runtime_attr_merge_labeled_vcfs = runtime_attr_merge_labeled_vcfs,
                runtime_attr_override_pairwise_pv_integration_PCRMINUS = runtime_attr_override_pairwise_pv_integration_PCRMINUS,
                runtime_attr_override_pairwise_pv_integration_PCRPLUS = runtime_attr_override_pairwise_pv_integration_PCRPLUS,
                runtime_attr_override_plus_minus_pv_integration = runtime_attr_override_plus_minus_pv_integration,
                runtime_attr_override_one_vs_all_integration_PCRMINUS = runtime_attr_override_one_vs_all_integration_PCRMINUS,
                runtime_attr_override_one_vs_all_integration_PCRPLUS = runtime_attr_override_one_vs_all_integration_PCRPLUS
        }
    }


   # Integrate the Batch Effect statistics and annotate vcf accordingly
    call MiniTasks.ConcatStats as concat_plus_vs_minus_pv{
           input:
               shard_bed_files = DetectBatchEffects.plus_vs_minus_pv,
               prefix = "~{prefix}.plus_vs_minus",
               sv_base_mini_docker = sv_base_mini_docker,
               runtime_attr_override = runtime_attr_override_concat_stat
    }

    call MiniTasks.ConcatStats as concat_pre_vs_post_pv{
           input:
               shard_bed_files = DetectBatchEffects.pre_vs_post_pv,
               prefix = "~{prefix}.pre_vs_post",
               sv_base_mini_docker = sv_base_mini_docker,
               runtime_attr_override = runtime_attr_override_concat_stat
    }

    call MiniTasks.ConcatStats as concat_pairwise_minus_pv{
           input:
               shard_bed_files = DetectBatchEffects.pairwise_minus_pv,
               prefix = "~{prefix}.pairwise_minus",
               sv_base_mini_docker = sv_base_mini_docker,
               runtime_attr_override = runtime_attr_override_concat_stat
    }

    call MiniTasks.ConcatStats as concat_pairwise_plus_pv{
           input:
               shard_bed_files = DetectBatchEffects.pairwise_plus_pv,
               prefix = "~{prefix}.pairwise_plus",
               sv_base_mini_docker = sv_base_mini_docker,
               runtime_attr_override = runtime_attr_override_concat_stat
    }

    call MiniTasks.ConcatStats as concat_one_vs_all_minus_pv{
           input:
               shard_bed_files = DetectBatchEffects.one_vs_all_minus_pv,
               prefix = "~{prefix}.one_vs_all_minus",
               sv_base_mini_docker = sv_base_mini_docker,
               runtime_attr_override = runtime_attr_override_concat_stat
    }

    call MiniTasks.ConcatStats as concat_one_vs_all_plus_pv{
           input:
               shard_bed_files = DetectBatchEffects.one_vs_all_plus_pv,
               prefix = "~{prefix}.one_vs_all_plus",
               sv_base_mini_docker = sv_base_mini_docker,
               runtime_attr_override = runtime_attr_override_concat_stat
    }

    call IntegrateBatchEffectPvalues{
        input:
            plus_vs_minus_pv = concat_plus_vs_minus_pv.merged_stat,
            pre_vs_post_pv = concat_pre_vs_post_pv.merged_stat,
            pairwise_minus_pv = concat_pairwise_minus_pv.merged_stat,
            pairwise_plus_pv = concat_pairwise_plus_pv.merged_stat,
            one_vs_all_minus_pv = concat_one_vs_all_minus_pv.merged_stat,
            one_vs_all_plus_pv = concat_one_vs_all_plus_pv.merged_stat,
            prefix = prefix,
            sv_benchmark_docker = sv_benchmark_docker,
            runtime_attr_override = runtime_attr_override_integrate_batch_effect_pvalues
    }

    scatter (i in range(length(vcf_list))) {
        call AnnotateBatchEffect{
            input:
                vcf = vcf_list[i],
                vcf_idx = vcf_idx_list[i],
                sample_batch_assignments = sample_batch_assignments,
                SVID_info = IntegrateBatchEffectPvalues.SVID_info,
                SVID_filter = IntegrateBatchEffectPvalues.SVID_filter,
                SVID_format = IntegrateBatchEffectPvalues.SVID_format,
                sv_benchmark_docker = sv_benchmark_docker,
                runtime_attr_override = runtime_attr_override_annotate_batch_effect
        }
    }

    output{
        Array[File] annotated_vcf = AnnotateBatchEffect.anno_vcf
        Array[File] annotated_vcf_ids = AnnotateBatchEffect.anno_vcf_idx
    }
}



# Integrate batch effect statistics from different comparisons
task IntegrateBatchEffectPvalues{
  input{
      File plus_vs_minus_pv
      File pre_vs_post_pv
      File pairwise_minus_pv
      File pairwise_plus_pv
      File one_vs_all_minus_pv
      File one_vs_all_plus_pv
      String prefix
      String sv_benchmark_docker
      RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 5,
    disk_gb: 5,
    boot_disk_gb: 5,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail

    Rscript  /src/identify_SVs_failed_BatchEffect.R \
    --stat_plus_vs_minus    ~{plus_vs_minus_pv} \
    --stat_pre_vs_post      ~{pre_vs_post_pv} \
    --stat_pairwise_minus   ~{pairwise_minus_pv} \
    --stat_pairwise_plus    ~{pairwise_plus_pv} \
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
    docker: sv_benchmark_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


# Annotate vcf with Batch Effect results
task AnnotateBatchEffect{
  input{
    File vcf
    File vcf_idx
    File sample_batch_assignments
    File SVID_info
    File SVID_filter
    File SVID_format
    String sv_benchmark_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 5,
    disk_gb: ceil(10.0 + size(vcf, "GiB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  String vcf_prefix = basename(vcf, '.vcf.gz')

  command <<<
    set -euo pipefail

    python /src/annotate_vcf_with_BatchEffect.py \
    --vcf ~{vcf} \
    --SVID_info ~{SVID_info} \
    --SVID_filter ~{SVID_filter} \
    --SVID_format ~{SVID_format} \
    --sample_batch ~{sample_batch_assignments} \
    -o ~{vcf_prefix}.BatchEffect.vcf.gz

    tabix -p vcf ~{vcf_prefix}.BatchEffect.vcf.gz
  >>>

  output {
    File anno_vcf = "~{vcf_prefix}.BatchEffect.vcf.gz"
    File anno_vcf_idx = "~{vcf_prefix}.BatchEffect.vcf.gz.idx"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_benchmark_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}



 
