##########################
## EXPERIMENTAL WORKFLOW
##########################

# This experimental workflow is a helper subworkflow for the second step in the batch
# effect detection workflow. See Module07DetectBatchEffectsStep2.wdl for more info


version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as tasks

workflow AnnotateBatchEffects {
  input {
    File vcf
    File vcf_idx
    Int sv_per_shard
    File sample_batch_assignments
    File SVID_info
    File SVID_filter
    File SVID_format
    Boolean drop_empty_records = true

    String sv_pipeline_docker
    String sv_pipeline_updates_docker
    String sv_benchmark_docker

    RuntimeAttr? runtime_attr_override_shard_vcfs
    RuntimeAttr? runtime_attr_override_annotate_batch_effect
    RuntimeAttr? runtime_attr_override_combine_sharded_vcfs
  }


  # Shard input VCF for annotation
  call tasks.ScatterVcf {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      prefix=basename(vcf, ".vcf.gz") + ".sharded",
      sv_pipeline_docker=sv_pipeline_updates_docker,
      records_per_shard=sv_per_shard,
      runtime_attr_override=runtime_attr_override_shard_vcfs
  }

  # Scatter over VCF shards
  Array[Array[File]] shards_with_indexes = transpose([ScatterVcf.shards, ScatterVcf.shards_idx])
  scatter ( shard_w_index in shards_with_indexes ) {
    # Annotate batch effects per shard
    call AnnotateBatchEffect as ScatteredAnnotation {
      input:
        vcf = shard_w_index[0],
        vcf_idx = shard_w_index[1],
        sample_batch_assignments = sample_batch_assignments,
        SVID_info = SVID_info,
        SVID_filter = SVID_filter,
        SVID_format = SVID_format,
        sv_benchmark_docker = sv_benchmark_docker,
        runtime_attr_override = runtime_attr_override_annotate_batch_effect
      }
  	}

  # Merge shards into single VCF
  call CombineShardedVcfs {
    input:
      vcfs=ScatteredAnnotation.anno_vcf,
      sv_pipeline_docker=sv_pipeline_docker,
      prefix=basename(vcf, ".vcf.gz") + ".batch_fx_labeled",
      drop_empty_records=drop_empty_records,
      runtime_attr_override=runtime_attr_override_combine_sharded_vcfs
  }

  # Final output
  output {
    File annotated_vcf = CombineShardedVcfs.vcf_out
    File annotated_vcf_idx = CombineShardedVcfs.vcf_out_idx
  }
}


# Annotate vcf with Batch Effect results
task AnnotateBatchEffect {
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
  mem_gb: 8,
  disk_gb: ceil(10.0 + (3 * size(vcf, "GiB"))),
  boot_disk_gb: 10,
  preemptible_tries: 0,
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

  tabix -p vcf -f ~{vcf_prefix}.BatchEffect.vcf.gz
  >>>

  output {
  File anno_vcf = "~{vcf_prefix}.BatchEffect.vcf.gz"
  File anno_vcf_idx = "~{vcf_prefix}.BatchEffect.vcf.gz.tbi"
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



# Merge VCF shards & drop records with zero remaining non-ref alleles
task CombineShardedVcfs {
  input{
    Array[File] vcfs
    String prefix
    String sv_pipeline_docker
    Boolean? drop_empty_records
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3,
    disk_gb: 20 + (4 * ceil(size(vcfs, "GB"))),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command {
    set -euo pipefail
    vcf-concat -f ~{write_lines(vcfs)} | vcf-sort | bgzip -c > merged.vcf.gz
    if [ ~{drop_empty_records} == "true" ]; then
      bcftools +fill-tags \
        -O v merged.vcf.gz \
        -- -t AC,AN \
      | bcftools view \
        -O z \
        -o "~{prefix}.vcf.gz" \
        --exclude 'AC==0 && INFO/SVTYPE!="CNV"'
    else
      cp merged.vcf.gz "~{prefix}.vcf.gz"
    fi
    tabix -p vcf -f "~{prefix}.vcf.gz"
  }


  output {
    File vcf_out = "~{prefix}.vcf.gz"
    File vcf_out_idx = "~{prefix}.vcf.gz.tbi"
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

