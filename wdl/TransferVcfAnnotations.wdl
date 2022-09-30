version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as tasks
import "Utils.wdl" as utils

workflow TransferVcfAnnotations {
  input {
    File vcf
    File vcf_with_annotations
    Int shard_size
    String? info_keys_list
    String? format_keys_list
    File? transfer_script
    String prefix
    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_shard
    RuntimeAttr? runtime_override_transfer
    RuntimeAttr? runtime_override_concat
  }

  call utils.ShardVcfPair {
    input:
      vcf_a=vcf,
      vcf_b=vcf_with_annotations,
      shard_size=shard_size,
      prefix_a=basename(vcf, ".vcf.gz"),
      prefix_b=basename(vcf_with_annotations, ".vcf.gz"),
      drop_a=false,
      drop_b=true,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_shard
  }

  scatter (i in range(length(ShardVcfPair.shards_a))) {
    call TransferVcfAnnotationsTask {
      input:
        vcf_to_annotate=ShardVcfPair.shards_a[i],
        vcf_with_annotations=ShardVcfPair.shards_b[i],
        info_keys_list=info_keys_list,
        format_keys_list=format_keys_list,
        script=transfer_script,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_transfer
    }
  }

  call tasks.ConcatVcfs {
    input:
      vcfs=TransferVcfAnnotationsTask.annotated_shard,
      vcfs_idx=TransferVcfAnnotationsTask.annotated_shard_index,
      naive=true,
      generate_index=true,
      outfile_prefix=prefix,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat
  }

  output {
    File annotated_vcf = ConcatVcfs.concat_vcf
    File annotated_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}


task TransferVcfAnnotationsTask {
  input {
    File vcf_to_annotate
    File vcf_with_annotations
    String? info_keys_list  # At least one of info_keys_list or format_keys_list must be provided
    String? format_keys_list
    File? script  # For debugging
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  String output_file_name = sub(sub(basename(vcf_to_annotate), ".gz$", ""), ".vcf$", "_annotated.vcf.gz")

  # Disk must be scaled proportionally to the size of the VCF
  # Memory may need to be increased as well, particularly if transferring FORMAT fields on large VCFs
  RuntimeAttr default_attr = object {
                               mem_gb: 7.5,
                               disk_gb: ceil(100.0 + size(vcf_to_annotate, "GB") * 2 + size(vcf_with_annotations, "GB")),
                               cpu_cores: 1,
                               preemptible_tries: 3,
                               max_retries: 1,
                               boot_disk_gb: 10
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -euo pipefail
    python ~{default="/opt/sv-pipeline/scripts/transfer_vcf_annotations.py" script} \
      ~{"--infos " + info_keys_list} \
      ~{"--formats " + format_keys_list} \
      --ann-vcf ~{vcf_with_annotations} \
      --out ~{output_file_name} \
      ~{vcf_to_annotate}
    tabix ~{output_file_name}
  >>>

  output {
    File annotated_shard = output_file_name
    File annotated_shard_index = output_file_name + ".tbi"
  }
}
