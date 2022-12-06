version 1.0

import "Structs.wdl"

workflow RemoveDuplicateEvents {
  input {
    File vcf
    File? vcf_index
    String prefix

    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_remove_duplicate_events_task
  }

  call RemoveDuplicateEventsTask {
    input:
      vcf = vcf,
      vcf_index = select_first([vcf_index, "~{vcf}.tbi"]),
      prefix = prefix,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_remove_duplicate_events_task
  }

  output {
    File deduplicated_vcf = RemoveDuplicateEventsTask.deduplicated_vcf
    File deduplicated_vcf_index = RemoveDuplicateEventsTask.deduplicated_vcf_index
    File duplicated_events_table = RemoveDuplicateEventsTask.duplicated_events_table
  }
}


task RemoveDuplicateEventsTask {
  input {
    File vcf
    File vcf_index
    String prefix

    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_vcf = "~{prefix}.duplicates_removed.vcf.gz"
  String output_table = "~{prefix}.duplicated_events.tsv"

  # Disk must be scaled proportionally to the size of the VCF
  Float input_size = size(vcf, "GiB")
  RuntimeAttr default_attr = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + (input_size * 2)),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    /opt/sv-pipeline/scripts/remove_duplicates.py ~{vcf} -o ~{output_vcf} -t ~{output_table}

    tabix ~{output_vcf}
  >>>

  output {
    File deduplicated_vcf = "~{output_vcf}"
    File deduplicated_vcf_index = "~{output_vcf}.tbi"
    File duplicated_events_table = "~{output_table}"
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
