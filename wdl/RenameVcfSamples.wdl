version 1.0

import "Structs.wdl"

workflow RenameVcfSamples {
  input {
    File vcf
    File vcf_index
    Array[String] current_sample_ids
    Array[String] new_sample_ids
    String cohort_name
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  call RenameVcfSamplesTask {
    input:
      vcf=vcf,
      vcf_idx=vcf_index,
      current_sample_list=write_lines(current_sample_ids),
      new_sample_list=write_lines(new_sample_ids),
      prefix="~{cohort_name}.cleaned.external_sample_ids",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_override
  }

  output {
    File vcf_out = RenameVcfSamplesTask.out
    File? vcf_out_index = RenameVcfSamplesTask.out_index
  }
}

task RenameVcfSamplesTask {
  input {
    File vcf
    File vcf_idx
    File current_sample_list
    File new_sample_list
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: 1.0,
    disk_gb: ceil(10 + size(vcf) * 2.0),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    paste ~{current_sample_list} ~{new_sample_list} > dict.tsv
    /opt/sv-pipeline/scripts/vcf_replacement_samples.py --vcf ~{vcf} --dict dict.tsv > reheader.list
    bcftools reheader --samples reheader.list -o ~{prefix}.vcf.gz ~{vcf}
    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
    File out_index = "~{prefix}.vcf.gz.tbi"
  }
}
