version 1.0

import "Structs.wdl"

workflow CheckIdenticalVcfs {
  input {
    File left_vcf
    File right_vcf
    File check_identical_vcfs_script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  call CheckIdenticalVcfsTask {
    input:
      left_vcf = left_vcf,
      right_vcf = right_vcf,
      check_identical_vcfs_script = check_identical_vcfs_script,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File check_output = CheckIdenticalVcfsTask.check_output
  }
}

task CheckIdenticalVcfsTask {
  input {
    File left_vcf
    File right_vcf
    File check_identical_vcfs_script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: ceil(10 + size([left_vcf, right_vcf, check_identical_vcfs_script], "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File check_output = "check_identical_vcfs.output.txt"
  }

  command <<<
    set -euo pipefail
    python ~{check_identical_vcfs_script} ~{left_vcf} ~{right_vcf} > check_identical_vcfs.output.txt
  >>>

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