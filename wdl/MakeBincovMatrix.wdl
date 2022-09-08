version 1.0

import "Structs.wdl"

workflow MakeBincovMatrix {
  input {
    Array[File] rd_files
    File? bincov_matrix
    File? reference_dict
    String batch
    Int? binsize
    Int? disk_overhead_gb
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  call MakeBincovMatrix {
    input:
      rd_files = rd_files,
      bincov_matrix = bincov_matrix,
      reference_dict = reference_dict,
      batch = batch,
      binsize = binsize,
      disk_overhead_gb = disk_overhead_gb,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File merged_bincov = MakeBincovMatrix.merged_bincov
    File merged_bincov_idx = MakeBincovMatrix.merged_bincov_idx
  }
}

task MakeBincovMatrix {
  input {
    Array[File] rd_files
    File? bincov_matrix
    File? reference_dict
    String batch
    Int? binsize
    Int? disk_overhead_gb
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    rd_files: {
      localization_optional: true
    }
    bincov_matrix: {
      localization_optional: true
    }
  }

  Int disk_gb = select_first([disk_overhead_gb, 10]) + ceil(size(rd_files, "GiB") + size(bincov_matrix, "GiB"))

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: disk_gb,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File merged_bincov = "~{batch}.rd.txt.gz"
    File merged_bincov_idx = "~{batch}.rd.txt.gz.tbi"
  }

  command <<<
    set -Eeu -o pipefail

    ln -s ~{write_lines(rd_files)} evidence.list
    gatk PasteDepthEvidence --java-options -Xmx3g \
      ~{"--sequence-dictionary " + reference_dict} \
      -F evidence.list ~{"-F " + bincov_matrix} \
      ~{"--bin-size " + binsize} \
      -O ~{batch}.rd.txt.gz
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + runtime_attr.disk_gb + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
