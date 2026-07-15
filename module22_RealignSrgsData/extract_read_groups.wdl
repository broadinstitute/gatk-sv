version 1.0

import "Structs.wdl"

workflow ExtractReadGroupsWorkflow {
  input {
    File fastq_gz
    String sample_id
    String docker = "ubuntu:22.04"
    RuntimeAttr? runtime_attr_extract
  }

  call ExtractReadGroups {
    input:
      fastq_gz              = fastq_gz,
      sample_id             = sample_id,
      docker                = docker,
      runtime_attr_override = runtime_attr_extract
  }

  output {
    Array[String] PU = ExtractReadGroups.PU
    Array[String] ID = ExtractReadGroups.ID
    String DT = ExtractReadGroups.DT
    String PL = ExtractReadGroups.PL
  }
}

task ExtractReadGroups {
  input {
    File fastq_gz
    String sample_id
    String docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores:        1,
    mem_gb:           ceil(2 + size(fastq_gz, "GB") * 2),
    disk_gb:          ceil(20 + size(fastq_gz, "GB") * 4),
    boot_disk_gb:     10,
    preemptible_tries: 2,
    max_retries:      1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # Extract flowcell + lane from every read header, then get unique combinations
    # Illumina header format: @<instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y> <read>:<filtered>:<control>:<index>
    zcat ~{fastq_gz} | awk 'NR % 4 == 1' | awk -F: '{print $3"."$4}' | sort -u > lane_combos.txt

    # Build PU (flowcell.lane) and ID (sample.flowcell.lane) for each unique combo
    > pu_list.txt
    > id_list.txt
    while IFS= read -r combo; do
      echo "${combo}" >> pu_list.txt
      echo "~{sample_id}.${combo}" >> id_list.txt
    done < lane_combos.txt

    # DT: date the task is run, in ISO8601 format
    date -u +"%Y-%m-%dT%H:%M:%SZ" > dt.txt
  >>>

  output {
    Array[String] PU = read_lines("pu_list.txt")
    Array[String] ID = read_lines("id_list.txt")
    String DT = read_string("dt.txt")
    String PL = "ILLUMINA"
  }

  runtime {
    cpu:           select_first([runtime_attr.cpu_cores,        default_attr.cpu_cores])
    memory:        select_first([runtime_attr.mem_gb,           default_attr.mem_gb]) + " GiB"
    disks:         "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb,   default_attr.boot_disk_gb])
    docker:        docker
    preemptible:   select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:    select_first([runtime_attr.max_retries,      default_attr.max_retries])
  }
}
