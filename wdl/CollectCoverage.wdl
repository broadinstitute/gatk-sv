version 1.0

import "Structs.wdl"

task CollectCounts {
  input {
    File intervals
    File bam
    File bam_idx
    String sample_id
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File? gatk4_jar_override
    Array[String]? disabled_read_filters

    # Runtime parameters
    String gatk_docker
    Float? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts
  }

  parameter_meta {
    bam: {
      localization_optional: true
    }
    bam_idx: {
      localization_optional: true
    }
  }

  Float mem_overhead_gb = 2.0
  Float machine_mem_gb = select_first([mem_gb, 12.0])
  Int command_mem_mb = floor((machine_mem_gb - mem_overhead_gb) * 1024)
  Array[String] disabled_read_filters_arr = if(defined(disabled_read_filters))
    then
      prefix(
        "--disable-read-filter ",
        select_first([disabled_read_filters])
      )
    else
      []

  command <<<
    set -euo pipefail
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

    gatk --java-options "-Xmx~{command_mem_mb}m" CollectReadCounts \
      -L ~{intervals} \
      --input ~{bam} \
      --read-index ~{bam_idx} \
      --reference ~{ref_fasta} \
      --format TSV \
      --interval-merging-rule OVERLAPPING_ONLY \
      --output ~{sample_id}.counts.tsv \
      ~{sep=' ' disabled_read_filters_arr}

    sed -ri "s/@RG\tID:GATKCopyNumber\tSM:.+/@RG\tID:GATKCopyNumber\tSM:~{sample_id}/g" ~{sample_id}.counts.tsv
    bgzip ~{sample_id}.counts.tsv
  >>>

  runtime {
    docker: gatk_docker
    memory: machine_mem_gb + " GiB"
    disks: "local-disk " + select_first([disk_space_gb, 50]) + if use_ssd then " SSD" else " HDD"
    cpu: select_first([cpu, 1])
    preemptible: select_first([preemptible_attempts, 5])
    maxRetries: 1
  }

  output {
    File counts = "~{sample_id}.counts.tsv.gz"
  }
}

task CondenseReadCounts {
  input {
    File counts
    String sample
    Int? num_bins
    Int? expected_bin_size
    File? gatk4_jar_override

    # Runtime parameters
    String condense_counts_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1.0,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float machine_mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int command_mem_mb = floor(machine_mem_gb*1000) - 500

  command <<<
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}
    gunzip -c ~{counts} > counts.tsv
    gatk --java-options "-Xmx~{command_mem_mb}m" CondenseReadCounts \
      -I counts.tsv \
      -O condensed_counts.~{sample}.tsv \
      --factor ~{select_first([num_bins, 20])} \
      --out-bin-length ~{select_first([expected_bin_size, 2000])}
    sed -ri "s/^@RG\tID:GATKCopyNumber\tSM:.+/@RG\tID:GATKCopyNumber\tSM:~{sample}/g" condensed_counts.~{sample}.tsv
    bgzip condensed_counts.~{sample}.tsv
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: machine_mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: condense_counts_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  output {
    File out = "condensed_counts.~{sample}.tsv.gz"
  }
}

task CountsToIntervals {
  input {
    File counts
    String output_name

    # Runtime parameters
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1.0,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    zgrep "^@" ~{counts} > ~{output_name}.interval_list
    zgrep -v "^@" ~{counts} | sed -e 1d | awk -F "\t" -v OFS="\t" '{print $1,$2,$3,"+","."}' >> ~{output_name}.interval_list
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  output {
    File out = "~{output_name}.interval_list"
  }
}
