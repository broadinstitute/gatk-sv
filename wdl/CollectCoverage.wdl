version 1.0

import "Structs.wdl"

task CondenseReadCounts {
  input {
    File rd_file
    File ref_dict
    Int? max_interval_size
    Int? min_interval_size
    File? gatk4_jar_override

    # Runtime parameters
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    rd_file: {
      localization_optional: true
    }
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.0,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

    ln -s /dev/stdout tmp.rd.txt
    gatk --java-options -Xmx2g CondenseDepthEvidence -F ~{rd_file} \
        --sequence-dictionary ~{ref_dict} \
        --max-interval-size ~{default=2000 max_interval_size} \
        --min-interval-size ~{default=101 min_interval_size} \
        -O tmp.rd.txt \
      | awk 'BEGIN{FS=OFS="\t"}
         {if(NR==1){print "@RG\tID:GATKCopyNumber\tSM:" $4; print "CONTIG\tSTART\tEND\tCOUNT"}
          else print $1, $2+1, $3, $4}' \
      | cat ~{ref_dict} - \
      | bgzip > condensed_counts.tsv.gz
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  output {
    File counts = "condensed_counts.tsv.gz"
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
    File interval_list = "~{output_name}.interval_list"
  }
}
