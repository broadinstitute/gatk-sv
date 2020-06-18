version 1.0

import "Structs.wdl"

workflow MergeCohortVcfs {
  input {
    Array[File] beds    # Filtered depth VCFs across batches
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_merge_list
  }

  call MergeList {
    input:
      regeno_beds = beds,
      prefix = "master_regeno",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_merge_list
  }

  output {
    File regeno = MergeList.master_regeno
  }
}

task MergeList{
    input{
    String prefix
    Array[File] regeno_beds
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command{
        cat ~{sep=' ' regeno_beds} |sort -k1,1V -k2,2n -k3,3n > ~{prefix}.bed
    }
    output{
        File master_regeno="master_regeno.bed"
    }
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
