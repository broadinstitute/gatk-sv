version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks

# Consumes a set of vcfs and distributes them into contig shards

workflow ReshardVcf {
  input {
    Array[File] vcfs  # Order does not matter but must be sorted and indexed
    File contig_list
    String prefix
    Boolean? use_ssd
    String sv_base_mini_docker
    RuntimeAttr? runtime_override_reshard
  }

  Array[String] contigs = transpose(read_tsv(contig_list))[0]

  scatter (i in range(length(vcfs))) {
    File vcf_indexes = vcfs[i] + ".tbi"
  }

  scatter (contig in contigs) {
    call ReshardContig {
      input:
        vcfs=vcfs,
        vcf_indexes=vcf_indexes,
        contig=contig,
        prefix="~{prefix}.~{contig}.resharded",
        use_ssd=use_ssd,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_reshard
    }
  }
  output {
    Array[File] resharded_vcfs = ReshardContig.out
    Array[File] resharded_vcf_indexes = ReshardContig.out_index
  }
}

task ReshardContig {
  input {
    Array[File] vcfs
    Array[File] vcf_indexes
    String contig
    String prefix
    String sv_base_mini_docker
    Boolean use_ssd = false
    RuntimeAttr? runtime_attr_override
  }

  String disk_type = if use_ssd then "SSD" else "HDD"

  RuntimeAttr runtime_default = object {
                                  mem_gb: 8,
                                  disk_gb: ceil(10 + size(vcfs, "GB") * 2),
                                  cpu_cores: 4,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])
  Int n_cpu = select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
  String threads_arg = if n_cpu > 1 then "--threads " + n_cpu else ""
  runtime {
    memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, runtime_default.disk_gb]) + " " + disk_type
    cpu: n_cpu
    preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    bcftools concat ~{threads_arg} --allow-overlaps --regions "~{contig}" -Oz -o ~{prefix}.vcf.gz ~{sep=" " vcfs}
    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
    File out_index = "~{prefix}.vcf.gz.tbi"
  }
}
