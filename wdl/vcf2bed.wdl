##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/01_pesr_clustering_single_algorithm/22/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

version 1.0

import "Structs.wdl"
import "TasksBenchmark.wdl" as tasks

workflow VcfToBed {
  input {
    File vcf
    File? vcf_idx
    String sv_base_mini_docker
  }

  if (!defined(vcf_idx)){
    call TabixVcf{
      input:
        vcf = vcf,
        sv_base_mini_docker = sv_base_mini_docker
    }
  }

  vcf_index = select_first([vcf_idx, TabixVcf.vcf_idx])

  call vcf2bed{
    input:
      vcf = vcf,
      vcf_index = vcf_index,
      sv_base_mini_docker = sv_base_mini_docker
  }

  output {
    File clustered_bed = vcf2bed.bed
  }
}


task vcf2bed {
  input {
    File vcf
    File? vcf_index
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 10, 
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  String filename = basename(vcf, ".vcf.gz")

  output {
    File bed = "${filename}.bed"
  }

  command <<<

    set -Eeuo pipefail
    svtk vcf2bed -i SVTYPE -i SVLEN ~{vcf} ~{filename}.bed

  >>>
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

task TabixVcf {
  input {
    File vcf
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output{
    File vcf_idx = "~{vcf}.idx"
  }

  command <<<
    tabix -p vcf ~{vcf}
  >>>

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




