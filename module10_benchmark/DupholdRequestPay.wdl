##########################################################################################


## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

## Copyright Broad Institute, 2020
## 
## This WDL pipeline implements Duphold 
##
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

version 1.0

import "Structs.wdl"
import "TasksBenchmark.wdl" as mini_tasks
workflow DupholdRequestPay{
  input{
    String prefix
    String bam_or_cram_file
    String bam_or_cram_index
    File vcf_file
    File ref_fasta
    File ref_fai
    File ref_dict
    File contig_list
    String duphold_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_duphold
    RuntimeAttr? runtime_attr_bcf2vcf
    RuntimeAttr? runtime_attr_LocalizeCram
    RuntimeAttr? runtime_attr_SplitVcf
    RuntimeAttr? runtime_attr_ConcatVcfs
  }

  call RunDuphold{
    input:
      prefix = prefix,
      bam_or_cram_file =  bam_or_cram_file,
      bam_or_cram_index = bam_or_cram_index,
      vcf_file = vcf_file,
      ref_fasta = ref_fasta,
      ref_fai = ref_fai,
      ref_dict = ref_dict,
      duphold_docker = duphold_docker,
      runtime_attr_override = runtime_attr_duphold
    }

  call Bcf2Vcf{
      input:
        prefix = prefix,
        bcf = RunDuphold.bcf,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_bcf2vcf
      }

    output{
      File vcf = Bcf2Vcf.vcf
      File vcf_idx = Bcf2Vcf.vcf_idx
    }
  }

task RunDupholdPerContig{
  input{
    String prefix
    String contig
    File bam_or_cram_file
    File bam_or_cram_index
    File vcf_file
    File vcf_index
    File ref_fasta
    File ref_fai
    File ref_dict
    String duphold_docker
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


  output {
    File bcf = "~{prefix}.~{contig}.bcf"
  }
  command <<<

    set -Eeuo pipefail
    
    duphold -t 4 \
    -v ~{vcf_file} \
    -b ~{bam_or_cram_file} \
    -f ~{ref_fasta} \
    -o ~{prefix}.~{contig}.bcf

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: duphold_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task Bcf2Vcf{
  input{
    String prefix
    File bcf
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 5,
    boot_disk_gb: 5,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output{
    File vcf = "~{prefix}.duphold.vcf.gz"
    File vcf_idx = "~{prefix}.duphold.vcf.gz.tbi"
  }
  command <<<
      set -Eeuo pipefail
      bcftools view ~{bcf} | bgzip > ~{prefix}.duphold.vcf.gz
      tabix -p vcf ~{prefix}.duphold.vcf.gz
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

task RunDuphold{
  input{
    String prefix
    File bam_or_cram_file
    File bam_or_cram_index
    File vcf_file
    File ref_fasta
    File ref_fai
    File ref_dict
    String duphold_docker
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


  output {
    File bcf = "~{prefix}.bcf"
  }
  command <<<

    set -Eeuo pipefail
    
    duphold -t 4 \
    -v ~{vcf_file} \
    -b ~{bam_or_cram_file} \
    -f ~{ref_fasta} \
    -o ~{prefix}.bcf

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: duphold_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}





