##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/00_pesr_preprocessing_MMDLW/15/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

version 1.0

import "Structs.wdl"
# Does perlim translocation resolve from raw manta calls
workflow MakeVcfTarBall{
  input{
    Array[File] vcf_list
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr
  }

  call GenerateVcfTarBall{
    input:
      vcf_list = vcf_list,
      prefix = prefix,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr
  }

  output{
    File tar_ball = GenerateVcfTarBall.out
  }
}

task GenerateVcfTarBall {
  input {
    Array[File] vcf_list
    String prefix
    String sv_pipeline_docker
      RuntimeAttr? runtime_attr_override
  }

  Int num_samples = length(vcf_list)

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 3
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    vcfs=(~{sep=" " vcf_list})
    mkdir out
    for (( i=0; i<~{num_samples}; i++ ));
    do
      vcf=${vcfs[$i]}
      cp $vcf out/
    done
    tar czf ~{prefix}.tar.gz -C out/ .
  >>>

  output {
    File out = "~{prefix}.tar.gz"
  }
  
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


