##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/RawVcfQC/6/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

## Copyright Broad Institute, 2017
## 
## This WDL pipeline implements quality check on vcfs from each SV calling algorithms (https://github.com/arq5x/lumpy-sv)

## two output will be generated from this script: 
##  XXX.low describes individuals with low number of SVs on certain chromosome, indicating FC failure while running  algorithm on the sample
##  XXX.high describes individuals with high number of SVs on certain chromosome,indicating potential library fail on those samples

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

workflow RawVcfQC {
  input {
    Array[File] vcfs
    String caller
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_qc
    RuntimeAttr? runtime_attr_outlier
  }

  scatter (vcf in vcfs) {
    call RunIndividualQC {
      input:
        vcf = vcf,
        caller = caller,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_qc
    }
  }

  call PickOutliers {
    input:
      stat_files = RunIndividualQC.stat,
      caller = caller,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_outlier
  }

  output {
    File low = PickOutliers.low
    File high = PickOutliers.high
  }
}

task RunIndividualQC {
  input {
    File vcf
    String caller
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  
  String sample_name = basename(vcf, ".vcf.gz")

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 1, 
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File stat = "${caller}.${sample_name}.QC.stat"
  }
  command <<<

    python /opt/sv-pipeline/pre_SVCalling_and_QC/raw_vcf_qc/calcu_num_SVs.by_type_chromo.py ~{vcf} ~{caller}.~{sample_name}.QC.stat
    
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

task PickOutliers {
  input {
    Array[File] stat_files
    String caller
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 7.5,
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File low = "${caller}.QC.outlier.low"
    File high = "${caller}.QC.outlier.high"
  }
  command <<<

    set -euo pipefail
    cat ~{sep=" " stat_files} > ~{caller}.QC.stat
    Rscript  /opt/sv-pipeline/pre_SVCalling_and_QC/raw_vcf_qc/calcu_num_SVs.pick_outlier.R -s ~{caller}.QC.stat -o ~{caller}.QC.outlier
    
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

