version 1.0

# Workflow definition for Calculating Median Coverage 

import "Structs.wdl"

workflow MedianCov {
  input {
    File bincov_matrix
    String cohort_id
    Float? mem_gb_override
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr
  }

  call CalcMedCov {
    input:
      bincov_matrix = bincov_matrix,
      cohort_id = cohort_id,
      mem_gb_override = mem_gb_override,
      sv_pipeline_qc_docker = sv_pipeline_qc_docker,
      runtime_attr_override = runtime_attr
  }

  output {
    File medianCov = CalcMedCov.medianCov
  }
}

task CalcMedCov {
  input {
    File bincov_matrix
    String cohort_id
    Float? mem_gb_override
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 80,   # Requires ~0.5Gb per sample
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File medianCov = "${cohort_id}_medianCov.transposed.bed"
  }
  command <<<

    set -euo pipefail
    zcat ~{bincov_matrix} > ~{cohort_id}_fixed.bed 
    Rscript /opt/WGD/bin/medianCoverage.R ~{cohort_id}_fixed.bed -H ~{cohort_id}_medianCov.bed
    Rscript -e "x <- read.table(\"~{cohort_id}_medianCov.bed\",check.names=FALSE); xtransposed <- t(x[,c(1,2)]); write.table(xtransposed,file=\"~{cohort_id}_medianCov.transposed.bed\",sep=\"\\t\",row.names=F,col.names=F,quote=F)"
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([mem_gb_override, runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_qc_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

