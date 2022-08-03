## Copyright Broad Institute, 2022
## 
##
## Consolidate boost scores per sample across all batches and write those scores
## directly into an input VCF
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


workflow FilterVcfWithRandomForestScores {
  input {
    File vcf
    File vcf_idx
    File cutoff_table
    String sv_benchmark_docker
    RuntimeAttr? runtime_attr_override_filter_vcf
  }

  # Scatter over tarballs of boost scores


  # Column-wise merge of all annotated VCFs
  call FilterVcf {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      cutoff_table=cutoff_table,
      sv_benchmark_docker=sv_benchmark_docker,
      runtime_attr_override = runtime_attr_override_filter_vcf
  }

  output {
    File annotated_vcf = FilterVcf.filtered_vcf
    File annotated_vcf_idx = FilterVcf.filtered_vcf_idx
  }
}


task FilterVcf {
  input {
    File vcf
    File vcf_idx
    File cutoff_table
    String sv_benchmark_docker
    RuntimeAttr? runtime_attr_override
  }


  Float input_size = size(vcf, "GB")
  Float base_disk_gb = 10.0
  RuntimeAttr runtime_default = object {
      mem_gb: 1,
      disk_gb: ceil(base_disk_gb + (input_size * 10.0)),
      cpu_cores: 1,
      preemptible_tries: 3,
      max_retries: 1,
      boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
      memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
      disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
      cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
      preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
      maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
      docker: sv_benchmark_docker
      bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  String prefix = basename(vcf, '.vcf.gz')
  command <<<
    set -eu -o pipefail

    python /src/filter_vcf_with_RF.py \
    ~{vcf} \
    ~{prefix}.filtered.vcf.gz \
    --stats ~{cutoff_table} 

    tabix -p vcf ~{prefix}.filtered.vcf.gz

  >>>

  output {
    File filtered_vcf = "~{prefix}.filtered.vcf.gz"
    File filtered_vcf_idx = "~{prefix}.filtered.vcf.gz.tbi"
  }
}


