version 1.0

import "Structs.wdl"

# Filters genotypes of samples with unassigned sex for depth-only calls on chromsomes X and Y
# This is a temporary solution until #658 is addressed (https://github.com/broadinstitute/gatk-sv/issues/658)
workflow FilterUnassignedSexGenotypes {
  input {
    File vcf
    String output_prefix
    File ped_file
    String? chr_x
    String? chr_y
    File? script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  call FilterUnassignedSexGenotypesTask {
    input:
      vcf = vcf,
      output_prefix = "~{output_prefix}.filter_unassigned_sex",
      ped_file = ped_file,
      chr_x = chr_x,
      chr_y = chr_y,
      script = script,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File filter_unassigned_sex_vcf = FilterUnassignedSexGenotypesTask.out
    File filter_unassigned_sex_vcf_index = FilterUnassignedSexGenotypesTask.out_index
  }
}

# Task to run collect-pesr on a single sample
task FilterUnassignedSexGenotypesTask {
  input {
    File vcf
    String output_prefix
    File ped_file
    String? chr_x
    String? chr_y
    File? script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(50 + size(vcf, "GiB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])

  output {
    File out = "~{output_prefix}.vcf.gz"
    File out_index = "~{output_prefix}.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail
    python ~{default="/opt/sv-pipeline/scripts/filter_unassigned_sex.py" script} \
      --vcf ~{vcf} \
      --ped-file ~{ped_file} \
      ~{"--chr-x " + chr_x} \
      ~{"--chr-y " + chr_y} \
      --out ~{output_prefix}.vcf.gz
    tabix ~{output_prefix}.vcf.gz
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

