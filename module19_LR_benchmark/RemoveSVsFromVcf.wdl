version 1.0

import "Structs.wdl"


workflow RemoveSVsFromVcf {
  input {
    Array[File] vcfs
    String sv_pipeline_base_docker

    RuntimeAttr? runtime_attr_remove_SV
  }

  scatter (vcf in vcfs) {
    call RemoveSVTYPE {
      input:
        vcf = vcf,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_remove_SV
    }
  }

  output {
    Array[File] filtered_vcfs = RemoveSVTYPE.out_vcf
    Array[File] filtered_idxs = RemoveSVTYPE.out_vcf_idx
  }
}

task RemoveSVTYPE {
  input {
    File vcf
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10,
    disk_gb: ceil(100 + size(vcf, "GB") *10),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  String prefix = basename(vcf, ".vcf.gz")

  command <<<
    set -euo pipefail

    # Write header without END
    bcftools view -h ~{vcf} | grep -v END > ~{prefix}.SNV_INDEL.vcf

    # Append non-SV records
    bcftools view -H -e 'INFO/SVTYPE' ~{vcf} >> ~{prefix}.SNV_INDEL.vcf
    bgzip ~{prefix}.SNV_INDEL.vcf
    tabix -p vcf ~{prefix}.SNV_INDEL.vcf.gz
  >>>

  output {
    File out_vcf = "~{prefix}.SNV_INDEL.vcf.gz"
    File out_vcf_idx = "~{prefix}.SNV_INDEL.vcf.gz.tbi"

  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

