version 1.0

import "Structs.wdl"

workflow SanitizeHeader {
  input {
    File vcf
    File? vcf_index
    String prefix
    String drop_fields
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_sanitize_header
  }

  call SanitizeHeaderTask {
    input:
      vcf=vcf,
      vcf_index=select_first([vcf_index, "~{vcf}.tbi"]),
      prefix=prefix,
      drop_fields=drop_fields,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_sanitize_header
  }

  output {
    File vcf_header_sanitized = SanitizeHeaderTask.out
    File vcf_header_sanitized_index = SanitizeHeaderTask.out_index

  }
}

task SanitizeHeaderTask {
  input {
    File vcf
    File vcf_index
    String prefix
    String drop_fields
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + size(vcf, "GiB") * 2.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    bcftools view --no-version -h ~{vcf} > header.vcf

    grep -v ^"##bcftools" header.vcf | sed 's/Minimum passing GQ for each biallelic non-refvariant/Minimum passing SL for each biallelic variant/g' > newheader.vcf

    bcftools reheader -h newheader.vcf ~{vcf} \
      | bcftools annotate -x ~{drop_fields} \
        --no-version
        -O z
        -o ~{prefix}.vcf.gz

    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
    File out_index = "~{prefix}.vcf.gz.tbi"
  }
}
