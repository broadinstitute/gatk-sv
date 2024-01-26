version 1.0

import "Structs.wdl"

workflow RemoveINSChr2End2Fields {
  input {
    Array[String] prefixes
    Array[File] vcfs

    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_remove_fields
  }

  scatter (i in range(length(vcfs))) {
    call RemoveFields {
      input:
        prefix = prefixes[i],
        vcf = vcfs[i],
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_remove_fields
    }
  }

  output {
    Array[File] revised_vcfs = RemoveFields.revised_vcf
    Array[File] revised_vcf_indexes = RemoveFields.revised_vcf_index
  }
}

task RemoveFields {
  input {
    String prefix
    File vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + (size(vcf, "GB") * 2)),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -euo pipefail

    python <<CODE
import pysam
with pysam.VariantFile("~{vcf}") as vcf_in, pysam.VariantFile("~{prefix}.vcf.gz", mode='w', header=vcf_in.header) as vcf_out:
  for record in vcf_in:
    if record.info['SVTYPE'] == 'INS':
      if 'CHR2' in record.info:
        record.info.pop('CHR2')
      if 'END2' in record.info:
        record.info.pop('END2')
    vcf_out.write(record)
CODE

    tabix ~{prefix}.vcf.gz

  >>>

  output {
    File revised_vcf = "~{prefix}.vcf.gz"
    File revised_vcf_index = "~{prefix}.vcf.gz.tbi"
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
