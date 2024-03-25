version 1.0

import "Structs.wdl"

workflow RemoveAnnotationsAcrossContigs {
  input {
    String prefix
    Array[File] vcfs
    Array[File] vcf_indexes
    String fields_to_remove  # supplied to bcftools annotate -x "<fields>"
    File primary_contigs_list

    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_remove_annotations
  }

  Array[String] contigs = read_lines(primary_contigs_list)

  scatter (i in range(length(vcfs))) {
    call RemoveAnnotations {
      input:
        prefix = "~{prefix}.~{contigs[i]}",
        vcf = vcfs[i],
        vcf_index = vcf_indexes[i],
        fields_to_remove = fields_to_remove,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_remove_annotations
    }
  }

  output {
    Array[File] annotations_removed_vcfs = RemoveAnnotations.annotations_removed_vcf
    Array[File] annotations_removed_vcf_indexes = RemoveAnnotations.annotations_removed_vcf_index
  }
}

task RemoveAnnotations {
  input {
    String prefix
    File vcf
    File? vcf_index
    String fields_to_remove

    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  # Disk must be scaled proportionally to the size of the VCF
  Float input_size = size(vcf, "GiB")
  RuntimeAttr default_attr = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + (input_size * 2)),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -euo pipefail

    bcftools annotate -x ~{fields_to_remove} --no-version -Oz -o ~{prefix}.vcf.gz ~{vcf}

    tabix ~{prefix}.vcf.gz

  >>>

  output {
    File annotations_removed_vcf = "~{prefix}.vcf.gz"
    File annotations_removed_vcf_index = "~{prefix}.vcf.gz.tbi"
  }

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
