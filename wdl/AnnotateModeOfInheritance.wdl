version 1.0

import "Structs.wdl"

# Annotate mode of inheritance (MOI) for trio de novo callsets.
#
# Compares the case sample's genotype at every variant with the (optional)
# mother and father genotypes and adds two INFO fields to every record:
#
#   MOI:            DE_NOVO | INHERITED_FROM_MOTHER | INHERITED_FROM_FATHER |
#                   INHERITED_FROM_BOTH | PARENT_ONLY | UNASSESSABLE
#   MOI_CONFIDENCE: CONFIRMED | UNCONFIRMED
#
# UNCONFIRMED indicates that one parent was not assayed (or its genotype was
# missing) at that variant, so the MOI call cannot be fully confirmed
# (e.g. a DE_NOVO call in a half-trio).
task AnnotateModeOfInheritance {
  input {
    File vcf          # bgzipped VCF
    File vcf_idx      # tabix index at <vcf>.tbi
    String prefix     # output prefix; <prefix>.vcf.gz / <prefix>.vcf.gz.tbi / <prefix>.moi_summary.tsv
    String case_sample
    String? mother_sample
    String? father_sample
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(20 + size(vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{prefix}.vcf.gz"
    File out_index = "~{prefix}.vcf.gz.tbi"
    File moi_summary = "~{prefix}.moi_summary.tsv"
  }

  command <<<
    set -euo pipefail

    python /opt/sv-pipeline/05_annotation/scripts/annotate_moi.py \
      ~{vcf} \
      ~{prefix} \
      --case ~{case_sample} \
      ~{if defined(mother_sample) then "--mother ~{mother_sample}" else ""} \
      ~{if defined(father_sample) then "--father ~{father_sample}" else ""}
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
