version 1.0

import "Structs.wdl"

workflow FilterAndMergeSnvIndelWithSvPerContig {
  input {
    # Per-contig SNV/indel VCFs and their indexes in contig order.
    Array[File] snv_indel_vcfs
    Array[File] snv_indel_vcf_indexes

    # Per-contig SV VCFs and their indexes in matching contig order.
    Array[File] sv_vcfs
    Array[File] sv_vcf_indexes

    String output_prefix = "merged_snv_indel_sv"
    String bcftools_docker = "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"

    RuntimeAttr? runtime_attr_validate
    RuntimeAttr? runtime_attr_merge
  }

  call ValidateVcfListLengths {
    input:
      snv_indel_count       = length(snv_indel_vcfs),
      sv_count              = length(sv_vcfs),
      snv_indel_index_count = length(snv_indel_vcf_indexes),
      sv_index_count        = length(sv_vcf_indexes),
      runtime_attr_override = runtime_attr_validate
  }

  scatter (idx in range(length(snv_indel_vcfs))) {
    String contig_label = basename(snv_indel_vcfs[idx], ".vcf.gz")

    call MergeOneContigSnvIndelAndSv {
      input:
        snv_indel_vcf         = snv_indel_vcfs[idx],
        snv_indel_vcf_index   = snv_indel_vcf_indexes[idx],
        sv_vcf                = sv_vcfs[idx],
        sv_vcf_index          = sv_vcf_indexes[idx],
        contig_label          = contig_label,
        output_prefix         = output_prefix,
        bcftools_docker       = bcftools_docker,
        runtime_attr_override = runtime_attr_merge
    }
  }

  output {
    File validated_inputs = ValidateVcfListLengths.validation_file
    Array[File] merged_vcfs = MergeOneContigSnvIndelAndSv.merged_vcf
    Array[File] merged_vcf_indexes = MergeOneContigSnvIndelAndSv.merged_vcf_tbi
  }
}


task ValidateVcfListLengths {
  input {
    Int snv_indel_count
    Int sv_count
    Int snv_indel_index_count
    Int sv_index_count
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    if [ "~{snv_indel_count}" -ne "~{sv_count}" ]; then
      echo "Input list length mismatch: snv_indel_vcfs=~{snv_indel_count}, sv_vcfs=~{sv_count}" >&2
      exit 1
    fi

    if [ "~{snv_indel_index_count}" -ne "~{snv_indel_count}" ]; then
      echo "Index list length mismatch: snv_indel_vcf_indexes=~{snv_indel_index_count}, snv_indel_vcfs=~{snv_indel_count}" >&2
      exit 1
    fi

    if [ "~{sv_index_count}" -ne "~{sv_count}" ]; then
      echo "Index list length mismatch: sv_vcf_indexes=~{sv_index_count}, sv_vcfs=~{sv_count}" >&2
      exit 1
    fi

    echo "ok" > validate_lists.ok
  >>>

  output {
    File validation_file = "validate_lists.ok"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: "ubuntu:22.04"
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


task MergeOneContigSnvIndelAndSv {
  input {
    File snv_indel_vcf
    File snv_indel_vcf_index
    File sv_vcf
    File sv_vcf_index
    String contig_label
    String output_prefix
    String bcftools_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 2,
    mem_gb: 32,
    disk_gb: ceil(100 + size(snv_indel_vcf, "GB") * 6 + size(sv_vcf, "GB") * 4),
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  String output_name = output_prefix + "." + contig_label + ".vcf.gz"

  command <<<
    set -euo pipefail

    # Filter SNV/indel VCF: strip records that carry SVTYPE= in their INFO field.
    # Pipe header + awk-filtered body through bcftools view -Oz to guarantee
    # well-formed BGZF output (avoids "Bad BCF record" errors from raw bgzip).
    ( bcftools view -h ~{snv_indel_vcf}; \
      bcftools view -H ~{snv_indel_vcf} \
        | awk -F'\t' '$8 !~ /(^|;)SVTYPE=/' ) \
      | bcftools view -Oz -o snv_indel.no_svtype.vcf.gz

    tabix -p vcf snv_indel.no_svtype.vcf.gz

    # Combine filtered SNV/indel records with SV records, then sort and index.
    bcftools concat -a snv_indel.no_svtype.vcf.gz ~{sv_vcf} -Ou \
      | bcftools sort -Oz -o ~{output_name}

    tabix -p vcf ~{output_name}
  >>>

  output {
    File merged_vcf = output_name
    File merged_vcf_tbi = output_name + ".tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: bcftools_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
