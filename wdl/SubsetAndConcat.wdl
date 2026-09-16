version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as TasksMakeCohortVcf

workflow SubsetAndConcat {
  input {
    Array[File] vcfs
    File contigs_list
    String prefix
    String bcftools_view_options
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_override_subset_vcf
    RuntimeAttr? runtime_override_concat_vcfs
    RuntimeAttr? runtime_override_concat_sites_only_vcfs
  }

  Array[String] contigs = read_lines(contigs_list)

  scatter (i in range(length(vcfs))) {
    call SubsetVcf {
      input:
        vcf=vcfs[i],
        contig=contigs[i],
        prefix=prefix,
        bcftools_view_options=bcftools_view_options,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_subset_vcf
    }
  }

  call TasksMakeCohortVcf.ConcatVcfs {
    input:
      vcfs=SubsetVcf.outvcf,
      vcfs_idx=SubsetVcf.outvcf_index,
      sites_only=false,
      outfile_prefix="~{prefix}",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_vcfs
  }

  call TasksMakeCohortVcf.ConcatVcfs as ConcatSitesOnlyVcfs {
    input:
      vcfs=SubsetVcf.outvcf,
      vcfs_idx=SubsetVcf.outvcf_index,
      sites_only=true,
      outfile_prefix="~{prefix}.sites_only",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_sites_only_vcfs
  }

  output {
    Array[File] subset_vcfs = SubsetVcf.outvcf
    Array[File] subset_vcf_indexes = SubsetVcf.outvcf_index
    File concat_vcf = ConcatVcfs.concat_vcf
    File concat_vcf_index = ConcatVcfs.concat_vcf_idx
    File sites_only_concat_vcf = ConcatSitesOnlyVcfs.concat_vcf
    File sites_only_concat_vcf_index = ConcatSitesOnlyVcfs.concat_vcf_idx
  }
}

task SubsetVcf {
  input {
    File vcf
    String contig
    String prefix
    String bcftools_view_options
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_prefix = prefix + "." + contig
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: ceil(10 + 2.0 * size(vcf, "GiB")),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  command <<<
    set -euo pipefail

    bcftools view \
      --no-version \
      ~{bcftools_view_options} \
      ~{vcf} \
      | bcftools +fill-tags -- -t AC,AF,AN \
      | bgzip -c \
      > ~{output_prefix}.vcf.gz
    tabix ~{output_prefix}.vcf.gz
  >>>

  output {
    File outvcf = "~{output_prefix}.vcf.gz"
    File outvcf_index = "~{output_prefix}.vcf.gz.tbi"
  }

  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }
}
