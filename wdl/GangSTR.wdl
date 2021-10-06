## Workflow to run GangSTR (https://github.com/gymreklab/GangSTR), a tool
## for computing genome-wide profile of short tandem repeats (STR) from
## short reads.

version 1.0

import "Structs.wdl"

workflow GangSTR {

  input {
    Array[File] bams_or_crams
    Array[File]? bams_or_crams_indexes
    File reference_fasta
    File? reference_fasta_index
    File regions
    String docker
    RuntimeAttr? runtime_attr
  }

  scatter (i in range(length(bams_or_crams))) {
    File bam_or_cram_ = bams_or_crams[i]
    Boolean is_bam =
      basename(bam_or_cram_, ".bam") + ".bam" == basename(bam_or_cram_)

    File bam_or_cram_index_ =
      if defined(bams_or_crams_indexes) then
        select_first([bams_or_crams_indexes])[i]
      else
        bam_or_cram_ + if is_bam then ".bai" else ".crai"
  }
  Array[File] bams_or_crams_indexes_ = select_all(bam_or_cram_index_)

  File reference_fasta_index_ = select_first([
    reference_fasta_index, reference_fasta + ".fai"])

  call CallGangSTR {
    input:
      bams_or_crams = bams_or_crams,
      bams_or_crams_indexes = bams_or_crams_indexes_,
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index_,
      regions = regions,
      docker = docker,
      runtime_attr_override = runtime_attr
  }

  output {
    File output_vcf = CallGangSTR.output_vcf
    File sample_stats = CallGangSTR.sample_stats
    File insdata = CallGangSTR.insdata
  }
}

task CallGangSTR {
  input {
    Array[File] bams_or_crams
    Array[File] bams_or_crams_indexes
    File reference_fasta
    File reference_fasta_index
    File regions
    String docker
    RuntimeAttr? runtime_attr_override
  }

  output {
    File output_vcf = "output.vcf"
    File sample_stats = "output.samplestats.tab"
    File insdata = "output.insdata.tab"
  }

  command <<<
    set -euxo pipefail

    joined_bams=(~{sep="," bams_or_crams})

    GangSTR \
      --bam $joined_bams \
      --ref ~{reference_fasta} \
      --regions ~{regions} \
      --out output
  >>>

  RuntimeAttr runtime_attr_str_profile_default = object {
    cpu_cores: 1,
    mem_gb: 4,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1,
    disk_gb: 10 + ceil(size([
      bams_or_crams,
      bams_or_crams_indexes,
      reference_fasta,
      reference_fasta_index], "GiB"))
  }
  RuntimeAttr runtime_attr = select_first([
    runtime_attr_override,
    runtime_attr_str_profile_default])

  runtime {
    docker: docker
    cpu: runtime_attr.cpu_cores
    memory: runtime_attr.mem_gb + " GiB"
    disks: "local-disk " + runtime_attr.disk_gb + " HDD"
    bootDiskSizeGb: runtime_attr.boot_disk_gb
    preemptible: runtime_attr.preemptible_tries
    maxRetries: runtime_attr.max_retries
  }
}
