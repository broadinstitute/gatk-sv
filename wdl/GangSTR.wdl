## Workflow to run GangSTR (https://github.com/gymreklab/GangSTR), a tool
## for computing genome-wide profile of short tandem repeats (STR) from
## short reads.

version 1.0

import "Structs.wdl"

workflow GangSTR {

  input {
    File bam_or_cram
    File? bam_or_cram_index
    File reference_fasta
    File? reference_fasta_index
    File regions
    String docker
    RuntimeAttr? runtime_attr
  }

  Boolean is_bam =
    basename(bam_or_cram, ".bam") + ".bam" == basename(bam_or_cram)

  File bam_or_cram_index_ =
    if defined(bam_or_cram_index) then
      select_first([bam_or_cram_index])
    else
      bam_or_cram + if is_bam then ".bai" else ".crai"

  File reference_fasta_index_ = select_first([
    reference_fasta_index, reference_fasta + ".fai"])

  call CallGangSTR {
    input:
      bam_or_cram = bam_or_cram,
      bam_or_cram_index = bam_or_cram_index,
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
    Array[File] bam_or_cram
    Array[File] bam_or_cram_index
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

    GangSTR \
      --bam ~{bam_or_cram} \
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
      bam_or_cram,
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
