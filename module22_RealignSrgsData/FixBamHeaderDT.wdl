version 1.0

import "Structs.wdl"
import "paired_fastq_to_ubam_per_readgroup.wdl" as ubam_utils

## Fixes malformed @RG fields in a BAM header: DT (e.g. "DT:HMMYHDSX7.1",
## a flowcell ID that leaked into the date field), PL, and CN. DT is set
## to the current date/time in SAM-spec ISO8601 format; PL and CN are
## both set to "ILLUMINA". Runs per uBAM in a list, then merges the
## fixed uBAMs into one.

workflow FixBamHeaderDT {
  input {
    Array[File] input_ubams
    String sample_name
    String merge_docker
    RuntimeAttr? runtime_attr_override
  }

  scatter (idx in range(length(input_ubams))) {
    call FixHeader {
      input:
        input_bam = input_ubams[idx],
        output_basename = sample_name + ".rg" + idx,
        runtime_attr_override = runtime_attr_override
    }
  }

  call ubam_utils.MergeSortedReadGroupUbams {
    input:
      sample_name = sample_name,
      input_ubams = FixHeader.fixed_bam,
      docker = merge_docker
  }

  output {
    File merged_ubam = MergeSortedReadGroupUbams.merged_ubam
    Array[File] fixed_ubams = FixHeader.fixed_bam
    Array[File] original_headers = FixHeader.original_header
    Array[File] fixed_headers = FixHeader.fixed_header
  }
}

task FixHeader {
  input {
    File input_bam
    String output_basename
    RuntimeAttr? runtime_attr_override
  }

  # Size-based defaults: reheadering is a lightweight, metadata-only
  # operation, so mem/cpu stay modest even for large BAMs. Disk needs
  # room for the input plus the reheadered copy (~2x) with headroom.
  Float input_size = size(input_bam, "GB")
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4.0,
    disk_gb: ceil(input_size * 2) + 20,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # Pull out the current header
    samtools view -H ~{input_bam} > original_header.sam

    # Current timestamp, SAM spec ISO8601 (DT tag)
    CURRENT_DATE=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    echo "Replacing DT with: ${CURRENT_DATE}"

    # Walk every @RG line and overwrite whatever DT:, PL:, and CN: fields
    # are there (malformed or not) with the current date / ILLUMINA.
    # All other RG fields (ID, SM, LB, PU, etc.) are left untouched.
    awk -v dt="$CURRENT_DATE" 'BEGIN{FS="\t"; OFS="\t"} {
      if ($1 == "@RG") {
        has_dt = 0; has_pl = 0; has_cn = 0
        for (i=1; i<=NF; i++) {
          if ($i ~ /^DT:/) { $i = "DT:" dt; has_dt = 1 }
          else if ($i ~ /^PL:/) { $i = "PL:ILLUMINA"; has_pl = 1 }
          else if ($i ~ /^CN:/) { $i = "CN:ILLUMINA"; has_cn = 1 }
        }
        if (!has_dt) { $0 = $0 "\tDT:" dt }
        if (!has_pl) { $0 = $0 "\tPL:ILLUMINA" }
        if (!has_cn) { $0 = $0 "\tCN:ILLUMINA" }
      }
      print
    }' original_header.sam > fixed_header.sam

    # Reheader in place (fast: text-only edit, no need for -P since
    # read/sort order is unaffected by a tag value change). No indexing
    # here since this is an unaligned/queryname-ordered uBAM.
    samtools reheader fixed_header.sam ~{input_bam} > ~{output_basename}.reheadered.bam
  >>>

  output {
    File fixed_bam = "~{output_basename}.reheadered.bam"
    File original_header = "original_header.sam"
    File fixed_header = "fixed_header.sam"
  }

  runtime {
    docker: "biocontainers/samtools:v1.9-4-deb_cv1"
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
