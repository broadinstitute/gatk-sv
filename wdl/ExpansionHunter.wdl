## Runs ExpansionHunter denovo

version 1.0

import "Structs.wdl"

workflow ComputeSTRProfiles {

  input {
    Array[File] case_reads_filenames
    Array[File] control_reads_filenames
    File manifest_filename
    File reference_filename
    String output_prefix
    Int min_anchor_mapq
    Int max_irr_mapq
    String ehdn_docker
    RuntimeAttr? runtime_attr
  }

  parameter_meta {
    case_reads_filenames: ""
    control_reads_filenames: ""
    manifest_filename: ""
    reference_filename: ""
    output_prefix: ""
    min_anchor_mapq: ""
    max_irr_mapq: ""
    ehdn_docker: ""
    runtime_attr: ""
  }

  scatter(reads_filename in case_reads_filenames) {
    call ComputeSTRProfile {
      input:
        reads_filename = reads_filename,
        reference_filename = reference_filename,
        output_prefix = output_prefix,
        min_anchor_mapq = min_anchor_mapq,
        max_irr_mapq = max_irr_mapq,
        ehdn_docker = ehdn_docker,
        runtime_attr_override = runtime_attr
    }
  }

  call Merge {
    input:
      reference_filename = reference_filename,
      manifest_filename = manifest_filename,
      output_prefix = output_prefix,
      ehdn_docker = ehdn_docker
  }

  output {
    Array[File] locus = ComputeSTRProfile.locus
    Array[File] motif = ComputeSTRProfile.motif
    Array[File] profile = ComputeSTRProfile.profile
    File multisample_profile = Merge.multisample_profile
  }
}

task ComputeSTRProfile {
  input {
    File reads_filename
    File reference_filename
    String output_prefix
    Int min_anchor_mapq
    Int max_irr_mapq
    String ehdn_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File locus = "${output_prefix}.locus.tsv"
    File motif = "${output_prefix}.motif.tsv"
    File profile = "${output_prefix}.str_profile.json"
  }
  command <<<

    ExpansionHunterDenovo profile \
      --reads ~{reads_filename} \
      -reference ~{reference_filename} \
      --output-prefix ~{output_prefix} \
      --min-anchor-mapq ~{min_anchor_mapq} \
      --max-irr-mapq ~{max_irr_mapq}
  >>>

  runtime {
    cpu: runtime_attr.cpu_cores
    memory: runtime_attr.mem_gb + " GiB"
    disks: "local-disk " + runtime_attr.disk_gb + " HDD"
    bootDiskSizeGb: runtime_attr.boot_disk_gb
    preemptible: runtime_attr.preemptible_tries
    maxRetries: runtime_attr.max_retries
    docker: ehdn_docker
  }
}

task Merge {
  input {
    File reference_filename
    File manifest_filename
    String output_prefix
    String ehdn_docker
  }

  output {
    File multisample_profile = "${output_prefix.multisample_profile.json}"
  }

  command <<<
    ExpansionHunterDenovo merge \
      --reference ~{reference_filename} \
      --manifest ~{manifest_filename} \
      --output-prefix ~{output_prefix}
  >>>

  runtime {
    docker: ehdn_docker
  }
}
