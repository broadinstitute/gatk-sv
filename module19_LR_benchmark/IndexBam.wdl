version 1.0

import "Structs.wdl" as tasks0506

# Index BAM and/or CRAM alignment files. Each input may be a .bam or a .cram;
# the matching index (.bai or .crai) is produced for each.
workflow IndexBams {
  input {
    Array[File] bam_files          # .bam and/or .cram files
    File reference_fasta
    File reference_index
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_index_bam
  }

  scatter (bam in bam_files) {
    call IndexBam {
      input:
        bam_file = bam,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_index_bam
    }
  }

  output {
    Array[File] bai_files = IndexBam.bam_index   # .bai for BAM inputs, .crai for CRAM inputs
  }
}

task IndexBam {
  input {
    File bam_file                  # .bam or .cram
    File reference_fasta
    File? reference_index
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    bam_file: {
      localization_optional: true
    }
  }

  # Detect input type from the extension. basename(f, ".cram") strips the suffix
  # only when present, so a changed value means the input ended in .cram.
  String input_name = basename(bam_file)
  Boolean is_cram = basename(bam_file, ".cram") != input_name
  String index_extension = if is_cram then ".crai" else ".bai"
  String index_name = input_name + index_extension

  File reference_index_file = select_first([reference_index, reference_fasta + ".fai"])

  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 4]) else 4
  Float mem_size_gb = num_cpu * 4.0

  # Indexing localizes the file as-is (no decompression), so disk tracks the
  # actual file size plus the reference and some overhead.
  Float disk_overhead = 20.0
  Float file_size = size(bam_file, "GiB")
  Float ref_size = size(reference_fasta, "GiB")
  Float ref_index_size = size(reference_index_file, "GiB")
  Int vm_disk_size = ceil(file_size + ref_size + ref_index_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: mem_size_gb,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File bam_index = index_name
  }
  command <<<

        set -Eeuo pipefail

        # localize and index the alignment file (BAM -> .bai, CRAM -> .crai)
        gsutil cp ~{bam_file} ./
        samtools index -@ ~{num_cpu} ~{input_name}
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}
