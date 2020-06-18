version 1.0

import "Structs.wdl"

workflow CramToBam {
  input {
    File cram_file
    File reference_fasta
    File? reference_index
    String samtools_cloud_docker
    Boolean requester_pays = false
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
      cram_file: ".bam or .cram file to search for SVs. bams are preferable, crams will be converted to bams."
      reference_fasta: ".fasta file with reference used to align bam or cram file"
      reference_index: "[optional] reference index file. If omitted, the WDL will look for an index by appending .fai to the .fasta file"
  }
  meta {
      author: "Ted Brookings"
      email: "tbrookin@broadinstitute.org"
  }

  if (requester_pays) {
    call RunCramToBamRequesterPays {
      input:
        cram_file = cram_file,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        samtools_cloud_docker = samtools_cloud_docker,
        runtime_attr_override = runtime_attr_override
    }
  }
  if (!requester_pays) {
    call RunCramToBam {
      input:
        cram_file = cram_file,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        samtools_cloud_docker = samtools_cloud_docker,
        runtime_attr_override = runtime_attr_override
    }
  }

  output {
    File bam_file = select_first([RunCramToBamRequesterPays.bam_file, RunCramToBam.bam_file])
    File bam_index = select_first([RunCramToBamRequesterPays.bam_index, RunCramToBam.bam_index])
  }
}

task RunCramToBam {
  input {
    File cram_file
    File reference_fasta
    File? reference_index
    String samtools_cloud_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    cram_file: {
      localization_optional: true
    }
  }

  String bam_file_name = basename(cram_file, ".cram") + ".bam"

  File reference_index_file = select_first([reference_index, reference_fasta + ".fai"])

  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 4]) else 4
  Float mem_size_gb = num_cpu * 4.0
  
  Float cram_inflate_ratio = 3.5
  Float disk_overhead = 20.0
  Float cram_size = size(cram_file, "GiB")
  Float bam_size = cram_inflate_ratio * cram_size
  Float ref_size = size(reference_fasta, "GiB")
  Float ref_index_size = size(reference_index_file, "GiB")
  Int vm_disk_size = ceil(bam_size + ref_size + ref_index_size + disk_overhead)

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
    File bam_file = bam_file_name
    File bam_index = bam_file_name + ".bai"
  }
  command <<<

        set -Eeuo pipefail

        # necessary for getting permission to read from google bucket directly
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
  
        # covert cram to bam
        samtools view \
                 -b \
                 -h \
                 -@ ~{num_cpu} \
                 -T "~{reference_fasta}" \
                 -o "~{bam_file_name}" \
                 "~{cram_file}"
        
        # index bam file
        samtools index -@ ~{num_cpu} "~{bam_file_name}"    
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: samtools_cloud_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

# This version localizes the cram since samtools does not support streaming from requester pays buckets
task RunCramToBamRequesterPays {
  input {
    File cram_file
    File reference_fasta
    File? reference_index
    String samtools_cloud_docker
    RuntimeAttr? runtime_attr_override
  }

  String bam_file_name = basename(cram_file, ".cram") + ".bam"

  File reference_index_file = select_first([reference_index, reference_fasta + ".fai"])

  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 4]) else 4
  Float mem_size_gb = num_cpu * 4.0

  Float cram_inflate_ratio = 3.5
  Float disk_overhead = 20.0
  Float cram_size = size(cram_file, "GiB")
  Float bam_size = cram_inflate_ratio * cram_size
  Float ref_size = size(reference_fasta, "GiB")
  Float ref_index_size = size(reference_index_file, "GiB")
  Int vm_disk_size = ceil(cram_size + bam_size + ref_size + ref_index_size + disk_overhead)

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
    File bam_file = bam_file_name
    File bam_index = bam_file_name + ".bai"
  }
  command <<<

        set -Eeuo pipefail

        # covert cram to bam
        samtools view \
                 -b \
                 -h \
                 -@ ~{num_cpu} \
                 -T "~{reference_fasta}" \
                 -o "~{bam_file_name}" \
                 "~{cram_file}"

        # index bam file
        samtools index -@ ~{num_cpu} "~{bam_file_name}"
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: samtools_cloud_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}
