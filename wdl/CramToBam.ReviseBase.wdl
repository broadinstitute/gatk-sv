version 1.0

import "Structs.wdl"

workflow CramToBamReviseBase {
  input {
    File cram_file
    File? cram_index  # required if cram is requester pays
    File reference_fasta
    File? reference_index
    File contiglist
    String samtools_cloud_docker
    RuntimeAttr? runtime_attr_split_cram
    RuntimeAttr? runtime_attr_revise_base
    RuntimeAttr? runtime_attr_concat_bam
  }

  Array[String] contigs = transpose(read_tsv(contiglist))[0]
  String base_name = basename(cram_file, ".cram")

  scatter (contig in contigs) {
    call SplitCramPerContig {
      input:
        cram_file = cram_file,
        contig = contig,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        samtools_cloud_docker = samtools_cloud_docker,
        runtime_attr_override = runtime_attr_split_cram
    }
    call ReviseBaseInBam {
      input:
        bam_file = SplitCramPerContig.bam_file,
        bam_index = SplitCramPerContig.bam_index,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        samtools_cloud_docker = samtools_cloud_docker,
        runtime_attr_override = runtime_attr_revise_base
    }
  }
 
  call ConcatBam {
    input:
      prefix = base_name,
      bam_files = ReviseBaseInBam.revised_bam_file,
      bam_indexes = ReviseBaseInBam.revised_bam_index, 
      samtools_cloud_docker = samtools_cloud_docker,
      runtime_attr_override = runtime_attr_concat_bam
  }

  output {
    File bam_file = ConcatBam.bam_file
    File bam_index = ConcatBam.bam_index
  }
}

task SplitCramPerContig {
  input {
    File cram_file
    File reference_fasta
    File? reference_index
    String contig
    String samtools_cloud_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    cram_file: {
      localization_optional: true
    }
  }

  String bam_file_name = basename(cram_file, ".cram")

  File reference_index_file = select_first([reference_index, reference_fasta + ".fai"])

  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 4]) else 4
  
  Float cram_inflate_ratio = 3.5
  Float disk_overhead = 10.0
  Float cram_size = size(cram_file, "GiB")
  Float bam_size = (cram_inflate_ratio * cram_size) / 10
  Float ref_size = size(reference_fasta, "GiB")
  Float ref_index_size = size(reference_index_file, "GiB")
  Int vm_disk_size = ceil(bam_size + ref_size + ref_index_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: 1.5, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File bam_file  = "~{bam_file_name}.~{contig}.bam"
    File bam_index = "~{bam_file_name}.~{contig}.bam.bai"
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
                 -o "~{bam_file_name}.~{contig}.bam" \
                 "~{cram_file}" \
                 "~{contig}"
        
        # index bam file
        samtools index -@ ~{num_cpu} "~{bam_file_name}.~{contig}.bam"    
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

task ReviseBaseInBam {
  input {
    File bam_file
    File bam_index
    File reference_fasta
    File? reference_index
    String samtools_cloud_docker
    RuntimeAttr? runtime_attr_override
  }

  String base_name = basename(bam_file, ".bam") 

  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 4]) else 4

  File reference_index_file = select_first([reference_index, reference_fasta + ".fai"])
  Float cram_inflate_ratio = 3.5
  Float disk_overhead = 30.0
  Float cram_size = size(bam_file, "GiB")
  Float bam_size = cram_inflate_ratio * cram_size
  Float ref_size = size(reference_fasta, "GiB")
  Float ref_index_size = size(reference_index_file, "GiB")
  Int vm_disk_size = ceil(cram_size + bam_size + ref_size + ref_index_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: 1.5,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    samtools view -H ~{bam_file} > ~{base_name}.revised.sam

    paste \
      <(samtools view ~{bam_file} | cut -f1-9) \
      <(samtools view ~{bam_file} | cut -f10 | sed -e "s/Y/N/g" | sed -e "s/R/N/g" | sed -e "s/W/N/g" | sed -e "s/S/N/g" | sed -e "s/K/N/g" | sed -e "s/M/N/g" | sed -e "s/D/N/g" | sed -e "s/H/N/g" | sed -e "s/V/N/g" | sed -e "s/B/N/g" | sed -e "s/X/N/g" ) \
      <(samtools view ~{bam_file} | cut -f11-) \
      >> ~{base_name}.revised.sam

    samtools view -Sb ~{base_name}.revised.sam -o ~{base_name}.revised.bam

    samtools index ~{base_name}.revised.bam

  >>>

  output{
    File revised_bam_file = "~{base_name}.revised.bam"
    File revised_bam_index = "~{base_name}.revised.bam.bai"
  }

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

task ConcatBam {
  input {
    String prefix
    Array[File] bam_files
    Array[File] bam_indexes
    String samtools_cloud_docker
    RuntimeAttr? runtime_attr_override
  }

  Int num_cpu = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 2]) else 2
  
  Float disk_overhead = 10.0
  Float bam_size = size(bam_files, "GiB")
  Float index_size = size(bam_indexes, "GiB")
  Int vm_disk_size = 2 * ceil(bam_size + index_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu,
    mem_gb: 1.5,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    samtools merge ~{prefix}.bam ~{sep=" "  bam_files} 
    samtools index ~{prefix}.bam
  >>>

  output{
    File bam_file = "~{prefix}.bam"
    File bam_index = "~{prefix}.bam.bai"
  }

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




