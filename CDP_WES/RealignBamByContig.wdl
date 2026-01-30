version 1.0

import "Structs.wdl"

workflow RealignBamByContig {
  input {
    File input_bam
    File input_bai
    Array[String] contig_list
    File reference_fasta
    File reference_fai
    String sv_base_mini_docker
    String sv_pipeline_base_docker
    String bwa_docker
  }

  call DetectContigs {
    input:
      bam = input_bam,
      docker_image = sv_pipeline_base_docker
  }

  scatter (i in range(length(contig_list))){
    call BamToFastq{
      input:
          bam = input_bam,
          bai = input_bai,
          contig = DetectContigs.contigs[i],
          docker_image = sv_pipeline_base_docker
    }
  }

  call ConcatAndBgzipFastq as concat_r1 {
    input:
      fastq_files = BamToFastq.fq_1,
      output_prefix = 'fq_1',
      docker_image = sv_base_mini_docker
    }

  call ConcatAndBgzipFastq as concat_r2 {
    input:
      fastq_files = BamToFastq.fq_2,
      output_prefix = 'fq_2',
      docker_image = sv_base_mini_docker
    }

    call RealignOneContig {
      input:
        fq_1 = concat_r1.fastq_gz,
        fq_2 = concat_r2.fastq_gz,
        contig_fa = reference_fasta,
        contig_fai = reference_fai,
        docker_image = sv_base_mini_docker
   }


  output {
    File realigned_bam = RealignOneContig.out_bam
    File realigned_bai = RealignOneContig.out_bai
  }
}

task DetectContigs {
  input {
    File bam
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  # Define runtime attributes
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: ceil(10 + size(bam, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    samtools idxstats ~{bam} | awk '$3+$4>0 {print $1}' > contigs.txt
  >>>

  output {
    Array[String] contigs = read_lines("contigs.txt")
  }

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

task BamToFastq {
  input {
    File bam
    File bai
    String contig
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  # Define runtime attributes
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: ceil(30 + size(bam, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # Extract reads for contig
    samtools view -b ~{bam} "~{contig}" > contig.bam

    # Convert to paired FASTQ
    samtools fastq \
      -1 r1.fq.gz \
      -2 r2.fq.gz \
      -0 /dev/null \
      -s /dev/null \
      -n \
      contig.bam
  >>>

  output {
    File fq_1 = "r1.fq.gz"
    File fq_2 = "r2.fq.gz"
  }

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

task ConcatAndBgzipFastq {
  input {
    Array[File] fastq_files   # List of input FASTQ files
    String output_prefix      # Prefix for the output file (without extension)
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  # Default runtime attributes
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 5,
    disk_gb: ceil(30 + size(fastq_files, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # Concatenate all FASTQ files into one
    cat ~{sep=' ' fastq_files} | bgzip -c > ~{output_prefix}.fastq.gz
  >>>

  output {
    File fastq_gz = "~{output_prefix}.fastq.gz"
  }

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

task SplitRefToContig {
  input {
    String contig
    File reference_fasta
    File reference_fai
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  # Define runtime attributes
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: ceil(10 + size(reference_fasta, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # Extract contig FASTA
    samtools faidx ~{reference_fasta} "~{contig}" > contig.fa
    samtools faidx contig.fa
  >>>

  output {
    File ref_fa = "contig.fa"
    File ref_fai = "contig.fa.fai"
  }

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

task RealignOneContig {
  input {
    File fq_1
    File fq_2
    File contig_fa
    File contig_fai
    String prefix
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  # Define runtime attributes
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: ceil(20 + size(fq_1, "GB") * 5),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    bwa index ~{contig_fa}

    # Align
    bwa mem ~{contig_fa} ~{fq_1} ~{fq_2} | \
      samtools sort -o "~{prefix}.realigned.bam"

    samtools index "~{prefix}.realigned.bam"
  >>>

  output {
    File out_bam = "~{prefix}.realigned.bam"
    File out_bai = "~{prefix}.realigned.bam.bai"
  }

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

task MergeBams {
  input {
    Array[File] bams
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  # Define runtime attributes
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: ceil(20 + size(bams, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    samtools merge merged.bam ~{sep=' ' bams}
    samtools index merged.bam
  >>>

  output {
    File merged_bam = "merged.bam"
    File merged_bai = "merged.bam.bai"
  }

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
