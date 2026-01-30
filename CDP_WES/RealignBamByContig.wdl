version 1.0

workflow RealignBamByContig {
  input {
    File input_bam
    File contig_list
    File input_bai
    File reference_fasta
    File reference_fai
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
	        contig = DetectContigs.contigs[i],
	        docker_image = sv_pipeline_base_docker
  	}

  	call SplitRefToContig {
  		input:
  			contig = contig_list[i],
	        reference_fasta = reference_fasta,
	        reference_fai = reference_fai,
	        docker_image = sv_pipeline_base_docker
  	}

    call RealignOneContig {
      input:
        fq_1 = BamToFastq.fq_1,
        fq_2 = BamToFastq.fq_2,
        contig_fa = SplitRefToContig.ref_fa,
        contig_fai = SplitRefToContig.ref_fai,
        docker_image = bwa_docker
    }
  }

  call MergeBams {
    input:
      bams = RealignOneContig.out_bam,
      docker_image = sv_pipeline_base_docker
  }

  output {
    File realigned_bam = MergeBams.merged_bam
    File realigned_bai = MergeBams.merged_bai
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

    # Align
    bwa mem ~{contig_fa} ~{fq_1} ~{fq_2} | \
      samtools sort -o realigned.bam

    samtools index realigned.bam
  >>>

  output {
    File out_bam = "realigned.bam"
    File out_bai = "realigned.bam.bai"
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
