version 1.0

import "Structs.wdl"

workflow ExtractRegionFromBam {
    input {
        Array[File] bam_list          # List of BAMs
        Array[File] bai_list          # List of BAIs (same order as BAM)
        Int start
        Int end
        String chrom
        String mid_fix
        String gatk_docker
        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_tabix_bam
        RuntimeAttr? runtime_attr_extract_region
        RuntimeAttr? runtime_attr_bam_to_fastq
    }

    scatter(i in range(length(bam_list))) {
        call ExtractRegion {
            input:
                bam = bam_list[i],
                bai = bai_list[i],
                start = start,
                end = end,
                chrom = chrom,
                mid_fix = mid_fix,
                docker_image = gatk_docker,
                runtime_attr_override = runtime_attr_extract_region
        }

        call BamToFastq {
            input:
                bam = ExtractRegion.regional_bam,
                bai = ExtractRegion.regional_bai,
                docker_image = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_bam_to_fastq
         }
    }

    output {
        Array[File] regional_bams = ExtractRegion.regional_bam
        Array[File] regional_bais = ExtractRegion.regional_bai
        Array[File] regional_fastq = BamToFastq.fastq
    }
}

task ExtractRegion {
    input {
        File bam
        File bai
        String chrom
        Int start
        Int end
        String mid_fix
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bam, ".bam")

    command <<<
        set -euo pipefail

        gatk PrintReads \
           -I ~{bam} \
           -L ~{chrom}:~{start}-~{end} \
           -O "~{prefix}.~{mid_fix}.bam"

        samtools index "~{prefix}.~{mid_fix}.bam"
    >>>

    output {
        File regional_bam = "~{prefix}.~{mid_fix}.bam"
        File regional_bai = "~{prefix}.~{mid_fix}.bam.bai"
    }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(bam, "GiB")*2),
    disk_gb: 15 + ceil(size(bam, "GiB")*2),
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

task TabixBam {
    input {
        File bam
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bam, ".bam")
    command <<<
        set -euo pipefail
        gsutil cp ~{bam} ./
        samtools index ~{prefix}.bam
 
    >>>

    output {
        File regional_bam = "~{bam}.bam"
        File regional_bam_bai = "~{bam}.bam.bai"
    }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(bam, "GiB")),
    disk_gb: 15 + ceil(size(bam, "GiB")),
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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
    File? bai
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String output_prefix = basename(bam, ".bam")
  command <<<
    set -euo pipefail

    # Convert BAM to FASTQ first
    samtools fastq ~{bam} | gzip  > ~{output_prefix}.fastq.gz

  >>>

  output {
    File fastq = "~{output_prefix}.fastq.gz"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(bam, "GiB")),
    disk_gb: 15 + ceil(size(bam, "GiB")),
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

