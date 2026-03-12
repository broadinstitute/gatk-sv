version 1.0

import "Structs.wdl"

workflow ExtractRegionFromBam {
    input {
        Array[File] bam_list
        Array[File?] bai_list          # Changed to optional Array
        Int start
        Int end
        String chrom
        String mid_fix
        String gatk_docker
        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_index_bam
        RuntimeAttr? runtime_attr_extract_region
        RuntimeAttr? runtime_attr_bam_to_fastq
    }

    scatter(i in range(length(bam_list))) {
        
        # Logic to index BAM if BAI is not provided
        if (!defined(bai_list[i])) {
            call IndexBam {
                input:
                    bam = bam_list[i],
                    docker_image = sv_pipeline_base_docker,
                    runtime_attr_override = runtime_attr_index_bam
            }
        }

        # Select the provided BAI or the newly generated one
        File final_bai = select_first([bai_list[i], IndexBam.bai])

        call ExtractRegion {
            input:
                bam = bam_list[i],
                bai = final_bai,
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
        Array[File] indexed_bais = select_all(IndexBam.bai) # Output newly generated bais
    }
}

task IndexBam {
    input {
        File bam
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        samtools index ~{bam} ~{basename(bam)}.bai
    >>>

    output {
        File bai = "~{basename(bam)}.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10 + ceil(size(bam, "GiB")),
        boot_disk_gb: 10,
        preemptible_tries: 3,
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

task ExtractRegion {
    input {
        File bam
        File bai
        Int start
        Int end
        String chrom
        String mid_fix
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String output_bam = basename(bam, ".bam") + "." + mid_fix + ".bam"

    command <<<
        set -euo pipefail
        # Create a localized version where BAI and BAM are in the same folder if necessary
        # However, GATK/Samtools usually handle paths if passed explicitly
        gatk PrintReads \
            -I ~{bam} \
            --read-index ~{bai} \
            -L ~{chrom}:~{start}-~{end} \
            -O ~{output_bam}
    >>>

    output {
        File regional_bam = "~{output_bam}"
        File regional_bai = "~{output_bam}.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10 + ceil(size(bam, "GiB")),
        boot_disk_gb: 10,
        preemptible_tries: 3,
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
        samtools fastq ~{bam} | bgzip > ~{output_prefix}.fastq.gz
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