version 1.0

import "Structs.wdl"

# A wrapper to minimap2 for mapping & aligning (groups of) sequences to a reference
task Minimap2 {
    input {
        File reads
        File ref_fasta
        String map_preset
        String prefix = "out"
        String long_read_align_docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        reads:      "query sequences to be mapped and aligned"
        ref_fasta:  "reference fasta"
        map_preset: "preset to be used for minimap2 parameter '-x'"
        prefix:     "[default-valued] prefix for output BAM"
    }

    Int disk_size = 1 + 3*ceil(size(reads, "GB") + size(ref_fasta, "GB"))

    command <<<
        set -euxo pipefail

        mem=$(grep '^MemTotal' /proc/meminfo | awk '{ print int($2/1000000) }')
        cpus=$(grep -c '^processor' /proc/cpuinfo | awk '{ print $1 }')

        MAP_PARAMS="-ayL --MD -x ~{map_preset} -t $cpus ~{ref_fasta}"
        SORT_PARAMS="-@ $cpus -m ${mem}G --no-PG -o ~{prefix}.bam"

        if [[ "~{reads}" =~ \.bam$ ]]; then
            samtools fastq ~{reads} > reads.fq
            FILE="reads.fq"
        else
            FILE="~{reads}"
        fi

        minimap2 $MAP_PARAMS ~{reads} | samtools sort $SORT_PARAMS

        samtools index ~{prefix}.bam
    >>>

    output {
        File aligned_bam = "~{prefix}.bam"
        File aligned_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             30,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 long_read_align_docker
    }

}

task Minimap2_simple {
    input {
        File reads
        File ref_fasta

        String map_preset

        String prefix = "out"
        String long_read_align_docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        reads:      "query sequences to be mapped and aligned"
        ref_fasta:  "reference fasta"
        map_preset: "preset to be used for minimap2 parameter '-x'"
        prefix:     "[default-valued] prefix for output BAM"
    }

    Int disk_size = 1 + 3*ceil(size(reads, "GB") + size(ref_fasta, "GB"))

    Int cpus = 4
    Int mem = 30

    command <<<
        set -euxo pipefail

        minimap2 -ayYL --MD -x ~{map_preset} -t ~{cpus} ~{ref_fasta} ~{reads} | samtools sort -@~{cpus} -m~{mem}G --no-PG -o ~{prefix}.bam
        samtools index ~{prefix}.bam
    >>>

    output {
        File aligned_bam = "~{prefix}.bam"
        File aligned_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             mem,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 long_read_align_docker
    }
}

# A simple task to covert SAM-formatted alignment to PAF format
task SAMtoPAF {
    input {
        File sam_formatted_file
        File? index

        RuntimeAttr? runtime_attr_override
        String long_read_align_docker
    }

    parameter_meta {
        sam_formatted_file: "SAM-formated input file to be converted to PAF (note currently we only support SAM or BAM, not CRAM)"
        index:              "[optional] index for sam_formatted_file"
    }

    String prefix = basename(basename(sam_formatted_file, ".bam"), ".sam") # we have hack like this because WDL stdlib doesn't provide endsWith stuff

    Int disk_size = 2*ceil(size(sam_formatted_file, "GB"))

    command <<<
        set -eu

        filename=$(basename -- ~{sam_formatted_file})
        extension="${filename##*.}"
        if [[ "$extension" == "sam" ]]; then
            /minimap2-2.17_x64-linux/k8 \
                /minimap2-2.17_x64-linux/paftools.js \
                sam2paf \
                -L \
                ~{sam_formatted_file} \
                > ~{prefix}".paf"
        elif [[ "$extension" == "bam" ]]; then
            samtools view -h ~{sam_formatted_file} | \
                /minimap2-2.17_x64-linux/k8 \
                /minimap2-2.17_x64-linux/paftools.js \
                sam2paf \
                -L \
                - \
                > ~{prefix}".paf"
        else
            echo "Currently we only support SAM or BAM (not CRAM)." && exit 1;
        fi
    >>>

    output {
        File pat_formatted_file = "~{prefix}.paf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 long_read_align_docker
    }
}
