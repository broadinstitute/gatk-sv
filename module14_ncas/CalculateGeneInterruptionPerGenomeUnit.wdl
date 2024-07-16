version 1.0

import "Structs.wdl"

workflow CalculateGeneInterruptionPerGenomeUnit {
    input{
        File AoU_bed
        File src
        Array[File] sample_list
        String sv_base_mini_docker

    }

    scatter(sample_file in sample_list){
        call CalcuSvGenePerGenome{
            input:
                bed = AoU_bed,
                src = src,
                sample_file = sample_file,
                sv_base_mini_docker = sv_base_mini_docker
        }
    }

    String prefix = basename(AoU_bed, '.bed.gz')
    
    call ConcatFile{
        input:
            file_list = ExtractSVsPerContig.output,
            prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker
    }


    output{
        File SV_gene_sample_stat = ConcatFile.Concat_file
    }
}



task CalcuSvGenePerGenome{
    input{
        File bed
        File src
        File sample_file
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String filebase = basename(bed, '.bed.gz')
    output{
        File output = "~{filebase}.tsv"
    }

    command <<<
        set -Eeuo pipefail

        gsutil cp ~{src} ./
        tar zxvf src.tar.gz

        Rscript ./src/calcu.SV_gene_interruption_per_genome.R -i ~{bed} -o ~{filebase}.tsv -s ~{sample_file}

   >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }    }


task ConcatFile {
    input {
        Array[File] file_list
        String prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    } 
    
    
    # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
    # be held in memory or disk while working, potentially in a form that takes up more space)
    Float input_size = size(file_list, "GB")
    RuntimeAttr runtime_default = object {
        mem_gb: 2.0,
        disk_gb: ceil(10.0 + input_size * 7.0),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -eux

        set -o pipefail

        while read SPLIT; do
          zcat $SPLIT
        done < ~{write_lines(file_list)} \
          | (grep -Ev "^#" || printf "") \
          | bgzip -c \
          > ~{prefix}.gz

    >>>

  output {
    File Concat_file = "~{prefix}.gz"
  }
}






