version 1.0

import "Structs.wdl"

workflow CalculateSaturationStat {
    input{
        File AoU_bed
        File AoU_bed_tbi
        File SV_count_stat
        File SVID_genomic_context
        File script_calcu_saturation
        File contig_list
        String sv_base_mini_docker

    }

    Array[Array[String]] contigs=read_tsv(contig_list)

    scatter(contig in contigs){
        call ExtractSVsPerContig{
            input:
                bed = AoU_bed,
                bed_tbi = AoU_bed_tbi,
                contig = contig[0],
                sv_base_mini_docker = sv_base_mini_docker
        }

        call CalcuSaturationStat as calculate_saturation_all_SVs{
            input:
                bed = ExtractSVsPerContig.split_bed,
                SV_count_stat = SV_count_stat,
                script_calcu_saturation = script_calcu_saturation,
                SVID_genomic_context = SVID_genomic_context,
                split_cate = "none",
                sv_base_mini_docker = sv_base_mini_docker
        }

        call CalcuSaturationStat as calculate_saturation_by_svtype{
            input:
                bed = ExtractSVsPerContig.split_bed,
                SV_count_stat = SV_count_stat,
                script_calcu_saturation = script_calcu_saturation,
                SVID_genomic_context = SVID_genomic_context,
                split_cate = "SVTYPE",
                sv_base_mini_docker = sv_base_mini_docker
        }

        call CalcuSaturationStat as calculate_saturation_by_AF{
            input:
                bed = ExtractSVsPerContig.split_bed,
                SV_count_stat = SV_count_stat,
                script_calcu_saturation = script_calcu_saturation,
                SVID_genomic_context = SVID_genomic_context,
                split_cate = "AF",
                sv_base_mini_docker = sv_base_mini_docker
        }

        call CalcuSaturationStat as calculate_saturation_by_genomic_context{
            input:
                bed = ExtractSVsPerContig.split_bed,
                SV_count_stat = SV_count_stat,
                script_calcu_saturation = script_calcu_saturation,
                SVID_genomic_context = SVID_genomic_context,
                split_cate = "GC",
                sv_base_mini_docker = sv_base_mini_docker
        }
    }


    output{
        Array[File] saturation_all = calculate_saturation_all_SVs.accumu_tables
        Array[File] saturation_by_svtype = calculate_saturation_by_svtype.accumu_tables
        Array[File] saturation_by_af = calculate_saturation_by_AF.accumu_tables
        Array[File] saturation_by_genomic_context = calculate_saturation_by_genomic_context.accumu_tables

    }
}



task ExtractSVsPerContig{
    input{
        File bed
        File bed_tbi
        String contig
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

    String filebase = basename(bed, 'bed.gz')
    output{
        File split_bed = "~{filebase}.~{contig}.bed.gz"
    }

    command <<<
        set -Eeuo pipefail

        cat <(zcat ~{bed} | head -1) <(tabix ~{bed} ~{contig}) | bgzip > ~{filebase}.~{contig}.bed.gz

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

task CalcuSaturationStat{
    input{
        File bed
        File SV_count_stat
        File script_calcu_saturation
        File SVID_genomic_context
        String split_cate
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

    String filebase = basename(bed, 'bed.gz')

    output{
        File accumu_tables = "~{filebase}.accumu_table.tar.gz"
    }

    command <<<
        set -Eeuo pipefail

        Rscript ~{script_calcu_saturation} -i ~{bed} -o ~{filebase}.accumu_table -t ~{SV_count_stat} -s ~{split_cate} --genomic_context ~{SVID_genomic_context}
        mkdir ~{filebase}.accumu_table.folder
        mv ~{filebase}.accumu_table.* ~{filebase}.accumu_table.folder
        tar czvf ~{filebase}.accumu_table.tar.gz ~{filebase}.accumu_table.folder

   >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}







