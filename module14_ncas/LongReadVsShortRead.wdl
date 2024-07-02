version 1.0

import "Structs.wdl"

workflow LongReadVsShortRead {
    input{
        File long_read_SVs
        File short_read_SVs
        File sample_list 
        File ref_header
        File query_header
        File src_tar
        String sv_base_mini_docker

    }


    call Vcf2Bed{
        input:
            vcf = long_read_SVs,
            sv_base_mini_docker = sv_base_mini_docker
    }

    Array[Array[String]] samples=read_tsv(sample_list)
    scatter(sample in samples){
        call ExtractQueryRef as extract_query_ref_LR{
            input:
                bed = Vcf2Bed.bed,
                sample = sample[0],
                ref_header = ref_header,
                query_header = query_header,
                sv_base_mini_docker = sv_base_mini_docker
        }

        call ExtractQueryRef as extract_query_ref_SR{
            input:
                bed = short_read_SVs,
                sample = sample[0],
                ref_header = ref_header,
                query_header = query_header,
                sv_base_mini_docker = sv_base_mini_docker
        }

        call SVComparison as LR_vs_SR{
            input:
                src = src_tar,
                query = extract_query_ref_LR.query,
                ref = extract_query_ref_SR.ref,
                prefix = "LR_vs_SR",
                sample = sample[0],
                sv_base_mini_docker = sv_base_mini_docker
        }

        call SVComparison as SR_vs_LR{
            input:
                src = src_tar,
                query = extract_query_ref_SR.query,
                ref = extract_query_ref_LR.ref,
                prefix = "SR_vs_LR",
                sample = sample[0],
                sv_base_mini_docker = sv_base_mini_docker
        }
    }

    output{
        Array[File] LR_vs_SR_output = LR_vs_SR.comparison
        Array[File] SR_vs_LR_output = SR_vs_LR.comparison
    }
}



task Vcf2Bed{
    input{
        File vcf
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

    output{
        File bed = "~{filebase}.bed"
    }

    String filebase = basename(SV_sites_file,".vcf.gz")

    command <<<
        set -Eeuo pipefail

        svtk vcf2bed -i SVTYPE -i SVLEN -i AF --include-filters ~{vcf} ~{filebase}.bed
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

task ExtractQueryRef{
    input{
        File bed
        File ref_header
        File query_header
        String sample
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

    output{
        File query = "~{sample}.query.gz"
        File ref = "~{sample}.ref.gz"
    }


    command <<<
        set -Eeuo pipefail

        grep ~{sample} ~{bed} | cut -f1-4,7,8 | cat ~{query_header} - | bgzip > ~{sample}.query.gz
        grep ~{sample} ~{bed} | cut -f1-4,7,8,9 | sed -e 's/$/\t~{sample}' | cat ~{ref_header} - | bgzip > ~{sample}.ref.gz

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

task SVComparison{
    input{
        File src_tar
        File query
        File ref
        String prefix
        String sample
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

    output{
        File comparison = "~{prefix}.~{sample}.bed"
    }


    command <<<
        set -Eeuo pipefail

        gsutil cp src_tar ./
        tar zxvf src.tar.gz
        bash src/compare_callsets_V2.sh -O ~{prefix}.~{sample}.bed -p ~{prefix}.~{sample} ~{query} ~{ref}
        bgzip ~{prefix}.~{sample}.bed
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















