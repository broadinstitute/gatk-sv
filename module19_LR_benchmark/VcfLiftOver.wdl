version 1.0
import "Structs.wdl"

workflow VcfLiftOver {
    input {
        Array[File] vcfs
        Array[File] vcf_idxes
        File chain_file   # hs1ToHg38.over.chain.gz
        File vcf2bed_script
        File UpdateVcfWithBed_script
        String sv_base_mini_docker
        String liftover_docker
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_liftover
        RuntimeAttr? runtime_attr_update_vcf

    }

    scatter(i in range(length(vcfs))){
        call Vcf2Bed {
            input: 
                vcf = vcfs[i],
                vcf_idx = vcf_idxes[i],
                vcf2bed_script = vcf2bed_script,
                docker_file = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_vcf2bed
            }

        call LiftOver{
            input:
                bed = Vcf2Bed.bed,
                chain = chain_file,
                docker_file = liftover_docker,
                runtime_attr_override = runtime_attr_liftover
        }

        call UpdateVcf {
            input: 
                bed = LiftOver.bed_hg38, 
                vcf = vcfs[i],
                vcf_idx = vcf_idxes[i],
                UpdateVcfWithBed_script = UpdateVcfWithBed_script
                docker_file = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_update_vcf
        }


    }


    output {
        Array[File] lifted_vcf = UpdateVcf.updated_vcf
        Array[File] lifted_vcf_idx = UpdateVcf.updated_vcf_tbi
    }
}

task Vcf2Bed {
    input {
        File vcf
        File vcf_idx
        File vcf2bed_script
        File docker_file
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 10, 
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }


    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(vcf, ".vcf.gz")


    command <<<
        python ~{vcf2bed_script} ~{vcf} ~{prefix}.bed
    >>>

    output {
        File bed = "~{prefix}.bed"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_file
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }    
}

task LiftOver {
    input {
        File bed
        File chain
        File docker_file
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 10, 
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }


    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(bed, ".bed")

    command <<<
        ./src/liftOver ~{bed} ~{chain} ~{prefix}.hg38.bed ~{prefix}.hg38.remain
    >>>

    output {
        File bed_hg38 = "~{prefix}.hg38.bed"
        File remain   = "~{prefix}.hg38.remain"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_file
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }    
}

task UpdateVcf {
    input {
        File bed
        File vcf
        File vcf_idx
        File UpdateVcfWithBed_script
        File docker_file
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 10, 
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }


    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(vcf, ".vcf.gz")

    command <<<
        python ~{UpdateVcfWithBed_script} ~{bed} ~{vcf} ~{prefix}.hg38.vcf.gz
        tabix -p vcf ~{prefix}.hg38.vcf.gz
    >>>

    output {
        File updated_vcf = "~{prefix}.hg38.vcf.gz"
        File updated_vcf_tbi = "~{prefix}.hg38.vcf.gz.tbi"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_file
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }    
}
