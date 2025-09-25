version 1.0
import "Structs.wdl"

workflow VcfLiftOver {
    input {
        File vcf
        File vcf_idx
        String sv_base_mini_docker
        String liftover_docker
        File update_vcf_py
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_liftover
        RuntimeAttr? runtime_attr_update_vcf

    }

    call Vcf2Bed {
        input: 
            vcf = vcf,
            vcf_idx = vcf_idx,
            docker_file = liftover_docker,
            runtime_attr_override = runtime_attr_vcf2bed
        }

    call LiftOver{
        input:
            bed = Vcf2Bed.bed,
            docker_file = liftover_docker,
            runtime_attr_override = runtime_attr_liftover
    }

    call UpdateVcf {
        input: 
            bed = LiftOver.bed_hg38, 
            vcf = vcf,
            vcf_idx = vcf_idx,
            update_vcf_py = update_vcf_py,
            docker_file = liftover_docker,
            runtime_attr_override = runtime_attr_update_vcf
    }


    output {
        File lifted_vcf = UpdateVcf.updated_vcf
        File lifted_vcf_idx = UpdateVcf.updated_vcf_tbi
    }
}

task Vcf2Bed {
    input {
        File vcf
        File vcf_idx
        String docker_file
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
        python /opt/xz_scripts/vcf2bed.mc_vcf.py ~{vcf} ~{prefix}.bed
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
        String docker_file
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
        /opt/liftOver ~{bed} /opt/chain_data/hs1ToHg38.over.chain.gz ~{prefix}.hg38.bed ~{prefix}.hg38.remain
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
        File update_vcf_py
        String docker_file
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
        python ~{update_vcf_py} ~{bed} ~{vcf} ~{prefix}.hg38.vcf.gz
        gunzip ~{prefix}.hg38.vcf.gz
        bgzip ~{prefix}.hg38.vcf
        bcftools sort ~{prefix}.hg38.vcf.gz -Oz -o ~{prefix}.hg38.sorted.vcf.gz
        tabix -p vcf ~{prefix}.hg38.sorted.vcf.gz
    >>>

    output {
        File updated_vcf = "~{prefix}.hg38.sorted.vcf.gz"
        File updated_vcf_tbi = "~{prefix}.hg38.sorted.vcf.gz.tbi"
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
