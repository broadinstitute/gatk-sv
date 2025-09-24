version 1.0
import "Structs.wdl"

workflow FilterVCFsByAf {
    input {
        Array[File] vcfs
        String sv_base_pipeline_docker
        RuntimeAttr? runtime_attr_filter_to_01_perc
        RuntimeAttr? runtime_attr_filter_to_05_perc
        RuntimeAttr? runtime_attr_filter_to_1_perc

    }

    scatter (vcf in vcfs) {
        call FilterVcfTo01Perc { 
        	input: 
        		vcf = vcf,
        		sv_base_pipeline_docker = sv_base_pipeline_docker,
        		runtime_attr_override = runtime_attr_filter_to_01_perc
        }

        call FilterVcfTo05Perc{
        	input: 
        		vcf = vcf,
        		sv_base_pipeline_docker = sv_base_pipeline_docker,
        		runtime_attr_override = runtime_attr_filter_to_05_perc
        }


        call FilterVcfTo1Perc{
        	input: 
        		vcf = vcf,
        		sv_base_pipeline_docker = sv_base_pipeline_docker,
        		runtime_attr_override = runtime_attr_filter_to_1_perc
        }

    }

    output {
        Array[File] af001_vcfs = FilterVcfTo01Perc.af001_vcf
        Array[File] af005_vcfs = FilterVcfTo05Perc.af005_vcf
        Array[File] af01_vcfs  = FilterVcfTo1Perc.af01_vcf
        Array[File] af001_vcfs_tbi = FilterVcfTo01Perc.af001_vcf_idx
        Array[File] af005_vcfs_tbi = FilterVcfTo05Perc.af005_vcf_idx
        Array[File] af01_vcfs_tbi  = FilterVcfTo1Perc.af01_vcf_idx
    }
}

task FilterVcfTo01Perc {
    input {
        File vcf
        String sv_base_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 20, 
        disk_gb: 200,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String fname = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail

        # base filename without .vcf.gz

        # AF > 0.001
        bcftools view -i 'INFO/AF>0.001' ~{vcf} -Oz -o ~{fname}.AF001.vcf.gz
        bcftools index -t ~{fname}.AF001.vcf.gz

    >>>

    output {
        File af001_vcf = "~{fname}.AF001.vcf.gz"
        File af001_vcf_idx = "~{fname}.AF001.vcf.gz.tbi"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }    
}

task FilterVcfTo05Perc {
    input {
        File vcf
        String sv_base_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 20, 
        disk_gb: 200,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String fname = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail

        # base filename without .vcf.gz

        # AF > 0.001
        bcftools view -i 'INFO/AF>0.005' ~{vcf} -Oz -o ~{fname}.AF005.vcf.gz
        bcftools index -t ~{fname}.AF005.vcf.gz

    >>>

    output {
        File af005_vcf = "~{fname}.AF005.vcf.gz"
        File af005_vcf_idx = "~{fname}.AF005.vcf.gz.tbi"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }    
}

task FilterVcfTo1Perc {
    input {
        File vcf
        String sv_base_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 20, 
        disk_gb: 200,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String fname = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail

        # base filename without .vcf.gz

        # AF > 0.001
        bcftools view -i 'INFO/AF>0.01' ~{vcf} -Oz -o ~{fname}.AF01.vcf.gz
        bcftools index -t ~{fname}.AF01.vcf.gz

    >>>

    output {
        File af01_vcf = "~{fname}.AF01.vcf.gz"
        File af01_vcf_idx = "~{fname}.AF01.vcf.gz.tbi"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }    
}

