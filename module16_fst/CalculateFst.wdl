version 1.0

import "Structs.wdl"

workflow CalculateFst {
    input{
        File vcf
        File vcf_idx
        File samp_pop
        String variant_type 
        String sv_fst_docker
        RuntimeAttr? runtime_attr_fst
    }

    if (variant_type=="SV"){
        call CalcuFstSv{
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                samp_pop = samp_pop,
                sv_fst_docker = sv_fst_docker,
                runtime_attr_override = runtime_attr_fst
            }
        }

    if (variant_type=="SNV"){
        call CalcuFstSnv{
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                samp_pop = samp_pop,
                sv_fst_docker = sv_fst_docker,
                runtime_attr_override = runtime_attr_fst
            }
        }

    File Fst_sites = select_first([CalcuFstSv.Fst_sv_sites, CalcuFstSnv.Fst_snv_sites])
    File Fst_pop = select_first([CalcuFstSv.Fst_sv_pop, CalcuFstSnv.Fst_snv_pop])

    output {
        File output_fst_sites = Fst_sites
        File output_fst_pop = Fst_pop
        }
    }


task CalcuFstSv{
    input{
        File vcf
        File vcf_idx
        File samp_pop
        String sv_fst_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File Fst_sv_sites = "~{filebase}.Fst.sites"
        File Fst_sv_pop = "~{filebase}.Fst.pop"
    }

    String filebase = basename(vcf,".vcf.gz")

    command <<<
        set -Eeuo pipefail

        python /src/Calcu_Fst_table.SV.py \
            -v ~{vcf} \
            -s ~{samp_pop} \
            -o ~{filebase}.Fst.sites \
            -p ~{filebase}.Fst.pop
   >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_fst_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }    
}

task CalcuFstSnv{
    input{
        File vcf
        File vcf_idx
        File samp_pop
        String sv_fst_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File Fst_snv_sites = "~{filebase}.Fst.sites"
        File Fst_snv_pop = "~{filebase}.Fst.pop"
    }

    String filebase = basename(vcf,".vcf.gz")

    command <<<
        set -Eeuo pipefail

        python /src/Calcu_Fst_table.SNV_Indels.py \
            -v ~{vcf} \
            -s ~{samp_pop} \
            -o ~{filebase}.Fst.sites \
            -p ~{filebase}.Fst.pop
   >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_fst_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }    
}
