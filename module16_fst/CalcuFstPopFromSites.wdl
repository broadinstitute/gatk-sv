version 1.0

import "Structs.wdl"
import "CalculateFstGnomad.wdl" as calculate_fst

workflow CalcuFstPopFromSites {
    input{
        File Fst_sites
        String sv_fst_docker
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_fst_pop_from_sites
    }

    call Unzip{
        input:
            file = Fst_sites,
            sv_base_mini_docker  = sv_base_mini_docker
    }

    call CalcuFstPop{
        input:
            Fst_sites = Unzip.output_file,
            sv_fst_docker = sv_fst_docker,
            runtime_attr_override = runtime_attr_fst_pop_from_sites
    }

    output{
        File Fst_pop =  CalcuFstPop.fst_pop
    }
}


task Unzip{
    input{
        File file
        String sv_base_mini_docker  
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

    String prefix = basename(file, '.gz')

    output{
        File output_file = "~{prefix} "
    }


    command <<<
        set -Eeuo pipefail

        zcat ~{file} > ~{prefix} 
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


task CalcuFstPop{
    input{
        File Fst_sites
        String sv_fst_docker
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

    output{
        File fst_pop = "~{filebase}.pop"
    }

    String filebase = basename(Fst_sites,".sites")

    command <<<
        set -Eeuo pipefail

        head -1 ~{Fst_sites}  > temp.sites

        set -o pipefail

         grep -v "_num" ~{Fst_sites}  >> temp.sites

        python /src/Calcu_Fst_pop_from_sites.py -i temp.sites -p ~{filebase}.pop
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

