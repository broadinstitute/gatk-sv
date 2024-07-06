version 1.0

import "Structs.wdl"
import "SVvsConservative.wdl" as SVvsConservative

workflow CalcuNcasStat {
    input{
        Array[Int] permutation_list
        Array[File] ncas_rdata_list
        File src_tar
        String prefix
        String sv_base_mini_docker
        File? Filter_SVID
        RuntimeAttr? runtime_attr_calcu_ncas_del
        RuntimeAttr? runtime_attr_calcu_ncas_dup
        RuntimeAttr? runtime_attr_calcu_ncas_inv    
        RuntimeAttr? runtime_attr_calcu_ncas_cpx    
        RuntimeAttr? runtime_attr_calcu_ncas_ins
        RuntimeAttr? runtime_attr_calcu_ncas_alu
        RuntimeAttr? runtime_attr_calcu_ncas_line1
        RuntimeAttr? runtime_attr_calcu_ncas_sva

    }


    scatter(i in range(length(permutation_list))){

        Boolean run_SVID_filter = defined(Filter_SVID)
        if(run_SVID_filter){
            call FilterSvSites{
                Filter_SVID = Filter_SVID,
                src_tar = src_tar,
                ncas_rdata = ncas_rdata_list[i],
                sv_base_mini_docker = sv_base_mini_docker

            }
        }

        ncas_rdata = select_first([FilterSvSites.filtered_rdata, ncas_rdata_list[i]])

        
        call CalcuNcasStat as calcu_ncas_del{
            input:
                permu = ncas_rdata,
                prefix = prefix,
                svtype = 'DEL',
                src_tar = src_tar,
                ncas_rdata = ncas_rdata_list[i],
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_calcu_ncas_del
        }

        call CalcuNcasStat as calcu_ncas_dup{
            input:
                permu = ncas_rdata,
                svtype = 'DUP',
                prefix = prefix,
                src_tar = src_tar,
                ncas_rdata = ncas_rdata_list[i],
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_calcu_ncas_dup
        }

        call CalcuNcasStat as calcu_ncas_inv{
            input:
                permu = ncas_rdata,
                svtype = 'INV',
                prefix = prefix,
                src_tar = src_tar,
                ncas_rdata = ncas_rdata_list[i],
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_calcu_ncas_inv
        }

        call CalcuNcasStat as calcu_ncas_cpx{
            input:
                permu = ncas_rdata,
                svtype = 'CPX',
                prefix = prefix,
                src_tar = src_tar,
                ncas_rdata = ncas_rdata_list[i],
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_calcu_ncas_cpx
        }

        call CalcuNcasStat as calcu_ncas_ins{
            input:
                permu = ncas_rdata,
                svtype = 'INS',
                prefix = prefix,
                src_tar = src_tar,
                ncas_rdata = ncas_rdata_list[i],
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_calcu_ncas_ins
        }

        call CalcuNcasStat as calcu_ncas_alu{
            input:
                permu = ncas_rdata,
                svtype = 'INS:ME:ALU',
                prefix = prefix,
                src_tar = src_tar,
                ncas_rdata = ncas_rdata_list[i],
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_calcu_ncas_alu
        }

        call CalcuNcasStat as calcu_ncas_line1{
            input:
                permu = ncas_rdata,
                svtype = 'INS:ME:LINE1',
                prefix = prefix,
                src_tar = src_tar,
                ncas_rdata = ncas_rdata_list[i],
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_calcu_ncas_line1
        }

        call CalcuNcasStat as calcu_ncas_sva{
            input:
                permu = ncas_rdata,
                svtype = 'INS:ME:SVA',
                prefix = prefix,
                src_tar = src_tar,
                ncas_rdata = ncas_rdata_list[i],
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_calcu_ncas_sva
        }

        call IntegrateNcasStat{
            input:
                permu = ncas_rdata,
                svtype = 'INS:ME:SVA',
                prefix = prefix,

                ncas_del = calcu_ncas_del.ncas_stat_unit,
                ncas_dup = calcu_ncas_dup.ncas_stat_unit,
                ncas_inv = calcu_ncas_inv.ncas_stat_unit,
                ncas_cpx = calcu_ncas_cpx.ncas_stat_unit,
                ncas_ins = calcu_ncas_ins.ncas_stat_unit,
                ncas_alu = calcu_ncas_alu.ncas_stat_unit,
                ncas_line1 = calcu_ncas_line1.ncas_stat_unit,
                ncas_sva = calcu_ncas_sva.ncas_stat_unit,

                src_tar = src_tar,
                sv_base_mini_docker = sv_base_mini_docker

        }

    }

    output{
        Array[File] ncas_stat = IntegrateNcasStat.ncas_stat
    }
}


task FilterSvSites{
    input{
        File? Filter_SVID
        File src_tar
        File ncas_rdata
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 25, 
        disk_gb: 40,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String filebase = basename(ncas_rdata,'.rData')

    output{
        File filtered_rdata = "~{filebase}.filtered.rData"
    }

    command <<<
        set -Eeuo pipefail

        gsutil cp ~{src_tar} ./
        tar zxvf src.tar.gz 

        Rscript src/filter_sv_sites.R -r ~{ncas_rdata} -f ~{Filter_SVID} -o ~{filebase}.filtered.rData
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



task CalcuNcasStat{
    input{
        String permu
        String prefix
        String svtype
        File src_tar
        File ncas_rdata
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 25, 
        disk_gb: 40,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File ncas_stat_unit = "~{prefix}.~{svtype}.noncoding.permu_~{permu}.stat"
    }

    command <<<
        set -Eeuo pipefail

        gsutil cp ~{src_tar} ./
        tar zxvf src.tar.gz 

        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t ~{svtype} -g noncoding -p ~{prefix}
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


task IntegrateNcasStat{
    input{
        String permu
        String prefix
        String svtype

        File ncas_del
        File ncas_dup
        File ncas_inv
        File ncas_cpx
        File ncas_ins
        File ncas_alu
        File ncas_line1
        File ncas_sva

        File src_tar
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

    output{
        File ncas_stat = "ncas_stat.permu_~{permu}.tar.gz"
    }

    command <<<
        set -Eeuo pipefail

        gsutil cp ~{src_tar} ./
        tar zxvf src.tar.gz 


        Rscript ./src/integrate_cwas_stat_across_svtype.R \
            --del   ~{ncas_del} \
            --dup   ~{ncas_dup} \
            --inv   ~{ncas_inv} \
            --cpx   ~{ncas_cpx} \
            --ins   ~{ncas_ins} \
            --alu   ~{ncas_alu} \
            --line1 ~{ncas_line1} \
            --sva   ~{ncas_sva} \
            --output "~{prefix}.ALL.noncoding.permu_~{permu}.stat"

        mkdir ncas_stat.permu_~{permu}
        mv *.stat ncas_stat.permu_~{permu}/

        gsutil cp ~{ncas_del} ./ncas_stat.permu_~{permu}/
        gsutil cp ~{ncas_dup} ./ncas_stat.permu_~{permu}/
        gsutil cp ~{ncas_inv} ./ncas_stat.permu_~{permu}/
        gsutil cp ~{ncas_cpx} ./ncas_stat.permu_~{permu}/
        gsutil cp ~{ncas_ins} ./ncas_stat.permu_~{permu}/
        gsutil cp ~{ncas_alu} ./ncas_stat.permu_~{permu}/
        gsutil cp ~{ncas_line1} ./ncas_stat.permu_~{permu}/
        gsutil cp ~{ncas_sva} ./ncas_stat.permu_~{permu}/

        tar czvf ncas_stat.permu_~{permu}.tar.gz ncas_stat.permu_~{permu}/
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





