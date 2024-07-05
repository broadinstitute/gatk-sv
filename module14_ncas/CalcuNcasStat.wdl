version 1.0

import "Structs.wdl"
import "SVvsConservative.wdl" as SVvsConservative

workflow CalcuNcasStat {
    input{
        Array[Int] permutation_list
        Array[File] ncas_rdata_list
        File src_tar
        File ref_tar
        String prefix
        File contig_file
        String sv_base_mini_docker
    }


    call CalcuNcasStat{
        input:
            permu = i,
            prefix = prefix,
            src_tar = src_tar,
            ncas_rdata = GenerateNcasMetrics.ncas_rdata,
            sv_base_mini_docker = sv_base_mini_docker
    }
    }
    output{
        Array[File] ncas_data_list = GenerateNcasMetrics.ncas_rdata
        Array[File] ncas_stat = CalcuNcasStat.ncas_stat
    }
}


task CalcuNcasStat{
    input{
        String permu
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
        File ncas_stat = "ncas_stat.permu_~{permu}.tar.gz"
    }

    command <<<
        set -Eeuo pipefail

        gsutil cp ~{src_tar} ./
        tar zxvf src.tar.gz 

        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t DEL -g noncoding -p ~{prefix}
        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t DUP -g noncoding -p ~{prefix}
        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t INV -g noncoding -p ~{prefix}
        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t CPX -g noncoding -p ~{prefix}
        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t INS -g noncoding -p ~{prefix}
        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t INS:ME:ALU -g noncoding -p ~{prefix}
        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t INS:ME:SVA -g noncoding -p ~{prefix}
        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t INS:ME:LINE1 -g noncoding -p ~{prefix}
        Rscript ./src/integrate_cwas_stat_across_svtype.R \
            --del   "~{prefix}.DEL.noncoding.permu_~{permu}.stat" \
            --dup   "~{prefix}.DUP.noncoding.permu_~{permu}.stat" \
            --inv   "~{prefix}.INV.noncoding.permu_~{permu}.stat" \
            --cpx   "~{prefix}.CPX.noncoding.permu_~{permu}.stat" \
            --ins   "~{prefix}.INS.noncoding.permu_~{permu}.stat" \
            --alu   "~{prefix}.INS:ME:ALU.noncoding.permu_~{permu}.stat" \
            --line1 "~{prefix}.INS:ME:LINE1.noncoding.permu_~{permu}.stat" \
            --sva   "~{prefix}.INS:ME:SVA.noncoding.permu_~{permu}.stat" \
            --output "~{prefix}.ALL.noncoding.permu_~{permu}.stat"

        mkdir ncas_stat.permu_~{permu}
        mv *.stat ncas_stat.permu_~{permu}/
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


