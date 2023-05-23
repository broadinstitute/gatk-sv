version 1.0

import "Structs.wdl"

import "Duphold.wdl" as duphold
import "RdPeSrAnno.wdl" as rdpesr

workflow AnnoRdPeSr {
    input{

        File pe_matrix
        File pe_index
        File sr_matrix
        File sr_index
        File rd_matrix
        File rd_index

        File contig_list

        File bed
        File bed_le_flank
        File bed_ri_flank
        String sample
        String prefix

        String rdpesr_benchmark_docker
        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_rdpesr
    }

    call RunRdPeSrAnnotation{
        input:
            prefix = prefix,
            bed = bed,
            bed_le_flank = bed_le_flank,
            bed_ri_flank = bed_ri_flank,
            pe_matrix = pe_matrix,
            pe_index = pe_index,
            sr_matrix = sr_matrix,
            sr_index = sr_index,
            rd_matrix = rd_matrix,
            rd_index = rd_index,
            rdpesr_benchmark_docker = rdpesr_benchmark_docker,
            runtime_attr_override = runtime_attr_rdpesr
    }
    
    output{
            File PesrAnno  = RunRdPeSrAnnotation.pesr_anno
            File RdAnno    = RunRdPeSrAnnotation.cov
            File RdAnno_le = RunRdPeSrAnnotation.cov_le_flank 
            File RdAnno_ri = RunRdPeSrAnnotation.cov_ri_flank
        }
    }
 
task RunRdPeSrAnnotation{
    input{
        String prefix
        File bed
        File bed_le_flank
        File bed_ri_flank
        File pe_matrix
        File pe_index
        File sr_matrix
        File sr_index
        File rd_matrix
        File rd_index
        String rdpesr_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 15, 
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File cov = "~{filebase}.bed.Rd.gz"
        File cov_ri_flank = "~{filebase}.ri_flank.Rd.gz"
        File cov_le_flank = "~{filebase}.le_flank.Rd.gz"
        File pesr_anno = "~{filebase}.bed.PeSr.gz"
    }

    String filebase = basename(bed,".bed")

    command <<<

        set -Eeuo pipefail


        zcat ~{rd_matrix} | grep -v '@' | grep -v CONTIG |bgzip >    bincov.tsv.gz
        Rscript /src/bincov_to_normCov.R -i bincov.tsv.gz
        bgzip normCov.tsv
        tabix -b 2 -e 2 normCov.tsv.gz

        python3 /src/add_RD_to_SVs.py ~{bed} normCov.tsv.gz ~{filebase}.bed.Rd
        python3 /src/add_RD_to_SVs.py ~{bed_le_flank} normCov.tsv.gz ~{filebase}.le_flank.Rd
        python3 /src/add_RD_to_SVs.py ~{bed_ri_flank} normCov.tsv.gz ~{filebase}.ri_flank.Rd
        python3 /src/add_SR_PE_to_PB_INS.V2.py ~{bed} ~{pe_matrix} ~{sr_matrix} ~{filebase}.bed.PeSr

        bgzip ~{filebase}.bed.Rd
        bgzip ~{filebase}.ri_flank.Rd
        bgzip ~{filebase}.le_flank.Rd
        bgzip ~{filebase}.bed.PeSr

    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: rdpesr_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}


