version 1.0

# Author: Xuefang Zhao <XZHAO12@mgh.harvard.edu>

import "Structs.wdl"

# Workflow to annotate vcf file with genomic context
workflow ShardSVsByGenomicContext {
    input {
        File vcf
        File vcf_index
        File SVID_GC
        String Genomic_Context
        String svtype

        String sv_base_mini_docker
        String sv_benchmark_docker

        # overrides for MiniTasks
        RuntimeAttr? runtime_attr_override_extract_svid
       }

    call Extract_SVID_by_GC {
        input:
            SVID_GC = SVID_GC, 
            vcf = vcf, 
            vcf_index = vcf_index,
            Genomic_Context = Genomic_Context,
            svtype = svtype,
            sv_benchmark_docker = sv_benchmark_docker,
            runtime_attr_override = runtime_attr_override_extract_svid
    }

    output{
        File sharded_vcf = Extract_SVID_by_GC.out
        File sharded_vcf_idx = Extract_SVID_by_GC.out_idx
    }
}

task Extract_SVID_by_GC{
    input {
        File SVID_GC
        File vcf
        File vcf_index
        String Genomic_Context
        String svtype
        String sv_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    Float vcf_size = size(vcf, "GiB")
    Int vm_disk_size = ceil(vcf_size * 2)

    RuntimeAttr runtime_default = object {
        mem_gb: 1,
        disk_gb: vm_disk_size,
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_benchmark_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(vcf, ".vcf.gz")

    command <<<
        set -euxo pipefail

        python /src/extract_SVs_by_svtype_genomiccontect.py \
        ~{vcf} \
        ~{prefix}.~{svtype}.~{Genomic_Context}.vcf.gz \
        ~{SVID_GC} \
        ~{svtype} \
        ~{Genomic_Context} 

        tabix -p vcf ~{prefix}.~{svtype}.~{Genomic_Context}.vcf.gz 

    >>>

    output {
        File out = "~{prefix}.~{svtype}.~{Genomic_Context}.vcf.gz"
        File out_idx = "~{prefix}.~{svtype}.~{Genomic_Context}.vcf.gz.tbi"
    } 
}


