version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks

# Reheader a list of vcfs with the header from another vcf

workflow HarmonizeHeaders {
    input {
        File header_vcf     # Vcf containing desired header
        Array[File] vcfs    # Vcfs to replace headers of
        String prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_override_reheader
        RuntimeAttr? runtime_override_pull_header
    }

    call PullHeader {
        input:
            vcf=header_vcf,
            prefix=prefix,
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_override_pull_header
    }

    scatter (i in range(length(vcfs))) {
        call MiniTasks.ReheaderVcf {
            input:
                vcf=vcfs[i],
                vcf_index=vcfs[i] + ".tbi",
                header=PullHeader.out,
                prefix="~{prefix}.reheadered",
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_override_reheader
        }
    }

    output {
        Array[File] out = ReheaderVcf.out
        Array[File] out_index = ReheaderVcf.out_index
    }
}

task PullHeader {
    input {
        File vcf
        String prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
                                      mem_gb: 2.0,
                                      disk_gb: ceil(10.0 + size(vcf, "GiB") ),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        bcftools view --header-only ~{vcf} > ~{prefix}.header
    >>>

    output {
        File out = "~{prefix}.header"
    }
}