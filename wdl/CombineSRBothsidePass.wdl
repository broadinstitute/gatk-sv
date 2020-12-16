version 1.0

import "Structs.wdl"

workflow CombineSRBothsidePass {
    input {
        Array[File] pesr_vcfs
        Array[File] raw_sr_bothside_pass_files
        String prefix

        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_get_non_ref_vids
        RuntimeAttr? runtime_attr_calculate_support_frac
    }

    scatter (i in range(length(pesr_vcfs))) {
        call GetNonRefVariantLists {
            input:
                vcf=pesr_vcfs[i],
                prefix="~{prefix}.non_ref_vids.shard_~{i}",
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_get_non_ref_vids
        }
    }

    call CalculateBothsideSupportFraction {
        input:
            non_ref_vid_lists=GetNonRefVariantLists.out,
            raw_sr_bothside_pass_files=raw_sr_bothside_pass_files,
            prefix="~{prefix}.sr_bothside_support",
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override=runtime_attr_calculate_support_frac
    }

    output {
        File out = CalculateBothsideSupportFraction.out
    }
}

task GetNonRefVariantLists {
    input {
        File vcf
        String prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf, "GB")
    RuntimeAttr runtime_default = object {
                                      mem_gb: 3.75,
                                      disk_gb: ceil(10.0 + input_size * 2.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        bcftools view -G -i 'SUM(AC)>0||SUM(FORMAT/SR_GT)>0' ~{vcf} | bcftools query -f '%ID\n' \
            > ~{prefix}.list
    >>>
    output {
        File out = "~{prefix}.list"
    }
}


task CalculateBothsideSupportFraction {
    input {
        Array[File] non_ref_vid_lists
        Array[File] raw_sr_bothside_pass_files
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(non_ref_vid_lists, "GB") + size(raw_sr_bothside_pass_files, "GB")
    RuntimeAttr runtime_default = object {
                                      mem_gb: 3.75,
                                      disk_gb: ceil(10.0 + input_size * 2.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        python /opt/sv-pipeline/04_variant_resolution/scripts/calculate_sr_bothside_support.py \
            ~{write_lines(non_ref_vid_lists)} \
            ~{write_lines(raw_sr_bothside_pass_files)} \
            > ~{prefix}.txt
    >>>
    output {
        File out = "~{prefix}.txt"
    }
}