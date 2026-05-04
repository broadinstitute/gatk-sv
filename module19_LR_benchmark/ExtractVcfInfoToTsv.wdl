version 1.0

import "Structs.wdl"

workflow ExtractVcfInfoToTsv {
    input {
        Array[File] vcf_list
        File extract_script          # extract_vcf_info_to_tsv.py
        String output_suffix = ".info.tsv"
        String python_docker
        RuntimeAttr? runtime_attr_override
    }

    scatter (vcf in vcf_list) {
        call ExtractInfoTask {
            input:
                vcf            = vcf,
                extract_script = extract_script,
                output_suffix  = output_suffix,
                docker         = python_docker,
                runtime_attr_override = runtime_attr_override
        }
    }

    output {
        Array[File] tsv_list = ExtractInfoTask.out_tsv
    }
}

task ExtractInfoTask {
    input {
        File   vcf
        File   extract_script
        String output_suffix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size   = size(vcf, "GB")
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
        mem_gb:            4,
        disk_gb:           ceil(base_disk_gb + input_size * 2.0),
        cpu_cores:         1,
        preemptible_tries: 2,
        max_retries:       1,
        boot_disk_gb:      10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    String prefix = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail

        python3 ~{extract_script} \
            --input-vcf  ~{vcf} \
            --output-txt ~{prefix}~{output_suffix}
    >>>

    output {
        File out_tsv = "~{prefix}~{output_suffix}"
    }

    runtime {
        docker:         docker
        memory:         select_first([runtime_override.mem_gb,            runtime_default.mem_gb])           + " GB"
        disks:          "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb])   + " HDD"
        cpu:            select_first([runtime_override.cpu_cores,          runtime_default.cpu_cores])
        preemptible:    select_first([runtime_override.preemptible_tries,  runtime_default.preemptible_tries])
        maxRetries:     select_first([runtime_override.max_retries,        runtime_default.max_retries])
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb,       runtime_default.boot_disk_gb])
    }
}
