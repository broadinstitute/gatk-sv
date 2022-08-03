version 1.0

# Author: Xuefang Zhao <XZHAO12@mgh.harvard.edu>

import "Structs.wdl"

# Workflow to annotate vcf file with genomic context
workflow TrainRandomForestModelsBySVTYPE {
    input {
        File training_data
        File inheritance
        File genomic_context

        String svtype
        Array[String] size_range_list

        String prefix

        String sv_benchmark_docker

        RuntimeAttr? runtime_attr_override_train_RF_model

    }

    scatter (size_cate in size_range_list){
        call train_random_forest_model{
            input:
                inheritance = inheritance,
                genomic_context = genomic_context,
                train = training_data,
                svtype = svtype,
                size_cate = size_cate,
                prefix = prefix,
                sv_benchmark_docker = sv_benchmark_docker,
                runtime_attr_override = runtime_attr_override_train_RF_model
        }
    }


    output{
        Array[File] pbsv_models = train_random_forest_model.model_pbsv
        Array[File] vapor_models = train_random_forest_model.model_vapor
    }
}

task train_random_forest_model{
    input{
        File inheritance
        File genomic_context
        File train
        
        String prefix
        String svtype
        String size_cate

        String sv_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }


    RuntimeAttr runtime_default = object {
        mem_gb: 10,
        disk_gb: 15,
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



    command <<<
        set -euo pipefail

        Rscript /src/TrainRandomForestModel.R \
        --inheri ~{inheritance} \
        --genomic_context ~{genomic_context} \
        --train ~{train} \
        --output_PBSV "PBSV_~{prefix}_~{svtype}_~{size_cate}.RDS" \
        --output_VaPoR "VaPoR_~{prefix}_~{svtype}_~{size_cate}.RDS" \
        --svtype ~{svtype} \
        --size_cate ~{size_cate}

    >>>

    output{
        File model_pbsv = "PBSV_~{prefix}_~{svtype}_~{size_cate}.RDS"
        File model_vapor = "VaPoR_~{prefix}_~{svtype}_~{size_cate}.RDS"
    }
}





