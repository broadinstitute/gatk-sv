##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

## Copyright Broad Institute, 2020
## 
## This WDL pipeline implements Duphold 
##
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

version 1.0

import "Structs.wdl"

workflow ApplyBoostFilterPerSample{
    input{
        File site_feature_folder
        File model_folder
        Array[File] input_folders
        Array[String] prefixes

        String sv_benchmark_docker
        RuntimeAttr? runtime_attr_boost_filter
    }

    scatter (i in range(length(input_folders))){
        call ApplyBoostFilter as apply_boost_filter{
            input:
                site_feature_folder = site_feature_folder,
                model_folder = model_folder,
                input_folder = input_folders[i],
                prefix = prefixes[i],
                sv_benchmark_docker = sv_benchmark_docker,
                runtime_attr_override = runtime_attr_boost_filter
        }
    }

    output{
        Array[File] filtered_results = apply_boost_filter.boost_results
    }
}



task ApplyBoostFilter{
    input{
        File site_feature_folder
        File model_folder
        File input_folder
        String prefix
        String sv_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 7, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File boost_results = "~{prefix}.tar.gz"
    }

    command <<<
        set -Eeuo pipefail

        mkdir site_feature_folder/
        mkdir model_folder/
        mkdir input/
        mkdir "~{prefix}/"

        tar zxvf ~{site_feature_folder} -C site_feature_folder/
        tar zxvf ~{model_folder} -C model_folder/
        tar zxvf ~{input_folder} -C input/
        
        Rscript /src/Boost_Step2.Apply_Boost_Model.R \
            --site_feature_path site_feature_folder/per_chr_anno/ \
            --model_path model_folder/lgb_models/ \
            --input_path input/ \
            --output_path ~{prefix}/ 
        
        tar zcvf ~{prefix}.tar.gz ~{prefix}/ 
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}



