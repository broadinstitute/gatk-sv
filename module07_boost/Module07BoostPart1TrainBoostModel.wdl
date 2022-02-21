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

workflow TrainBoostModel{
    input{
        File site_feature_folder
        File sample_feature_folder
        File LD_feature_folder
        File lg_cnv_folder
        File sv_colors

        String prefix
        String sv_benchmark_docker
        RuntimeAttr? runtime_attr_boost_model
    }

    call GenerateBoostModel as generate_boost_model{
        input:
            site_feature_folder = site_feature_folder,
            sample_feature_folder = sample_feature_folder,
            LD_feature_folder = LD_feature_folder,
            lg_cnv_folder = lg_cnv_folder,
            sv_colors = sv_colors,
            prefix = prefix,
            sv_benchmark_docker = sv_benchmark_docker,
            runtime_attr_override = runtime_attr_boost_model
        }

    output{
        File filtered_results = generate_boost_model.boost_models
    }
}



task GenerateBoostModel{
    input{
        File site_feature_folder
        File sample_feature_folder
        File LD_feature_folder
        File lg_cnv_folder
        File sv_colors
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
        File boost_models = "~{prefix}.tar.gz"
    }

    command <<<
        set -Eeuo pipefail


        mkdir site_feature_folder/
        mkdir sample_feature_folder/
        mkdir LD_feature_folder/
        mkdir lg_cnv_folder/
        mkdir "~{prefix}/"
        mkdir "~{prefix}/lgb_models/"
        mkdir "~{prefix}/qc_plots/"

        tar zxvf ~{site_feature_folder} -C site_feature_folder/
        tar zxvf ~{sample_feature_folder} -C sample_feature_folder/
        tar zxvf ~{LD_feature_folder} -C LD_feature_folder/
        tar zxvf ~{lg_cnv_folder}     -C lg_cnv_folder/
        
        Rscript /src/Boost_Step1.Train_Boost_Model.R \
            --site_feature_path site_feature_folder/per_chr_anno/ \
            --sample_feature_path sample_feature_folder/per_sample_anno/ \
            --LD_feature_path LD_feature_folder/LD_training/ \
            --lg_cnv_path lg_cnv_folder/array_training/ \
            --model_path "~{prefix}/lgb_models/" \
            --qc_plots "~{prefix}/qc_plots/" \
            --sv_colors ~{sv_colors}

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



