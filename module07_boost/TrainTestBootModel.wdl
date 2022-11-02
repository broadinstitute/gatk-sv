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
import "TasksBenchmark.wdl" as mini_tasks
import "SplitPerSiteVCF.wdl" as split_per_site_vcf
import "SplitPerSampleGTGQandBEDPerSampleList.wdl" as split_per_sample_gtgq_and_bed
import "AnnotateILFeaturesBed.wdl" as annotate_il_features
import "AnnotateTrainingFeatures.wdl" as annotate_training_features

workflow TrainTestBootModel{
    input{
        Array[File] training_per_sample
        Array[File] testing_per_sample
        File site_anno

        String prefix
        String sv_types
        String size_ranges
        String af_ranges

        String sv_base_mini_docker
        String sv_benchmark_docker

        RuntimeAttr? runtime_attr_override_organize_training_data
        RuntimeAttr? runtime_attr_override_train_boost_model
        RuntimeAttr? runtime_attr_override_apply_boost_model
        RuntimeAttr? runtime_attr_override_generate_tarball
    }

    call OrganizeTrainingData{
        input:
            training_per_sample = training_per_sample,
            prefix = prefix,
            sv_types = sv_types,
            size_ranges = size_ranges,
            af_ranges = af_ranges,
            sv_benchmark_docker = sv_benchmark_docker,
            runtime_attr_override = runtime_attr_override_organize_training_data
    }

    scatter(train_data in OrganizeTrainingData.training_data){
        call TrainBoostModel{
            input:
                training_data = train_data,
                sv_benchmark_docker = sv_benchmark_docker,
                runtime_attr_override = runtime_attr_override_train_boost_model
        }
    }

    call mini_tasks.FilesToTarredFolder as FilesToTarredFolder{
        input:
            in_files = testing_per_sample,
            folder_name = "merged",
            tarball_prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_override_generate_tarball
    }

    call TestBoostModel{
            input:
                trained_models = TrainBoostModel.trained_model,
                test_data = FilesToTarredFolder.tarball,
                site_anno = site_anno,
                sv_benchmark_docker = sv_benchmark_docker,
                runtime_attr_override = runtime_attr_override_apply_boost_model
        }

    output{
        File bs_filtered_tarballs = TestBoostModel.bs_filtered_tar
    }

}



task OrganizeTrainingData{
    input{
        Array[File] training_per_sample
        String prefix
        String sv_types
        String size_ranges
        String af_ranges
        String sv_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 10,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    command<<<
        set -eux

        # note head -n1 stops reading early and sends SIGPIPE to zcat,
        # so setting pipefail here would result in early termination
        
        zcat ~{training_per_sample[0]} | head -n1 > integrated_training.tsv

        # no more early stopping
        set -o pipefail

        while read SPLIT; do
          zcat $SPLIT | tail -n+2
        done < ~{write_lines(training_per_sample)} \
          >> integrated_training.tsv

        python /src/split_ref_support.py \
            integrated_training.tsv \
            ~{prefix} \
            -t ~{sv_types} \
            -s ~{size_ranges} \
            -a ~{af_ranges}

        >>>

    output{
        File integrated_training = "integrated_training.tsv"
        Array[File] training_data = glob("~{prefix}.*.training.gz")
    }

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

task TrainBoostModel{
    input{
        File training_data
        String sv_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 20, 
        disk_gb: 25,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(training_data, '.training.gz')
    command<<<
        set -Eeuo pipefail
        Rscript /src/TrainRandomForestModel.R \
            --train_data ~{training_data} \
            --output ~{prefix}.RDS

        >>>

    output{
        File trained_model = "~{prefix}.RDS"
    }

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

task TestBoostModel{
    input{
        Array[File] trained_models
        File test_data
        File site_anno
        String sv_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 10,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(test_data, '.tar.gz')
    command<<<
        set -Eeuo pipefail

        tar zxvf ~{test_data}
        ls merged/*.gz > test_data_list.tsv

        Rscript /src/apply_boost_models.R \
            --input test_data_list.tsv \
            --site ~{site_anno} \
            --Boost_model_del_sz1_af1 ~{trained_models[0]} \
            --Boost_model_del_sz1_af2 ~{trained_models[1]} \
            --Boost_model_del_sz1_af3 ~{trained_models[2]} \
            --Boost_model_del_sz1_af4 ~{trained_models[3]} \
            --Boost_model_del_sz2_af1 ~{trained_models[4]} \
            --Boost_model_del_sz2_af2 ~{trained_models[5]} \
            --Boost_model_del_sz2_af3 ~{trained_models[6]} \
            --Boost_model_del_sz2_af4 ~{trained_models[7]} \
            --Boost_model_del_sz3_af1 ~{trained_models[8]} \
            --Boost_model_del_sz3_af2 ~{trained_models[9]} \
            --Boost_model_del_sz3_af3 ~{trained_models[10]} \
            --Boost_model_del_sz3_af4 ~{trained_models[11]} \
            --Boost_model_del_sz4_af1 ~{trained_models[12]} \
            --Boost_model_del_sz4_af2 ~{trained_models[13]} \
            --Boost_model_del_sz4_af3 ~{trained_models[14]} \
            --Boost_model_del_sz4_af4 ~{trained_models[15]} \
            --Boost_model_dup_sz1_af1 ~{trained_models[16]} \
            --Boost_model_dup_sz1_af2 ~{trained_models[17]} \
            --Boost_model_dup_sz1_af3 ~{trained_models[18]} \
            --Boost_model_dup_sz1_af4 ~{trained_models[19]} \
            --Boost_model_dup_sz2_af1 ~{trained_models[20]} \
            --Boost_model_dup_sz2_af2 ~{trained_models[21]} \
            --Boost_model_dup_sz2_af3 ~{trained_models[22]} \
            --Boost_model_dup_sz2_af4 ~{trained_models[23]} \
            --Boost_model_dup_sz3_af1 ~{trained_models[24]} \
            --Boost_model_dup_sz3_af2 ~{trained_models[25]} \
            --Boost_model_dup_sz3_af3 ~{trained_models[26]} \
            --Boost_model_dup_sz3_af4 ~{trained_models[27]} \
            --Boost_model_dup_sz4_af1 ~{trained_models[28]} \
            --Boost_model_dup_sz4_af2 ~{trained_models[29]} \
            --Boost_model_dup_sz4_af3 ~{trained_models[30]} \
            --Boost_model_dup_sz4_af4 ~{trained_models[31]} \
            --Boost_model_ins_sz1 ~{trained_models[32]} \
            --Boost_model_ins_sz2 ~{trained_models[33]} \
            --Boost_model_ins_sz3 ~{trained_models[34]} \
            --Boost_model_ins_sz4 ~{trained_models[35]} \
            --Boost_model_ins_alu ~{trained_models[36]} \
            --Boost_model_ins_l1 ~{trained_models[37]} \
            --Boost_model_ins_sva ~{trained_models[38]} \
            --Boost_model_ins_me ~{trained_models[39]} \
            --Boost_model_inv ~{trained_models[40]}

        mkdir bs_filtered/
        mv merged/*.bs bs_filtered/
        tar -czvf ~{prefix}.bs_filtered.tar.gz bs_filtered/
        >>>

    output{
        File bs_filtered_tar = "~{prefix}.bs_filtered.tar.gz"
    }

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
