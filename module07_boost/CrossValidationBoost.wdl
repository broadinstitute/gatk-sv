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
import "TrainTestBootModel.wdl" as train_test_boost_model

workflow CrossValidationBoost{
    input{

        Array[File] training_anno
        Int cross_validation_fold
        File site_anno

        String sv_types
        String size_ranges
        String af_ranges

        String sv_benchmark_docker
        String sv_base_mini_docker

        RuntimeAttr? runtime_attr_override_prepare_cross_validation
        RuntimeAttr? runtime_attr_override_organize_training_data
        RuntimeAttr? runtime_attr_override_generate_tarball
        RuntimeAttr? runtime_attr_override_apply_boost_model
        RuntimeAttr? runtime_attr_override_train_boost_model
        RuntimeAttr? runtime_attr_override_seek_for_optimal_cffs
    }


    call SplitTrainingIntoCrossValidation{
        input:
            training_per_sample = training_anno,
            cross_validation_fold  = cross_validation_fold,
            sv_benchmark_docker = sv_benchmark_docker,
            runtime_attr_override = runtime_attr_override_prepare_cross_validation

    }

    scatter(i in range(length(SplitTrainingIntoCrossValidation.CV_train))){
        Array[File] CV_train_files = read_lines(SplitTrainingIntoCrossValidation.CV_train[i])
        Array[File] CV_test_files = read_lines(SplitTrainingIntoCrossValidation.CV_test[i])
        String prefix_unit = basename(SplitTrainingIntoCrossValidation.CV_train[i],'.train')

        call train_test_boost_model.TrainTestBootModel as TrainTestBootModel{
            input:
                training_per_sample = CV_train_files,
                testing_per_sample = CV_test_files,
                site_anno = site_anno,
                prefix = prefix_unit,
                sv_types = sv_types,
                size_ranges = size_ranges,
                af_ranges = af_ranges,
                sv_benchmark_docker = sv_benchmark_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override_organize_training_data = runtime_attr_override_organize_training_data,
                runtime_attr_override_train_boost_model = runtime_attr_override_train_boost_model,
                runtime_attr_override_apply_boost_model = runtime_attr_override_apply_boost_model
        }
    }


    call SeekForOptimalCutoffs{
        input:
            bs_tarball = TrainTestBootModel.bs_filtered_tarballs,
            sv_benchmark_docker = sv_benchmark_docker,
            runtime_attr_override = runtime_attr_override_seek_for_optimal_cffs
    }

    output{
        File cff_performance_table = SeekForOptimalCutoffs.cff_table
    }
}



task SplitTrainingIntoCrossValidation{
    input{
        Array[File] training_per_sample
        Int cross_validation_fold
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

    File training_list = write_lines(training_per_sample)

    command<<<
        set -Eeuo pipefail
        python3 <<CODE
        fin=open("~{training_list}")
        all_list = []
        for line in fin:
            pin=line.strip().split()
            all_list+=pin
        fin.close()

        train_len = len(all_list)
        cv_len = int(train_len/~{cross_validation_fold})
        for i in range(~{cross_validation_fold}):
            test_list = all_list[(i*cv_len):(i*cv_len+cv_len+1)]
            train_list = [j for j in all_list if not j in test_list]
            fo=open('CV_cross_'+str(i)+'.train','w')
            for j in train_list:
                print(j, file=fo)
            fo.close()
            fo=open('CV_cross_'+str(i)+'.test','w')
            for j in test_list:
                print(j, file=fo)
            fo.close()
        CODE
    >>>

    output{
        Array[File] CV_train = glob("CV_cross_*.train")
        Array[File] CV_test = glob("CV_cross_*.test")
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

task SeekForOptimalCutoffs{
    input{
        Array[File] bs_tarball
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
        set -Eeuo pipefail
        mkdir tmp/
        while read SPLIT; do
            tar zxvf $SPLIT
            mv bs_filtered/* tmp/
        done < ~{write_lines(bs_tarball)}

        Rscript /src/SeekForOptimalCff.R \
        --input_folder tmp/ \
        --output cff_table.tsv
    >>>

    output{
        File cff_table = 'cff_table.tsv'
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





