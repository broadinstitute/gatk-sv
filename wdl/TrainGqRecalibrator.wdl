version 1.0

import "Utils.wdl" as Utils

workflow TrainGqRecalibrator {
    input {
        File train_vcf
        File train_vcf_index
        File truth_json
        Array[File] genome_tracks
        File? gq_recalibrator_model_file
        Array[String] train_args = []
        String gatk_docker
        String samtools_cloud_docker
        Int? preemptible_tries
    }

    call Utils.GetVcfSize {
        input:
            vcf=train_vcf,
            vcf_index=train_vcf_index,
            samtools_cloud_docker=samtools_cloud_docker
    }

    call TrainGqRecalibratorTask {
        input:
            train_vcf=train_vcf,
            train_vcf_index=train_vcf_index,
            truth_file=truth_json,
            genome_tracks=genome_tracks,
            gq_recalibrator_model_file=gq_recalibrator_model_file,
            train_args=train_args,
            gatk_docker=gatk_docker,
            num_entries=GetVcfSize.num_entries,
            preemptible_tries=preemptible_tries
    }

    output {
        File output_gq_recalibrator_model_file = TrainGqRecalibratorTask.output_gq_recalibrator_model_file
    }
}

task TrainGqRecalibratorTask {
    input {
        File train_vcf
        File train_vcf_index
        File? ped_file
        File truth_file
        Array[File] genome_tracks
        File? gq_recalibrator_model_file # can be passed to do extra rounds of training on existing model
        Array[String] train_args = []
        String gatk_docker
        Float? num_entries
        Float mem_scale_vcf_size = 25.2
        Float mem_scale_num_entries = "3.7e-7"
        Float mem_gb_overhead = 1.5
        Int preemptible_tries = 3
    }

    Int disk_gb = round(1000 + size([train_vcf, train_vcf_index, ped_file, truth_file], "GiB") +
                        size(genome_tracks, "GiB") + size(gq_recalibrator_model_file, "GiB"))
    Float mem_gb_java = if defined(num_entries)
        then 3.0 + mem_scale_num_entries * select_first([num_entries])
        else 3.0 + mem_scale_vcf_size * size(train_vcf, "GiB")
    Float max_mem_gb = 624  # current max n1-highmem, which is what cromwell is willing to allocate
    Float mem_gb = if mem_gb_java + mem_gb_overhead < max_mem_gb then mem_gb_java + mem_gb_overhead else max_mem_gb
    String model_file_name = if defined(gq_recalibrator_model_file)
        then basename(select_first([gq_recalibrator_model_file]))
        else "gq_recalibrator.model"

    runtime {
        docker: gatk_docker
        cpu: 1
        preemptible: preemptible_tries
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -euo pipefail
        ln -s ~{train_vcf} .
        ln -s ~{train_vcf_index} .
        if ~{defined(gq_recalibrator_model_file)}; then
            # if defined, copy exiting model file to working directory
            GQ_RECALIBRATOR_MODEL_FILE="~{if defined(gq_recalibrator_model_file) then select_first([gq_recalibrator_model_file]) else ""}"
            cp $GQ_RECALIBRATOR_MODEL_FILE ./~{model_file_name}
        fi

        mem_kb_java_actual=$(grep -m1 MemTotal /proc/meminfo \
                             | awk '{printf("%.0f\n", $2 - ~{mem_gb_overhead} * 1048576)}')

        NUM_PHYSICAL_CPUS=$(grep ^"core id" /proc/cpuinfo | sort -u | wc -l)

        gatk --java-options "-Xmx${mem_kb_java_actual}K" XGBoostMinGqVariantFilter \
          --mode "Train" \
          --variant ./$(basename ~{train_vcf}) \
          ~{"--pedigree " + ped_file} \
          --truth-file ~{truth_file} \
          --genome-track ~{sep=" --genome-track " genome_tracks} \
          --model-file ~{model_file_name} \
          --n-threads-xgboost $NUM_PHYSICAL_CPUS \
          ~{sep=' ' train_args}
    >>>

    output {
        File output_gq_recalibrator_model_file = model_file_name
    }
}
