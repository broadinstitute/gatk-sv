version 1.0

import "TrainGqRecalibrator.wdl" as TrainGqRecalibrator

workflow RecalibrateGq {
    input {
        File vcf
        File vcf_index
        Array[File] genome_tracks
        File gq_recalibrator_model_file
        Boolean standardize_vcf = true
        Array[String] standardize_vcf_args = []
        Array[String] recalibrate_gq_args = []
        String gatk_docker
        String sv_utils_docker
    }

    if(standardize_vcf) {
        call TrainGqRecalibrator.StandardizeVcfForGatk {
            input:
                vcf=vcf,
                standardize_vcf_args=standardize_vcf_args,
                sv_utils_docker=sv_utils_docker
        }
    }

    File vcf_ = select_first([StandardizeVcfForGatk.fixed_vcf, vcf])
    File vcf_index_ = select_first([StandardizeVcfForGatk.fixed_vcf_index, vcf_index])

    call RecalibrateGqTask {
        input:
            vcf=vcf_,
            vcf_index=vcf_index_,
            genome_tracks=genome_tracks,
            gq_recalibrator_model_file=gq_recalibrator_model_file,
            recalibrate_gq_args = recalibrate_gq_args,
            gatk_docker=gatk_docker
    }

    output {
        File filtered_vcf = RecalibrateGqTask.filtered_vcf
        File filtered_vcf_index = RecalibrateGqTask.filtered_vcf_index
    }
}


task RecalibrateGqTask {
    input {
        File vcf
        File vcf_index
        Array[File] genome_tracks
        File gq_recalibrator_model_file
        Array[String] recalibrate_gq_args = []
        String gatk_docker
        Float mem_gb_java = 9.0
        Float mem_gb_overhead = 1.5
    }

    Int disk_gb = round(1000 + size([vcf, vcf_index, gq_recalibrator_model_file], "GiB") + size(genome_tracks, "GiB"))
    Float mem_gb = mem_gb_java + mem_gb_overhead
    String filtered_vcf_name = sub(sub(basename(vcf), ".gz", ""), ".vcf", "_gq_recalibrated.vcf.gz")

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/tbrookin/gatk:900ac4539"
        cpu: 1
        preemptible: 3
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -euo pipefail

        ln -s ~{vcf} .
        ln -s ~{vcf_index} .

        mem_kb_java_actual=$(grep -m1 MemTotal /proc/meminfo \
                             | awk '{printf("%.0f\n", $2 - ~{mem_gb_overhead} * 1048576)}')

        gatk --java-options "-Xmx${mem_kb_java_actual}K" XGBoostMinGqVariantFilter \
          --mode "Filter" \
          --variant ~{basename(vcf)} \
          ~{if length(genome_tracks) > 0 then "--genome-track" else ""} ~{sep=" --genome-track " genome_tracks} \
          --model-file ~{gq_recalibrator_model_file} \
          --output ~{filtered_vcf_name} \
          ~{sep=' ' recalibrate_gq_args}

        # gatk indices still have problems, overwrite with bcftools index
        bcftools index --force --tbi --threads 2 ~{filtered_vcf_name}
    >>>

    output {
        File filtered_vcf = filtered_vcf_name
        File filtered_vcf_index = filtered_vcf_name + ".tbi"
    }
}
