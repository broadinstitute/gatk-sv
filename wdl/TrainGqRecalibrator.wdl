version 1.0

import "GetTruthOverlap.wdl" as GetTruthOverlap
import "Utils.wdl" as Utils

workflow TrainGqRecalibrator {
    input {
        File train_vcf
        File train_vcf_index
        File? annotations_vcf
        File? annotations_vcf_index
        Array[String]? annotations_to_transfer
        Array[File] truth_vcfs
        Array[File] truth_vcf_indices
        Array[String]? vapor_sample_ids
        Array[File]? vapor_files
        File ped_file
        Array[File] genome_tracks
        File? optimal_overlap_cutoffs
        File? gq_recalibrator_model_file
        Array[String] standardize_vcf_args = []
        Array[String] train_args = []
        Array[String] get_truth_overlap_args = []
        Boolean standardize_vcf = true
        String sv_utils_docker
        String gatk_docker
        String samtools_cloud_docker
    }

    if(standardize_vcf) {
        call StandardizeVcfForGatk {
            input:
                vcf=train_vcf,
                standardize_vcf_args = standardize_vcf_args,
                sv_utils_docker=sv_utils_docker
        }
    }

    if(defined(annotations_vcf)) {
        call Utils.TransferVcfAnnotations {
            input:
                vcf_to_annotate=select_first([StandardizeVcfForGatk.fixed_vcf, train_vcf]),
                vcf_to_annotate_index=select_first([StandardizeVcfForGatk.fixed_vcf_index, train_vcf_index]),
                vcf_with_annotations=select_first([annotations_vcf]),
                vcf_with_annotations_index=select_first([annotations_vcf_index]),
                annotations_to_transfer=select_first([annotations_to_transfer]),
                samtools_cloud_docker=samtools_cloud_docker
        }
    }

    File train_vcf_ = select_first([TransferVcfAnnotations.annotated_vcf, StandardizeVcfForGatk.fixed_vcf, train_vcf])
    File train_vcf_index_ = select_first([TransferVcfAnnotations.annotated_vcf_index,
                                          StandardizeVcfForGatk.fixed_vcf_index, train_vcf_index])

    call GetTruthOverlap.GetTruthOverlap {
        input:
            test_vcfs=[train_vcf_],
            test_vcf_indices=[train_vcf_index_],
            truth_vcfs=truth_vcfs,
            truth_vcf_indices=truth_vcf_indices,
            vapor_sample_ids=vapor_sample_ids,
            vapor_files=vapor_files,
            ped_files=[ped_file],
            optimal_overlap_cutoffs=optimal_overlap_cutoffs,
            get_truth_overlap_args=get_truth_overlap_args,
            sv_utils_docker=sv_utils_docker,
            samtools_cloud_docker=samtools_cloud_docker
    }

    call Utils.GetVcfSize {
        input:
            vcf=train_vcf_,
            vcf_index=train_vcf_index_,
            samtools_cloud_docker=samtools_cloud_docker
    }

    call TrainGqRecalibratorTask {
        input:
            train_vcf=train_vcf_,
            train_vcf_index=train_vcf_index_,
            ped_file=ped_file,
            truth_file=GetTruthOverlap.truth_overlap_info,
            genome_tracks=genome_tracks,
            gq_recalibrator_model_file=gq_recalibrator_model_file,
            train_args=train_args,
            gatk_docker=gatk_docker,
            num_entries=GetVcfSize.num_entries
    }

    output {
        File clean_vcf=train_vcf_
        File clean_vcf_index=train_vcf_index_
        File truth_overlap_info = GetTruthOverlap.truth_overlap_info
        File output_optimal_overlap_cutoffs = GetTruthOverlap.output_optimal_overlap_cutoffs
        File output_gq_recalibrator_model_file = TrainGqRecalibratorTask.output_gq_recalibrator_model_file
    }
}


task StandardizeVcfForGatk {
    input {
        File vcf
        Array[String] standardize_vcf_args
        String sv_utils_docker
    }

    String fixed_vcf_name = sub(sub(basename(vcf), ".gz$", ""), ".vcf$", "_fixed.vcf.gz")
    String index_file_name = fixed_vcf_name + ".tbi"

    Int disk_gb = 100 + round(2 * size(vcf, "GiB"))
    Int mem_gb_overhead = 2
    Float mem_scale = 2.0
    Float mem_gb = mem_gb_overhead + mem_scale * size(vcf, "GiB")

    runtime {
        docker: sv_utils_docker
        cpu: 1
        preemptible: 3
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -euo pipefail

        sv-utils fix-vcf ~{vcf} ~{fixed_vcf_name} --index-output-vcf true ~{sep=' ' standardize_vcf_args}
    >>>

    output {
        File fixed_vcf = fixed_vcf_name
        File fixed_vcf_index = index_file_name
    }
}


task TrainGqRecalibratorTask {
    input {
        File train_vcf
        File train_vcf_index
        File ped_file
        File truth_file
        Array[File] genome_tracks
        File? gq_recalibrator_model_file # can be passed to do extra rounds of training on existing model
        Array[String] train_args = []
        String gatk_docker
        Int? num_entries
        Float mem_scale_vcf_size = 25.2
        Float mem_scale_num_entries = "3.7e-7"
        Float mem_gb_overhead = 1.5
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
        preemptible: 3
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
          --pedigree ~{ped_file} \
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
