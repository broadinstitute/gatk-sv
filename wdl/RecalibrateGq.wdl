version 1.0

import "TransferVcfAnnotations.wdl" as transfer_annot
import "TrainGqRecalibrator.wdl" as TrainGqRecalibrator
import "TasksMakeCohortVcf.wdl" as tasks_cohort
import "Utils.wdl" as Utils

workflow RecalibrateGq {
    input {
        File vcf
        File vcf_index
        File? annotations_vcf
        File? annotations_vcf_index
        String? info_keys_to_transfer
        String? format_keys_to_transfer
        Int? recalibrate_records_per_shard
        File? transfer_script
        Int? transfer_annotations_shard_size
        Array[File] genome_tracks
        File gq_recalibrator_model_file
        Boolean standardize_vcf = true
        Array[String] standardize_vcf_args = []
        Array[String] recalibrate_gq_args = []
        String samtools_cloud_docker
        String gatk_docker
        String sv_base_mini_docker
        String sv_pipeline_docker
        String sv_utils_docker
        Float recalibrate_gq_mem_gb_java = 9.0
        Float recalibrate_gq_mem_gb_overhead = 1.5
    }

    if(standardize_vcf) {
        call TrainGqRecalibrator.StandardizeVcfForGatk {
            input:
                vcf=vcf,
                standardize_vcf_args=standardize_vcf_args,
                sv_utils_docker=sv_utils_docker
        }
    }

    if(defined(annotations_vcf)) {
        call transfer_annot.TransferVcfAnnotations {
            input:
                vcf=select_first([StandardizeVcfForGatk.fixed_vcf, vcf]),
                vcf_with_annotations=select_first([annotations_vcf]),
                shard_size=select_first([transfer_annotations_shard_size, 20000]),
                info_keys_list=info_keys_to_transfer,
                format_keys_list=format_keys_to_transfer,
                transfer_script=transfer_script,
                prefix=basename(vcf, ".vcf.gz") + ".transfer_annot",
                sv_base_mini_docker=sv_base_mini_docker,
                sv_pipeline_docker=sv_pipeline_docker
        }
    }

    File vcf_ = select_first([TransferVcfAnnotations.annotated_vcf, StandardizeVcfForGatk.fixed_vcf, vcf])
    File vcf_index_ = select_first([TransferVcfAnnotations.annotated_vcf_index, StandardizeVcfForGatk.fixed_vcf_index,
                                    vcf_index])

    call tasks_cohort.ScatterVcf {
        input:
            vcf=select_first([TransferVcfAnnotations.annotated_vcf, StandardizeVcfForGatk.fixed_vcf, vcf]),
            records_per_shard=select_first([recalibrate_records_per_shard, 20000]),
            prefix=basename(vcf, ".vcf.gz") + ".scatter",
            sv_pipeline_docker=sv_pipeline_docker
    }

    scatter ( i in range(length(ScatterVcf.shards)) ) {
        call RecalibrateGqTask {
            input:
                vcf=ScatterVcf.shards[i],
                genome_tracks=genome_tracks,
                gq_recalibrator_model_file=gq_recalibrator_model_file,
                recalibrate_gq_args=recalibrate_gq_args,
                gatk_docker=gatk_docker,
                mem_gb_java=recalibrate_gq_mem_gb_java,
                mem_gb_overhead=recalibrate_gq_mem_gb_overhead
        }
    }

    call tasks_cohort.ConcatVcfs {
        input:
            vcfs=RecalibrateGqTask.filtered_vcf,
            naive=true,
            outfile_prefix=basename(vcf, ".vcf.gz") + ".gq_recalibrated",
            sv_base_mini_docker=sv_pipeline_docker
    }

    output {
        File filtered_vcf = ConcatVcfs.concat_vcf
        File filtered_vcf_index = ConcatVcfs.concat_vcf_idx
    }
}


task RecalibrateGqTask {
    input {
        File vcf
        Array[File] genome_tracks
        File gq_recalibrator_model_file
        Array[String] recalibrate_gq_args = []
        String gatk_docker
        Float mem_gb_java = 9.0
        Float mem_gb_overhead = 1.5
    }

    Int disk_gb = round(100 + 2 * size([vcf, gq_recalibrator_model_file], "GiB") + size(genome_tracks, "GiB"))
    Float mem_gb = mem_gb_java + mem_gb_overhead
    String input_vcf_filename = basename(vcf)
    String filtered_vcf_name = sub(sub(input_vcf_filename, ".gz", ""), ".vcf", "_gq_recalibrated.vcf.gz")

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

        ln -s ~{vcf} ~{input_vcf_filename}
        tabix ~{input_vcf_filename}

        mem_kb_java_actual=$(grep -m1 MemTotal /proc/meminfo \
                             | awk '{printf("%.0f\n", $2 - ~{mem_gb_overhead} * 1048576)}')

        gatk --java-options "-Xmx${mem_kb_java_actual}K" XGBoostMinGqVariantFilter \
          --mode "Filter" \
          --variant ~{input_vcf_filename} \
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
