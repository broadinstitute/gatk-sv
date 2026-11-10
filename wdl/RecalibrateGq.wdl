version 1.0

import "Structs.wdl"
import "TrainGqRecalibrator.wdl" as TrainGqRecalibrator
import "TasksMakeCohortVcf.wdl" as tasks_cohort
import "Utils.wdl" as Utils

workflow RecalibrateGq {
    input {
        File vcf
        File vcf_index
        Int? recalibrate_records_per_shard
        Array[File] genome_tracks = []
        File gq_recalibrator_model_file
        Array[String] recalibrate_gq_args = []
        String gatk_docker
        String sv_base_mini_docker
        String sv_pipeline_docker
        Float recalibrate_gq_mem_gb_java = 9.0
        Float recalibrate_gq_mem_gb_overhead = 1.5
        RuntimeAttr? runtime_attr_recalibrate_gq
    }

    call tasks_cohort.ScatterVcf {
        input:
            vcf=vcf,
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
                mem_gb_overhead=recalibrate_gq_mem_gb_overhead,
                runtime_attr_override=runtime_attr_recalibrate_gq
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
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: mem_gb_java + mem_gb_overhead,
        disk_gb: round(100 + 2 * size([vcf, gq_recalibrator_model_file], "GiB") + size(genome_tracks, "GiB")),
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String input_vcf_filename = basename(vcf)
    String filtered_vcf_name = sub(sub(input_vcf_filename, ".gz", ""), ".vcf", "_gq_recalibrated.vcf.gz")

    runtime {
        docker: gatk_docker
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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

        # gatk indices still have problems, overwrite with tabix
        tabix -f ~{filtered_vcf_name}
    >>>

    output {
        File filtered_vcf = filtered_vcf_name
        File filtered_vcf_index = filtered_vcf_name + ".tbi"
    }
}
