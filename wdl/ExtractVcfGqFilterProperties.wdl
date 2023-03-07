version 1.0

import "Utils.wdl" as Utils
import "TasksMakeCohortVcf.wdl" as TasksMakeCohortVcf

workflow ExtractVcfGqFilterProperties {
    input {
        File vcf
        File vcf_index
        Array[File] genome_tracks
        Boolean collect_whole_vcf_parquet = true
        Array[String] tarred_properties_to_parquet_args = []
        String samtools_cloud_docker
        String sv_utils_docker
        String gq_recalibrator_docker
        Int records_per_shard = 10000
    }

    call TasksMakeCohortVcf.ScatterVcf {
        input:
            vcf=vcf,
            vcf_index=vcf_index,
            records_per_shard=records_per_shard,
            sv_pipeline_docker=samtools_cloud_docker
    }

    scatter(vcf_shard in ScatterVcf.shards) {
        call Utils.IndexVcf {
            input:
                vcf=vcf_shard,
                samtools_cloud_docker=samtools_cloud_docker
        }
        call ExtractVcfProperties {
            input:
                vcf=vcf_shard,
                vcf_index=vcf_index,
                genome_tracks=genome_tracks,
                sv_utils_docker=sv_utils_docker
        }

        File tsv_tar_shard = ExtractVcfProperties.properties_tar_file
        call TsvPropertiesToParquet as ShardParquet {
            input:
                properties_tar_files=[tsv_tar_shard],
                output_base_name=sub(basename(tsv_tar_shard), ".properties.tar$", ""),
                tarred_properties_to_parquet_args=tarred_properties_to_parquet_args,
                gq_recalibrator_docker=gq_recalibrator_docker
        }
    }

    if(collect_whole_vcf_parquet) {
        # get parquet file for whole data set
        call ConcatenateParquetTars as ParquetVcf {
            input:
                parquet_tars=ShardParquet.properties_parquet_tar,
                output_base_name=sub(sub(basename(vcf), ".gz$", ""), ".vcf$", ""),
                gq_recalibrator_docker=gq_recalibrator_docker
        }
    }

    output {
        Array[File] vcf_shards = ScatterVcf.shards
        Array[File] properties_parquet_shards = ShardParquet.properties_parquet_tar
        File? properties_parquet_tar = ParquetVcf.combined_parquet_tar
    }
}

task ExtractVcfProperties {
    input {
        File vcf
        File vcf_index
        Array[File] genome_tracks
        String sv_utils_docker
        Float mem_gb_java = 3.5
        Float mem_gb_overhead = 1.5
    }

    Float mem_gb = mem_gb_java + mem_gb_overhead
    Int disk_gb = round(50 + 2 * size([vcf, vcf_index], "GiB") + size(genome_tracks, "GiB"))

    runtime {
        docker: sv_utils_docker
        cpu: 1
        preemptible: 3
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    String output_folder_name = sub(sub(basename(vcf), ".gz$", ""), ".vcf$", ".properties")
    String output_file_name = output_folder_name + ".tar"
    String properties_summary_json = "properties_summary.json"

    command <<<
        set -euo pipefail

        mem_kb_java_actual=$(grep -m1 MemTotal /proc/meminfo \
                             | awk '{printf("%.0f\n", $2 - ~{mem_gb_overhead} * 1048576)}')

        NUM_PHYSICAL_CPUS=$(grep ^"core id" /proc/cpuinfo | sort -u | wc -l)

        gatk --java-options "-Xmx${mem_kb_java_actual}K"  ExtractSvProperties \
            --variant "~{vcf}" \
            -O "~{output_folder_name}" \
            ~{if len(genome_tracks) > 0 then "--genome-track" else ""} ~{sep=" --genome-track " genome_tracks}

        mv "~{output_folder_name}/~{properties_summary_json}" .

        tar --create --remove-files \
            --file "~{output_file_name}" \
            "~{output_folder_name}"
    >>>


    output {
        File properties_tar_file = output_file_name
        File properties_summary = properties_summary_json
    }
}

task TsvPropertiesToParquet {
    input {
        Array[File] properties_tar_files
        String output_base_name
        Array[String] tarred_properties_to_parquet_args = []
        String gq_recalibrator_docker
        Float mem_gb = 8
    }

    Float decompression_factor = 10
    Int disk_gb = round(50 + (2 + decompression_factor) * size(properties_tar_files, "GiB"))

    runtime {
        docker: gq_recalibrator_docker
        cpu: 1
        preemptible: 3
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    String properties_parquet_filename = output_base_name + ".parquet.tar"
    String properties_scaling_json_filename = output_base_name + ".scaling.json"

    command <<<
        set -euo pipefail

        gq-recalibrator tarred-properties-to-parquet \
            --input-tar ~{sep=" --input-tar " properties_tar_files} \
            --output-parquet "~{properties_parquet_filename}" \
            --temp-dir "$(mktemp -d --tmpdir=.)" \
            ~{sep=' ' tarred_properties_to_parquet_args}
    >>>

    output {
        File properties_parquet_tar = properties_parquet_filename
    }
}

task ConcatenateParquetTars {
    input {
        Array[File] parquet_tars
        String output_base_name
        String gq_recalibrator_docker
        Float mem_gb = 8
    }

    Int disk_gb = round(50 + 2 * size(parquet_tars, "GiB"))

    runtime {
        docker: gq_recalibrator_docker
        cpu: 1
        preemptible: 3
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    String combined_parquet_filename = output_base_name + ".parquet.tar"

    command <<<
        set -euo pipefail

        gq-recalibrator concat-parquet-shards \
            --input-parquet ~{sep="," parquet_tars} \
            --output-parquet "~{combined_parquet_filename}" \
            --temp-dir "$(mktemp -d --tmpdir=.)" \
    >>>

    output {
        File combined_parquet_tar = combined_parquet_filename
    }
}

