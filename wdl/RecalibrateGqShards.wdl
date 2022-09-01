version 1.0

import "TasksMakeCohortVcf.wdl" as TasksMakeCohortVcf
import "ExtractVcfGqFilterProperties.wdl" as ExtractVcfGqFilterProperties

workflow RecalibrateGqShards {
    input {
        Array[File] vcf_shards
        Array[File] properties_parquet_shards
        File gq_recalibrator_model
        String recalibrated_vcf_basename  # file name with no .vcf.gz suffix
        Int num_samples
        Array[String] recalibrate_gq_args = []
        Array[String] annotate_gq_args = []
        String sv_base_mini_docker
        String sv_utils_docker
    }

    scatter(shard_index in range(length(vcf_shards))) {
        File vcf_shard = vcf_shards[shard_index]
        File properties_parquet_shard = properties_parquet_shards[shard_index]
        call RecalibrateGqTask {
            input:
                properties_parquet_tar=properties_parquet_shard,
                gq_recalibrator_model=gq_recalibrator_model,
                recalibrate_gq_args=recalibrate_gq_args,
                sv_utils_docker=sv_utils_docker,
                num_samples=num_samples
        }
        call AnnotateRecalibratedGqs {
            input:
                vcf=vcf_shard,
                scores_parquet_tar=RecalibrateGqTask.recalibrated_scores_parquet_tar,
                annotate_gq_args=annotate_gq_args,
                sv_utils_docker=sv_utils_docker
        }
    }

    call ExtractVcfGqFilterProperties.ConcatenateParquetTars {
        input:
            parquet_tars=RecalibrateGqTask.recalibrated_scores_parquet_tar,
            output_base_name=recalibrated_vcf_basename,
            sv_utils_docker=sv_utils_docker,
            compute_weights=false
    }

    call TasksMakeCohortVcf.ConcatVcfs {
        input:
            vcfs=AnnotateRecalibratedGqs.recalibrated_vcf,
            naive=true,
            generate_index=true,
            outfile_prefix=recalibrated_vcf_basename,
            sv_base_mini_docker=sv_base_mini_docker
    }

    output {
        File recalibrated_scores_parquet = ConcatenateParquetTars.combined_parquet_tar
        File recalibrated_vcf = ConcatVcfs.concat_vcf
        File recalibrated_vcf_index = ConcatVcfs.concat_vcf_idx
    }
}

task RecalibrateGqTask {
    input {
        File properties_parquet_tar
        File gq_recalibrator_model
        Array[String] recalibrate_gq_args = []
        String sv_utils_docker
        Int num_samples
        Float mem_baseline_gb = 1.0
        Float mem_gb_per_sample = 0.005
    }

    Float mem_gb = mem_baseline_gb + mem_gb_per_sample * num_samples
    Int disk_gb = round(50 + 3 * size([properties_parquet_tar], "GiB"))
    String recalibrated_scores_name = sub(
        basename(properties_parquet_tar), ".tar$", "-recalibrated.tar"
    )

    runtime {
        docker: sv_utils_docker
        cpu: 4
        preemptible: 3
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -euo pipefail

        gq-recalibrator recalibrate-gq \
            --properties "~{properties_parquet_tar}" \
            --model "~{gq_recalibrator_model}" \
            --output-parquet "~{recalibrated_scores_name}" \
            --temp-dir "$(mktemp -d --tmpdir=.)"
            ~{sep=' ' recalibrate_gq_args}
    >>>

    output {
        File recalibrated_scores_parquet_tar = recalibrated_scores_name
    }
}

task AnnotateRecalibratedGqs {
    input {
        File vcf
        File scores_parquet_tar
        Array[String] annotate_gq_args = []
        String sv_utils_docker
        Float mem_gb = 24
    }

    Int disk_gb = round(50 + 2 * size([scores_parquet_tar], "GiB"))
    String annotated_vcf_name = sub(sub(basename(vcf), ".gz$", ""), ".vcf$", "-recalibrated.vcf.gz")

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

        gq-recalibrator annotate-recalibrated-gq \
            -i "~{vcf}" \
            -a "~{scores_parquet_tar}" \
            -o "~{annotated_vcf_name}" \
            --temp-dir "$(mktemp -d --tmpdir=.)"
            ~{sep=' ' annotate_gq_args}
    >>>

    output {
        File recalibrated_vcf = annotated_vcf_name
        File recalibrated_vcf_index = annotated_vcf_name + ".tbi"
    }
}
