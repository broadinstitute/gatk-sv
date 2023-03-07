version 1.0


task TrainGqRecalibratorTask {
    input {
        File properties_parquet_tar
        File truth_json
        File? pretrained_gq_recalibrator_model # can be passed to do extra rounds of training on existing model
        Array[String] train_args = []
        String gq_recalibrator_docker
        Float mem_gb = 8
        Float mem_gb_overhead = 1.5
    }

    Int disk_gb = round(50 + 2 * size([properties_parquet_tar], "GiB") + size(truth_json, "GiB"))
    String model_file_name = if defined(pretrained_gq_recalibrator_model)
        then basename(select_first([pretrained_gq_recalibrator_model]))
        else "gq_recalibrator.model"

    runtime {
        docker: gq_recalibrator_docker
        cpu: 1
        preemptible: 3
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -euo pipefail
        if ~{defined(pretrained_gq_recalibrator_model)}; then
            INPUT_ARGS="-i ~{pretrained_gq_recalibrator_model}"
        else
            INPUT_ARGS=""
        fi

        gq-recalibrator train-sv-gq-recalibrator \
            --properties "~{properties_parquet_tar}" \
            --truth-json "$TRUTH_FILE" \
            --torch-device cuda \
            --temp-dir "$(mktemp -d --tmpdir=.)"
            $INPUT_ARGS \
            -m "~{model_file_name}" \
            ~{sep=' ' train_args}
    >>>

    output {
        File gq_recalibrator_model = model_file_name
    }
}

task RecalibrateGqTask {
    input {
        File properties_parquet_tar
        File gq_recalibrator_model
        Array[String] recalibrate_gq_args = []
        String gq_recalibrator_docker
        Float mem_gb = 8
        Float mem_gb_overhead = 1.5
    }

    Int disk_gb = round(50 + 2 * size([properties_parquet_tar], "GiB"))
    String recalibrated_properties_name = sub(
        basename(properties_parquet_tar), ".tar$", "-recalibrated.tar"
    )

    runtime {
        docker: gq_recalibrator_docker
        cpu: 1
        preemptible: 3
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -euo pipefail

        gq-recalibrator recalibrate-gq \
            --properties "~{properties_parquet_tar}" \
            -m "~{gq_recalibrator_model}" \
            -o "~{recalibrated_properties_name}" \
            --temp-dir "$(mktemp -d --tmpdir=.)"
            ~{sep=' ' recalibrate_gq_args}
    >>>

    output {
        File recalibrated_properties_parquet_tar = recalibrated_properties_name
    }
}

task AnnotateRecalibratedGqs {
    input {
        File vcf
        File scores_parquet_tar
        Array[String] annotate_recalibrated_gq_args = []
        String gq_recalibrator_docker
        Float mem_gb = 8
        Float mem_gb_overhead = 1.5
    }

    Int disk_gb = round(50 + 2 * size([properties_parquet_tar], "GiB"))
    String annotated_vcf_name = sub(sub(basename(vcf), ".gz$", ""), ".vcf$", "-recalibrated.vcf.gz")

    runtime {
        docker: gq_recalibrator_docker
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
            ~{sep=' ' annotate_recalibrated_gq_args}
    >>>

    output {
        File recalibrated_vcf = annotated_vcf_name
    }
}
