version 1.0

import "Utils.wdl" as Utils
import "ExtractVcfGqFilterProperties.wdl" as ExtractVcfGqFilterProperties
import "RecalibrateGqShards.wdl" as RecalibrateGqShards
import "BenchmarkGqFilter.wdl" as BenchmarkGqFilter

workflow CrossValidateGqRecalibrator {
    input {
        File train_vcf
        File train_vcf_index
        String train_vcf_label
        File? pretrained_gq_recalibrator_model
        File training_truth_json
        File benchmark_truth_json
        File? ped_file  # used only for benchmarking
        Array[File] genome_tracks
        Int num_splits = 5
        String sv_utils_docker
        String sv_base_mini_docker
        String samtools_cloud_docker
        String gatk_recalibrator_docker
        # optional arguments for overriding default tool behaviors
        Array[String] standardize_vcf_args = []
        Array[String] train_args = []
        Array[String] recalibrate_gq_args = []
        Array[String] annotate_gq_args = []
        Array[String] benchmark_args = []
        Int new_pipeline_passing_score = 1
        String new_pipeline_score_property = "SL"
        Int old_pipeline_passing_score = 1
        String old_pipeline_score_property = "GQ"
        # any additional scores to compare these training results to:
        Array[ScoresDataSet] comparison_scores = []
        # optional values that may be passed from previous runs to decrease runtime while iterating on a solution:
        CrossValidationVcfs? pre_computed_cross_validation_vcfs
    }

    call ExtractVcfGqFilterProperties.ExtractVcfGqFilterProperties as ExtractWholeCohort {
        input:
            vcf=train_vcf,
            vcf_index=train_vcf_index,
            genome_tracks=genome_tracks,
            samtools_cloud_docker=samtools_cloud_docker,
            gatk_recalibrator_docker=gatk_recalibrator_docker,
            sv_utils_docker=sv_utils_docker
    }
    call TrainGqRecalibratorTask as TrainWholeCohort {
        input:
            properties_parquet_tar=select_first([ExtractWholeCohort.properties_parquet_tar]),
            truth_json=training_truth_json,
            pretrained_gq_recalibrator_model=pretrained_gq_recalibrator_model,
            train_args=train_args,
            num_samples=ExtractWholeCohort.num_samples,
            sv_utils_docker=sv_utils_docker
    }
    call RecalibrateGqShards.RecalibrateGqShards as RecalibrateWholeCohort {
        input:
            vcf_shards=ExtractWholeCohort.vcf_shards,
            properties_parquet_shards=ExtractWholeCohort.properties_parquet_shards,
            gq_recalibrator_model=TrainWholeCohort.gq_recalibrator_model,
            recalibrated_vcf_basename=sub(
                sub(basename(train_vcf), ".gz$", ""),
                ".vcf$",
                "_recalibrated"
            ),
            num_samples=ExtractWholeCohort.num_samples,
            recalibrate_gq_args=recalibrate_gq_args,
            annotate_gq_args=annotate_gq_args,
            sv_base_mini_docker=sv_base_mini_docker,
            sv_utils_docker=sv_utils_docker
    }

    if(!defined(pre_computed_cross_validation_vcfs)) {
        call MakeCrossValidationVcfs {
            input:
                vcf=train_vcf,
                truth_json=training_truth_json,
                num_splits=num_splits,
                sv_utils_docker=sv_utils_docker
        }
    }
    CrossValidationVcfs cross_validation_vcfs_ = select_first([
        MakeCrossValidationVcfs.cross_validation_vcfs, pre_computed_cross_validation_vcfs
    ])

    scatter(fold in range(length(cross_validation_vcfs_.train_vcfs))) {
        File train_fold_vcf = cross_validation_vcfs_.train_vcfs[fold]
        File train_fold_vcf_idx = cross_validation_vcfs_.train_vcf_indices[fold]
        File test_fold_vcf = cross_validation_vcfs_.test_vcfs[fold]
        File test_fold_idx = cross_validation_vcfs_.test_vcfs[fold]

        call ExtractVcfGqFilterProperties.ExtractVcfGqFilterProperties as ExtractFoldTrain {
            input:
                vcf=train_fold_vcf,
                vcf_index=train_fold_vcf_idx,
                genome_tracks=genome_tracks,
                samtools_cloud_docker=samtools_cloud_docker,
                gatk_recalibrator_docker=gatk_recalibrator_docker,
                sv_utils_docker=sv_utils_docker
        }
        call ExtractVcfGqFilterProperties.ExtractVcfGqFilterProperties as ExtractFoldTest {
            input:
                vcf=test_fold_vcf,
                vcf_index=test_fold_idx,
                collect_whole_vcf_parquet=false,
                genome_tracks=genome_tracks,
                samtools_cloud_docker=samtools_cloud_docker,
                gatk_recalibrator_docker=gatk_recalibrator_docker,
                sv_utils_docker=sv_utils_docker
        }
        call TrainGqRecalibratorTask as TrainFold {
            input:
                properties_parquet_tar=select_first([ExtractFoldTrain.properties_parquet_tar]),
                truth_json=training_truth_json,
                pretrained_gq_recalibrator_model=pretrained_gq_recalibrator_model,
                train_args=train_args,
                num_samples=ExtractFoldTrain.num_samples,
                sv_utils_docker=sv_utils_docker
         }
        call RecalibrateGqShards.RecalibrateGqShards as RecalibrateFold {
            input:
                vcf_shards=ExtractFoldTest.vcf_shards,
                properties_parquet_shards=ExtractFoldTest.properties_parquet_shards,
                gq_recalibrator_model=TrainFold.gq_recalibrator_model,
                recalibrated_vcf_basename=sub(
                    sub(basename(train_vcf), ".gz$", ""), ".vcf$", "_cross_validated"
                ),
                num_samples=ExtractFoldTest.num_samples,
                recalibrate_gq_args=recalibrate_gq_args,
                annotate_gq_args=annotate_gq_args,
                sv_base_mini_docker=sv_base_mini_docker,
                sv_utils_docker=sv_utils_docker
        }
    }

    call MergeCrossvalidatedVcf {
        input:
            filtered_vcfs=RecalibrateFold.recalibrated_vcf,
            filtered_vcf_indices=RecalibrateFold.recalibrated_vcf_index,
            merged_name=sub(
                sub(basename(train_vcf), ".gz$", ""), ".vcf$", "_cross_validated.vcf.gz"
            ),
            sv_base_mini_docker=sv_base_mini_docker
    }

    Array[ScoresDataSet] comparison_scores_ = flatten(
        [
            [
                {
                    "label": "original",
                    "vcf": train_vcf,
                    "vcf_index": train_vcf_index,
                    "score_property": old_pipeline_score_property,
                    "passing_score": old_pipeline_passing_score
                },

                {
                    "label": "recalibrated",
                    "vcf": RecalibrateWholeCohort.recalibrated_vcf,
                    "vcf_index": RecalibrateWholeCohort.recalibrated_vcf_index,
                    "score_property": new_pipeline_score_property,
                    "passing_score": new_pipeline_passing_score
                },
                {
                    "label": "cross-validated",
                    "vcf": MergeCrossvalidatedVcf.merged_vcf,
                    "vcf_index": MergeCrossvalidatedVcf.merged_vcf_index,
                    "score_property": new_pipeline_score_property,
                    "passing_score": new_pipeline_passing_score
                }
            ],
            comparison_scores
        ]
    )

    # actually call BenchmarkFilter workflow
    call BenchmarkGqFilter.BenchmarkGqFilter {
        input:
            data_label=train_vcf_label,
            comparison_scores=comparison_scores_,
            benchmark_truth_json=benchmark_truth_json,
            ped_file=ped_file,
            benchmark_args=benchmark_args,
            samtools_cloud_docker=samtools_cloud_docker,
            sv_utils_docker=sv_utils_docker,
            sv_utils_docker=sv_utils_docker
    }

    output {
        File original_scores = BenchmarkGqFilter.comparison_scores_parquet_tar[0]
        File gq_recalibrator_model = TrainWholeCohort.gq_recalibrator_model
        File filtered_vcf = RecalibrateWholeCohort.recalibrated_vcf
        File filtered_vcf_index = RecalibrateWholeCohort.recalibrated_vcf_index
        File filtered_scores = BenchmarkGqFilter.comparison_scores_parquet_tar[1]
        File cross_validated_vcf = MergeCrossvalidatedVcf.merged_vcf
        File cross_validated_vcf_index = MergeCrossvalidatedVcf.merged_vcf_index
        File cross_validated_scores = BenchmarkGqFilter.comparison_scores_parquet_tar[2]
        File benchmark_figure = BenchmarkGqFilter.benchmark_figure
        Array[File] comparison_scores_parquet_tars = BenchmarkGqFilter.comparison_scores_parquet_tar
        CrossValidationVcfs cross_validation_vcfs = cross_validation_vcfs_
    }
}


task TrainGqRecalibratorTask {
    input {
        File properties_parquet_tar
        File truth_json
        File? pretrained_gq_recalibrator_model # can be passed to do extra rounds of training on existing model
        Array[String] train_args = []
        String sv_utils_docker
        Int num_samples
        Float mem_gb_base = 1.5
        Float mem_gb_per_sample = 0.01
        Int num_cpu = 4
    }

    Float mem_gb = mem_gb_base + mem_gb_per_sample * num_samples
    Int disk_gb = round(50 + 2 * size([properties_parquet_tar], "GiB") + size(truth_json, "GiB"))
    String model_file_name = if defined(pretrained_gq_recalibrator_model)
        then basename(select_first([pretrained_gq_recalibrator_model]))
        else "gq_recalibrator.model"
    # save intermediate results to checkpoint file, allowing resume from preemption to incur
    # little cost in wall-clock time
    String checkpoint_file_name = "gq_recalibrator_checkpoint.model"

    runtime {
        docker: sv_utils_docker
        cpu: num_cpu
        preemptible: 10  # using checkpointFile, preemption has little cost
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
        checkpointFile: checkpoint_file_name
    }

    command <<<
        set -euo pipefail

        # in order of preference: resume from checkpoint, resume from pre-trained model, start new
        if [[ -s ~{checkpoint_file_name} ]]; then
            echo "Resuming from non-empty checkpoint file"
            INPUT_ARGS="--input-model ~{checkpoint_file_name}"
        elif ~{defined(pretrained_gq_recalibrator_model)}; then
            echo "Starting from pretrained model: ~{pretrained_gq_recalibrator_model}"
            INPUT_ARGS="--input-model ~{pretrained_gq_recalibrator_model}"
        else
            echo "Starting training from scratch"
            INPUT_ARGS=""
        fi

        # train, periodically saving to checkpoint
        gq-recalibrator train-sv-gq-recalibrator \
            --properties "~{properties_parquet_tar}" \
            --truth-json "~{truth_json}" \
            --model "~{checkpoint_file_name}" \
            --torch-device cuda \
            --temp-dir "$(mktemp -d --tmpdir=.)" \
            $INPUT_ARGS \
            ~{sep=' ' train_args}

        # copy checkpoint file to final model file
        cp "~{checkpoint_file_name}" "~{model_file_name}"
    >>>

    output {
        File gq_recalibrator_model = model_file_name
    }
}

struct CrossValidationVcfs {
    Array[File] train_vcfs
    Array[File] train_vcf_indices
    Array[File] test_vcfs
    Array[File] test_vcf_indices
}

task MakeCrossValidationVcfs {
    input {
        File vcf
        File truth_json
        Int num_splits = 5
        String sv_utils_docker
    }

    String fixed_vcf_name = sub(sub(basename(vcf), ".gz$", ""), ".vcf$", "_fixed.vcf.gz")
    String index_file_name = fixed_vcf_name + ".tbi"

    Int disk_gb = 50 + round((1 + num_splits) * size(vcf, "GiB") + size(truth_json, "GiB"))
    Float mem_gb = 2.0

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

        # create TRUTH_SAMPLES_FILE, with one sample ID per line that has truth data (duplicates are okay)
        # The cross-validated VCFs will attempt to distribute these truth samples evenly across batches
        { python - <<'CODE'
import json
input_json="~{truth_json}"
with open(input_json, 'r') as f_in:
    truth_data = json.load(f_in)
    for sample_id in truth_data.keys():
        print(sample_id)
CODE
} > "$TRUTH_SAMPLES_FILE"

        sv-utils make-cross-validation-vcfs ~{vcf} \
            --truth-samples-file "$TRUTH_SAMPLES_FILE" \
            --num-splits ~{num_splits} \
            --index-output-vcf true
    >>>

    output {
        CrossValidationVcfs cross_validation_vcfs = {
            "train_vcfs": glob("train_*.vcf.gz"),
            "train_vcf_indices": glob("train_*.vcf.gz.tbi"),
            "test_vcfs": glob("test_*.vcf.gz"),
            "test_vcf_indices": glob("test_*.vcf.gz.tbi")
        }
    }
}


task MergeCrossvalidatedVcf {
    input {
        Array[File] filtered_vcfs
        Array[File] filtered_vcf_indices
        String merged_name
        String sv_base_mini_docker
    }

    Int disk_gb = 50 + 2 * round(size(filtered_vcfs, "GiB") + size(filtered_vcf_indices, "GiB"))
    Float mem_gb = 4.0
    String merged_index = merged_name + ".tbi"

    runtime {
        docker: sv_base_mini_docker
        cpu: 1
        preemptible: 3
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -euo pipefail

        bcftools merge --threads $(nproc) -m id -O z -o ~{merged_name} ~{sep=" " filtered_vcfs}

        bcftools index --tbi -o ~{merged_index} ~{merged_name}
    >>>

    output {
        File merged_vcf = merged_name
        File merged_vcf_index = merged_index
    }
}
