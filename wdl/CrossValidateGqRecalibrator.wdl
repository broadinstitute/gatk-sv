version 1.0

import "Utils.wdl" as Utils
import "ExtractVcfGqFilterProperties.wdl" as ExtractVcfGqFilterProperties
import "RecalibrateGq.wdl" as RecalibrateGq
import "TrainGqRecalibrator.wdl" as TrainGqRecalibrator
import "BenchmarkGqFilter.wdl" as BenchmarkGqFilter

workflow CrossValidateGqRecalibrator {
    input {
        File train_vcf
        File train_vcf_index
        File? annotations_vcf
        File? annotations_vcf_index
        Array[String]? annotations_to_transfer
        File? pretrained_gq_recalibrator_model
        File training_truth_json
        File benchmark_truth_json
        String train_vcf_label
        Array[File] truth_vcfs = []
        Array[File] truth_vcf_indices = []
        Array[String]? vapor_sample_ids
        Array[File]? vapor_files
        File? ped_file
        Array[File] genome_tracks
        File? optimal_overlap_cutoffs
        Boolean standardize_vcf = true
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
        Array[String] get_truth_overlap_args = []
        Array[String] benchmark_args = []
        Int new_pipeline_passing_score = 1
        String new_pipeline_score_property = "SL"
        Int old_pipeline_passing_score = 1
        String old_pipeline_score_property = "OGQ"
        # any additional scores to compare these training results to:
        Array[ScoresDataSet] comparison_scores = []
        # optional values that may be passed from previous runs to decrease runtime while iterating on a solution:
        CrossValidationVcfs? pre_computed_cross_validation_vcfs
    }

    if(standardize_vcf) {
        call TrainGqRecalibrator.StandardizeVcfForGatk {
            input:
                vcf=train_vcf,
                standardize_vcf_args=standardize_vcf_args,
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

    File train_vcf_ = select_first([
        TransferVcfAnnotations.annotated_vcf, StandardizeVcfForGatk.fixed_vcf, train_vcf
    ])
    File train_vcf_index_ = select_first([
        TransferVcfAnnotations.annotated_vcf_index,
        StandardizeVcfForGatk.fixed_vcf_index,
        train_vcf_index
    ])
    call ExtractVcfGqFilterProperties.ExtractVcfGqFilterProperties as ExtractWholeCohort {
        input:
            vcf=train_vcf_,
            vcf_index=train_vcf_index_,
            genome_tracks=genome_tracks,
            samtools_cloud_docker=samtools_cloud_docker,
            gatk_recalibrator_docker=gatk_recalibrator_docker,
            sv_utils_docker=sv_utils_docker
    }
    call GqRecalibratorTasks.TrainGqRecalibratorTask as TrainWholeCohort {
        input:
            properties_parquet_tar=select_first([ExtractWholeCohort.properties_parquet_tar]),
            truth_json=training_truth_json,
            pretrained_gq_recalibrator_model=pretrained_gq_recalibrator_model,
            train_args=train_args,
            sv_utils_docker=sv_utils_docker
    }
    call RecalibrateGq.RecalibrateGq as RecalibrateWholeCohort {
        input:
            vcf_shards=ExtractWholeCohort.vcf_shards,
            properties_parquet_shards=ExtractWholeCohort.properties_parquet_shards,
            gq_recalibrator_model=TrainWholeCohort.gq_recalibrator_model,
            recalibrated_vcf_basename=sub(
                sub(basename(train_vcf_), ".gz$", ""), ".vcf$", "-recalibrated.vcf.gz"
            ),
            recalibrate_gq_args=recalibrate_gq_args,
            annotate_gq_args=annotate_gq_args,
            sv_base_mini_docker=sv_base_mini_docker,
            sv_utils_docker=sv_utils_docker
    }

    if(!defined(pre_computed_cross_validation_vcfs)) {
        call MakeCrossValidationVcfs {
            input:
                vcf=train_vcf_,
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
                collect_shard_parquets=false,
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
        call GqRecalibratorTasks.TrainGqRecalibratorTask as TrainFold {
            input:
                properties_parquet_tar=select_first([ExtractFoldTrain.properties_parquet_tar]),
                truth_json=training_truth_json,
                pretrained_gq_recalibrator_model=pretrained_gq_recalibrator_model,
                train_args=train_args,
                sv_utils_docker=sv_utils_docker
        }
        call RecalibrateGq.RecalibrateGq as RecalibrateFold {
            input:
                vcf_shards=ExtractFoldTest.vcf_shards,
                properties_parquet_shards=ExtractFoldTest.properties_parquet_shards,
                gq_recalibrator_model=TrainFold.gq_recalibrator_model,
                recalibrated_vcf_basename=sub(
                    sub(basename(train_vcf_), ".gz$", ""), ".vcf$", "-recalibrated.vcf.gz"
                ),
                recalibrate_gq_args=recalibrate_gq_args,
                annotate_gq_args=annotate_gq_args,
                sv_base_mini_docker=sv_base_mini_docker,
                sv_utils_docker=sv_utils_docker
        }
    }

    SHOULD WRITE MERGE RECALIBRATED TEST PARQUET or similar name
    Note this is concatenating along axis 1!

    call MergeRecalibratedTestVcfs {
        input:
            filtered_vcfs=RecalibrateFold.recalibrated_vcf,
            filtered_vcf_indices=RecalibrateFold.recalibrated_vcf_index,
            merged_name=sub(
                sub(basename(train_vcf), ".gz$", ""), ".vcf$", "_cross_validated.vcf.gz"
            ),
            sv_base_mini_docker=sv_base_mini_docker
    }

    ScoresDataSet original_scores = {
        "label": "original",
        "vcf": train_vcf_,
        "vcf_index": train_vcf_index_,
        "scores_parquet_tar": RecalibrateWholeCohort.recalibrated_scores_parquet,
        "property": old_pipeline_score_property,
        "passing_score": old_pipeline_passing_score
    }


    Array[ScoresDataSet] comparison_scores_ = flatten(
        [
            [
                {
                    "label": "recalibrated",
                    "vcf": DirectGqRecalibrator.filtered_vcf,
                    "vcf_index": DirectGqRecalibrator.filtered_vcf_index,
                    "scores_parquet_tar": RecalibrateWholeCohort.recalibrated_scores_parquet
                    "property": new_pipeline_score_property,
                    "passing_score": new_pipeline_passing_score
                },
                {
                    "label": "cross-validated",
                    "vcf": MergeRecalibratedTestVcfs.merged_vcf,
                    "vcf_index": MergeRecalibratedTestVcfs.merged_vcf_index,
                    "property": new_pipeline_score_property,
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
            original_scores=original_scores,
            variant_properties_parquet_tar=select_first([ExtractWholeCohort.properties_parquet_tar]),
            comparison_scores=comparison_scores_,
            benchmark_truth_json=benchmark_truth_json,
            ped_file=ped_file,
            benchmark_args=benchmark_args,
            sv_utils_docker=sv_utils_docker
    }

    output {
        File clean_vcf = TrainGqRecalibrator.clean_vcf
        File clean_vcf_index = TrainGqRecalibrator.clean_vcf_index
        File truth_overlap_info = truth_overlap_info_
        Array[File] output_optimal_overlap_cutoffs = TrainGqRecalibrator.output_optimal_overlap_cutoffs
        File output_gq_recalibrator_model_file = TrainGqRecalibrator.output_gq_recalibrator_model_file
        File filtered_vcf = DirectGqRecalibrator.filtered_vcf
        File filtered_vcf_index = DirectGqRecalibrator.filtered_vcf_index
        File cross_validated_filtered_vcf = MergeRecalibratedTestVcfs.merged_vcf
        File cross_validated_filtered_vcf_index = MergeRecalibratedTestVcfs.merged_vcf_index
        File benchmark_figure = BenchmarkGqFilter.benchmark_figure
        File variant_properties = BenchmarkGqFilter.variant_properties
        File pickled_original_scores = BenchmarkGqFilter.pickled_original_scores
        Array[File] pickled_comparison_scores = BenchmarkGqFilter.pickled_comparison_scores
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
        Float mem_gb = 8
        Float mem_gb_overhead = 1.5
    }

    Int disk_gb = round(50 + 2 * size([properties_parquet_tar], "GiB") + size(truth_json, "GiB"))
    String model_file_name = if defined(pretrained_gq_recalibrator_model)
        then basename(select_first([pretrained_gq_recalibrator_model]))
        else "gq_recalibrator.model"

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

    Int disk_gb = 50 + round(
        (1 + num_splits) * size(vcf, "GiB") + size(truth_vcfs, "GiB") + size(ped_file, "GiB")
    )
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
    truth_overlap = json.load(f_in)
    for sample_id in truth_overlap.keys():
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


task MergeRecalibratedTestVcfs {
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


    output {
        File merged_parquet_tar = output_base_name + ".parquet.tar"
    }
}
