version 1.0

import "TrainGqRecalibrator.wdl" as TrainGqRecalibrator
import "RecalibrateGq.wdl" as RecalibrateGq
import "BenchmarkGqFilter.wdl" as BenchmarkGqFilter
import "PickleVcfProperties.wdl" as PickleVcfProperties
import "Utils.wdl" as Utils

workflow CrossValidateGqRecalibrator {
    input {
        File train_vcf
        File train_vcf_index
        File? annotations_vcf
        File? annotations_vcf_index
        Array[String]? annotations_to_transfer
        String train_vcf_label
        Array[File] truth_vcfs
        Array[File] truth_vcf_indices
        Array[String]? vapor_sample_ids
        Array[File]? vapor_files
        File ped_file
        Array[File] genome_tracks
        File? optimal_overlap_cutoffs
        Boolean standardize_vcf = true
        Int num_splits = 5
        String sv_utils_docker
        String gatk_docker
        String sv_base_mini_docker
        String samtools_cloud_docker
        # optional arguments for overriding default tool behaviors
        Array[String] standardize_vcf_args = []
        Array[String] train_args = []
        Array[String] recalibrate_gq_args = []
        Array[String] get_truth_overlap_args = []
        Array[String] benchmark_args = []
        Int new_pipeline_passing_score = 1
        String new_pipeline_score_property = "sl"
        Int old_pipeline_passing_score = 1
        String old_pipeline_score_property = "gq"
        # any additional scores to compare these training results to:
        Array[ScoresDataSet] comparison_scores = []
        # optional values that may be passed from previous runs to decrease runtime while iterating on a solution:
        File? pre_computed_pickled_variant_properties
        File? pre_computed_pickled_original_scores
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

    File train_vcf_ = select_first([TransferVcfAnnotations.annotated_vcf, StandardizeVcfForGatk.fixed_vcf, train_vcf])
    File train_vcf_index_ = select_first([TransferVcfAnnotations.annotated_vcf_index,
                                          StandardizeVcfForGatk.fixed_vcf_index, train_vcf_index])

    call TrainGqRecalibrator.TrainGqRecalibrator {
        input:
            train_vcf=train_vcf_,
            train_vcf_index=train_vcf_index_,
            truth_vcfs=truth_vcfs,
            truth_vcf_indices=truth_vcf_indices,
            vapor_sample_ids=vapor_sample_ids,
            vapor_files=vapor_files,
            ped_file=ped_file,
            genome_tracks=genome_tracks,
            optimal_overlap_cutoffs=optimal_overlap_cutoffs,
            standardize_vcf=false,
            train_args=train_args,
            get_truth_overlap_args=get_truth_overlap_args,
            sv_utils_docker=sv_utils_docker,
            gatk_docker=gatk_docker,
            samtools_cloud_docker=samtools_cloud_docker
    }

    call RecalibrateGq.RecalibrateGq as DirectGqRecalibrator {
        input:
            vcf=train_vcf_,
            vcf_index=train_vcf_index_,
            genome_tracks=genome_tracks,
            gq_recalibrator_model_file=TrainGqRecalibrator.output_gq_recalibrator_model_file,
            standardize_vcf=false,
            recalibrate_gq_args=recalibrate_gq_args,
            gatk_docker=gatk_docker,
            sv_utils_docker=sv_utils_docker,
            samtools_cloud_docker=samtools_cloud_docker
    }

    if(!defined(pre_computed_cross_validation_vcfs)) {
        call MakeCrossValidationVcfs {
            input:
                vcf=train_vcf_,
                ped_file=ped_file,
                truth_vcfs=truth_vcfs,
                vapor_sample_ids=if defined(vapor_files) then select_first([vapor_sample_ids]) else [],
                num_splits=num_splits,
                sv_utils_docker=sv_utils_docker
        }
    }
    CrossValidationVcfs cross_validation_vcfs_ = select_first([MakeCrossValidationVcfs.cross_validation_vcfs,
                                                               pre_computed_cross_validation_vcfs])

    scatter(fold in range(length(cross_validation_vcfs_.train_vcfs))) {
        call TrainGqRecalibrator.TrainGqRecalibrator as CrossTrainGqRecalibrator {
            input:
                train_vcf=cross_validation_vcfs_.train_vcfs[fold],
                train_vcf_index=cross_validation_vcfs_.train_vcf_indices[fold],
                truth_vcfs=truth_vcfs,
                truth_vcf_indices=truth_vcf_indices,
                vapor_sample_ids=vapor_sample_ids,
                vapor_files=vapor_files,
                ped_file=ped_file,
                genome_tracks=genome_tracks,
                optimal_overlap_cutoffs=optimal_overlap_cutoffs,
                standardize_vcf=false,
                train_args=train_args,
                get_truth_overlap_args=get_truth_overlap_args,
                sv_utils_docker=sv_utils_docker,
                gatk_docker=gatk_docker,
                samtools_cloud_docker=samtools_cloud_docker
        }

        call RecalibrateGq.RecalibrateGq as CrossGqRecalibrator {
            input:
                vcf=cross_validation_vcfs_.test_vcfs[fold],
                vcf_index=cross_validation_vcfs_.test_vcf_indices[fold],
                genome_tracks=genome_tracks,
                gq_recalibrator_model_file=CrossTrainGqRecalibrator.output_gq_recalibrator_model_file,
                standardize_vcf=false,
                recalibrate_gq_args=recalibrate_gq_args,
                gatk_docker=gatk_docker,
                sv_utils_docker=sv_utils_docker,
                samtools_cloud_docker=samtools_cloud_docker
        }
    }

    call MergeRecalibratedTestVcfs {
        input:
            filtered_vcfs=CrossGqRecalibrator.filtered_vcf,
            filtered_vcf_indices=CrossGqRecalibrator.filtered_vcf_index,
            merged_name=sub(sub(basename(train_vcf), ".gz$", ""), ".vcf$", "_cross_validated.vcf.gz"),
            sv_base_mini_docker=sv_base_mini_docker
    }

    # cromwell doesn't launch tasks until an entire workflow is ready, so can speed up wall-clock time by getting
    # variant properties and original scores as soon as train_vcf_ is resolved
    if(!defined(pre_computed_pickled_variant_properties)) {
        call PickleVcfProperties.PickleVcfProperties as PickleVariantData {
            input:
                vcf=train_vcf_,
                vcf_index=train_vcf_index_,
                wanted_properties=["svtype", "svlen", "ac", "is_autosome"],
                samtools_cloud_docker=samtools_cloud_docker,
                sv_utils_docker=sv_utils_docker
        }
    }
    File pickled_variant_properties_ = select_first([pre_computed_pickled_variant_properties,
                                                     PickleVariantData.pickled_properties])

    if(!defined(pre_computed_pickled_original_scores)) {
        call PickleVcfProperties.PickleVcfProperties as PickleTrainScores {
            input:
                vcf=train_vcf_,
                vcf_index=train_vcf_index_,
                wanted_properties=[old_pipeline_score_property],
                samtools_cloud_docker=samtools_cloud_docker,
                sv_utils_docker=sv_utils_docker
        }
    }
    File pickled_original_scores_ = select_first([pre_computed_pickled_original_scores,
                                                  PickleTrainScores.pickled_properties])

    ScoresDataSet original_scores = {
        "label": "original",
        "vcf": train_vcf_,
        "vcf_index": train_vcf_index_,
        "pickled_scores_file": pickled_original_scores_,
        "property": old_pipeline_score_property,
        "passing_score": old_pipeline_passing_score
    }

    Array[ScoresDataSet] scores_0 = [
        {
            "label": "recalibrated",
            "vcf": DirectGqRecalibrator.filtered_vcf,
            "vcf_index": DirectGqRecalibrator.filtered_vcf_index,
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
    ]
    Array[ScoresDataSet] comparison_scores_ = flatten([scores_0, comparison_scores])

    # actually call BenchmarkFilter workflow
    call BenchmarkGqFilter.BenchmarkGqFilter {
        input:
            data_label=train_vcf_label,
            original_scores=original_scores,
            pickled_variant_properties=pickled_variant_properties_,
            comparison_scores=comparison_scores_,
            truth_overlap_info=TrainGqRecalibrator.truth_overlap_info,
            ped_file=ped_file,
            benchmark_args=benchmark_args,
            sv_utils_docker=sv_utils_docker,
            samtools_cloud_docker=samtools_cloud_docker
    }

    output {
        File clean_vcf = TrainGqRecalibrator.clean_vcf
        File clean_vcf_index = TrainGqRecalibrator.clean_vcf_index
        File truth_overlap_info = TrainGqRecalibrator.truth_overlap_info
        File output_optimal_overlap_cutoffs = TrainGqRecalibrator.output_optimal_overlap_cutoffs
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

struct CrossValidationVcfs {
    Array[File] train_vcfs
    Array[File] train_vcf_indices
    Array[File] test_vcfs
    Array[File] test_vcf_indices
}

task MakeCrossValidationVcfs {
    input {
        File vcf
        File ped_file
        Array[File] truth_vcfs
        Array[String] vapor_sample_ids
        Int num_splits = 5
        String sv_utils_docker
    }

    String fixed_vcf_name = sub(sub(basename(vcf), ".gz$", ""), ".vcf$", "_fixed.vcf.gz")
    String index_file_name = fixed_vcf_name + ".tbi"

    Int disk_gb = 1000 + round((1 + num_splits) * size(vcf, "GiB") + size(truth_vcfs, "GiB") + size(ped_file, "GiB"))
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

        # create TRUTH_SAMPLES_FILE, with one sample ID per line that has VaPoR/PacBio data (duplicates are okay)
        # The cross-validated VCFs will attempt to distribute these truth samples evenly across batches
        TRUTH_SAMPLES_FILE=~{write_lines(vapor_sample_ids)}
        cat ~{write_lines(truth_vcfs)} | while read TRUTH_VCF; do
            zgrep -m1 ^#[^#] "$TRUTH_VCF" | cut -f10- | tr '\t' '\n'
        done >> "$TRUTH_SAMPLES_FILE"

        sv-utils make-cross-validation-vcfs ~{vcf} \
            --ped-file ~{ped_file} \
            --truth-samples-file "$TRUTH_SAMPLES_FILE" \
            --num_splits ~{num_splits} \
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

    Int disk_gb = 1000 + 2 * round(size(filtered_vcfs, "GiB") + size(filtered_vcf_indices, "GiB"))
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

        bcftools merge --threads $(nproc) -m both -O z -o ~{merged_name} ~{sep=" " filtered_vcfs}

        bcftools index --tbi -o ~{merged_index} ~{merged_name}
    >>>

    output {
        File merged_vcf = merged_name
        File merged_vcf_index = merged_index
    }
}
