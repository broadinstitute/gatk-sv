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
        String train_vcf_label
        File truth_json
        Array[File] genome_tracks
        Int num_splits = 5
        String sv_utils_docker
        String gatk_docker
        String sv_pipeline_docker
        String sv_base_mini_docker
        String samtools_cloud_docker
        # optional arguments for overriding default tool behaviors
        Array[String] train_args = []
        Array[String] recalibrate_gq_args = []
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

    call TrainGqRecalibrator.TrainGqRecalibrator {
        input:
            train_vcf=train_vcf,
            train_vcf_index=train_vcf_index,
            truth_json=truth_json,
            genome_tracks=genome_tracks,
            train_args=train_args,
            gatk_docker=gatk_docker,
            samtools_cloud_docker=samtools_cloud_docker
    }

    call RecalibrateGq.RecalibrateGq as DirectGqRecalibrator {
        input:
            vcf=train_vcf,
            vcf_index=train_vcf_index,
            genome_tracks=genome_tracks,
            gq_recalibrator_model_file=TrainGqRecalibrator.output_gq_recalibrator_model_file,
            recalibrate_gq_args=recalibrate_gq_args,
            gatk_docker=gatk_docker,
            sv_base_mini_docker=sv_base_mini_docker,
            sv_pipeline_docker=sv_pipeline_docker
    }

    if(!defined(pre_computed_cross_validation_vcfs)) {
        call MakeCrossValidationVcfs {
            input:
                vcf=train_vcf,
                truth_vcfs=[],
                truth_json=truth_json,
                vapor_sample_ids=[],
                num_splits=num_splits,
                sv_utils_docker=sv_utils_docker
        }
    }
    CrossValidationVcfs cross_validation_vcfs_ = select_first([MakeCrossValidationVcfs.cross_validation_vcfs,
                                                               pre_computed_cross_validation_vcfs])

    scatter(fold in range(length(cross_validation_vcfs_.train_vcfs))) {
        call TrainGqRecalibrator.TrainGqRecalibrator as CrossTrainGqRecalibratorFold {
            input:
                train_vcf=cross_validation_vcfs_.train_vcfs[fold],
                train_vcf_index=cross_validation_vcfs_.train_vcf_indices[fold],
                truth_json=truth_json,
                genome_tracks=genome_tracks,
                train_args=train_args,
                gatk_docker=gatk_docker,
                samtools_cloud_docker=samtools_cloud_docker
        }

        call RecalibrateGq.RecalibrateGq as CrossGqRecalibrator {
            input:
                vcf=cross_validation_vcfs_.test_vcfs[fold],
                vcf_index=cross_validation_vcfs_.test_vcf_indices[fold],
                genome_tracks=genome_tracks,
                gq_recalibrator_model_file=CrossTrainGqRecalibratorFold.output_gq_recalibrator_model_file,
                recalibrate_gq_args=recalibrate_gq_args,
                gatk_docker=gatk_docker,
                sv_base_mini_docker=sv_base_mini_docker,
                sv_pipeline_docker=sv_pipeline_docker
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
                vcf=train_vcf,
                vcf_index=train_vcf_index,
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
                vcf=train_vcf,
                vcf_index=train_vcf_index,
                wanted_properties=[old_pipeline_score_property],
                samtools_cloud_docker=samtools_cloud_docker,
                sv_utils_docker=sv_utils_docker
        }
    }
    File pickled_original_scores_ = select_first([pre_computed_pickled_original_scores,
                                                  PickleTrainScores.pickled_properties])

    ScoresDataSet original_scores = {
        "label": "original",
        "vcf": train_vcf,
        "vcf_index": train_vcf_index,
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
            truth_overlap_info=truth_json,
            benchmark_args=benchmark_args,
            sv_utils_docker=sv_utils_docker,
            samtools_cloud_docker=samtools_cloud_docker
    }

    output {
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
        File? ped_file
        Array[File] truth_vcfs
        Array[String] vapor_sample_ids
        File? truth_json
        Int num_splits = 5
        String sv_utils_docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 1000 + round((1 + num_splits) * size(vcf, "GiB") + size(truth_vcfs, "GiB") + size(ped_file, "GiB"))
    RuntimeAttr default_attr = object {
                                   cpu_cores: 1,
                                   mem_gb: 2,
                                   disk_gb: disk_gb,
                                   boot_disk_gb: 10,
                                   preemptible_tries: 3,
                                   max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String fixed_vcf_name = sub(sub(basename(vcf), ".gz$", ""), ".vcf$", "_fixed.vcf.gz")
    String index_file_name = fixed_vcf_name + ".tbi"

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_utils_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    command <<<
        set -euo pipefail

        # create TRUTH_SAMPLES_FILE, with one sample ID per line that has VaPoR/PacBio data (duplicates are okay)
        # The cross-validated VCFs will attempt to distribute these truth samples evenly across batches
        TRUTH_SAMPLES_FILE=~{write_lines(vapor_sample_ids)}
        sed -e 's/[[:space:]]*#.*// ; /^[[:space:]]*$/d' ~{write_lines(truth_vcfs)} \
            | while read TRUTH_VCF; do
                zgrep -m1 ^#[^#] "$TRUTH_VCF" | cut -f10- | tr '\t' '\n'
            done >> "$TRUTH_SAMPLES_FILE"

        { python - <<'CODE'
import json
input_json="~{truth_json}"
if input_json:
    with open(input_json, 'r') as f_in:
        truth_overlap = json.load(f_in)
        for sample_id in truth_overlap.keys():
            print(sample_id)
CODE
} >> "$TRUTH_SAMPLES_FILE"

        sv-utils make-cross-validation-vcfs ~{vcf} \
            ~{"--ped-file " + ped_file} \
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
        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 1000 + 2 * round(size(filtered_vcfs, "GiB") + size(filtered_vcf_indices, "GiB"))
    RuntimeAttr default_attr = object {
                                   cpu_cores: 1,
                                   mem_gb: 4,
                                   disk_gb: disk_gb,
                                   boot_disk_gb: 10,
                                   preemptible_tries: 3,
                                   max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String merged_index = merged_name + ".tbi"

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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
