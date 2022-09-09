version 1.0

import "Utils.wdl" as Utils
import "PickleVcfProperties.wdl" as PickleVcfProperties


workflow BenchmarkGqFilter {
    input {
        String data_label
        ScoresDataSet original_scores
        File? pickled_variant_properties
        Array[ScoresDataSet] comparison_scores
        File truth_overlap_info
        File ped_file
        Array[String] benchmark_args = []
        String sv_utils_docker
        String samtools_cloud_docker
    }

    # extract properties from VCF and pickle them. It allows extraction on a cheaper (lower memory) machine in parallel
    # the pickled properties can be loaded very quickly. Some of them can be saved and (optionally) passed to the WDL
    # for faster production of benchmark figures when prototyping.
    if(!defined(pickled_variant_properties)) {
        call PickleVcfProperties.PickleVcfProperties as PickleVariantData {
            input:
                vcf=original_scores.vcf,
                vcf_index=original_scores.vcf_index,
                wanted_properties=["svtype", "svlen", "gt"],
                samtools_cloud_docker=samtools_cloud_docker,
                sv_utils_docker=sv_utils_docker
        }
    }
    File pickled_variant_properties_ = select_first([pickled_variant_properties, PickleVariantData.pickled_properties])

    # for similar reasons, pickle the actual scores of each ScoresDataSet
    if(!defined(original_scores.pickled_scores_file)) {
        call PickleVcfProperties.PickleVcfProperties as PickleOriginalScores {
            input:
                vcf=original_scores.vcf,
                vcf_index=original_scores.vcf_index,
                wanted_properties=[original_scores.property],
                samtools_cloud_docker=samtools_cloud_docker,
                sv_utils_docker=sv_utils_docker
        }
    }
    File pickled_original_scores_ = select_first([original_scores.pickled_scores_file,
                                                  PickleOriginalScores.pickled_properties])
    PickledScoresDataSet original_scores_ = {
        "label": original_scores.label,
        "scores_file": pickled_original_scores_,
        "score_property": original_scores.property,
        "passing_score": original_scores.passing_score
    }

    scatter(comparison_data in comparison_scores) {
        if(!defined(comparison_data.pickled_scores_file)) {
            call PickleVcfProperties.PickleVcfProperties as PickleComparisonScores {
                input:
                    vcf=comparison_data.vcf,
                    vcf_index=comparison_data.vcf_index,
                    wanted_properties=[comparison_data.property],
                    samtools_cloud_docker=samtools_cloud_docker,
                    sv_utils_docker=sv_utils_docker
            }
        }
        File pickled_comparison_scores_ = select_first([comparison_data.pickled_scores_file,
                                                        PickleComparisonScores.pickled_properties])
        PickledScoresDataSet comparison_scores_ = {
            "label": comparison_data.label,
            "scores_file": pickled_comparison_scores_,
            "score_property": comparison_data.property,
            "passing_score": comparison_data.passing_score
        }
    }

    call Utils.GetVcfSize {
        input:
            vcf=original_scores.vcf,
            vcf_index=original_scores.vcf_index,
            samtools_cloud_docker=samtools_cloud_docker
    }

    Array[PickledScoresDataSet] scores_data_sets_ = flatten([[original_scores_,], comparison_scores_])
    Float pickled_files_size = size(
        flatten([[pickled_variant_properties_, pickled_original_scores_], pickled_comparison_scores_]),
        "GiB"
    )

    call BenchmarkFilter {
        input:
            data_label=data_label,
            variant_properties=pickled_variant_properties_,
            scores_data_sets=scores_data_sets_,
            truth_overlap_info=truth_overlap_info,
            ped_file=ped_file,
            benchmark_args=benchmark_args,
            sv_utils_docker=sv_utils_docker,
            num_entries=GetVcfSize.num_entries,
            pickled_files_size=pickled_files_size
    }

    output {
        File benchmark_figure = BenchmarkFilter.benchmark_figure
        File variant_properties=pickled_variant_properties_
        File pickled_original_scores=pickled_original_scores_
        Array[File] pickled_comparison_scores=pickled_comparison_scores_
    }
}

# for when pickled scores may or may not be available, used to pass data to the workflow
struct ScoresDataSet {
    String label
    File vcf
    File vcf_index
    File? pickled_scores_file
    String property
    Int passing_score
}

# for when the pickled scores file is definitely available, only used from within the Workflow
struct PickledScoresDataSet {
    String label
    File scores_file
    String score_property
    Int passing_score
}

task BenchmarkFilter {
    input {
        String data_label
        File variant_properties
        Array[PickledScoresDataSet] scores_data_sets
        File truth_overlap_info
        File ped_file
        Array[String] benchmark_args = []
        String sv_utils_docker
        Int num_entries
        Float pickled_files_size
    }

    # High disk size for large throughput. A large proportion of run time is loading data from large files. Disk is cheap.
    Int disk_gb = round(1000 + pickled_files_size)

    Float mem_scale_vcf_entries = "2.5e-8"
    Float mem_gb_overhead = 2.0
    Float mem_gb = mem_gb_overhead + mem_scale_vcf_entries * num_entries

    String scores_data_json = "scores_data.json"
    String benchmark_figure_filename = "quality-benchmark.pdf"

    runtime {
        docker: sv_utils_docker
        cpu: 1
        preemptible: 1
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -euo pipefail

        # transform scores_data_sets into input format expected by benchmark_variant_filter
        python <<CODE
        import json
        with open("~{write_json(scores_data_sets)}", 'r') as f_in:
          input_json = json.load(f_in)
        output_json={
          "~{data_label}" : {
            "vcf" : "~{variant_properties}",
            "scores_source": {
              scores_source["label"] : {
                key: value for key, value in scores_source.items()
                if key != "label"
              }
              for scores_source in input_json
            }
          }
        }
        with open("~{scores_data_json}", 'w') as f_out:
          json.dump(output_json, f_out, indent=2)
        CODE

        # just for debugging:
        echo "~{scores_data_json}:"
        cat ~{scores_data_json}
        echo

        sv-utils benchmark-variant-filter \
            --overlap-results ~{truth_overlap_info} \
            --ped-file ~{ped_file} \
            --scores-data-json ~{scores_data_json} \
            --figure-save-file ~{benchmark_figure_filename} \
            ~{sep=' ' benchmark_args}
    >>>

    output {
        File benchmark_figure = benchmark_figure_filename
    }
}
