version 1.0

import "Utils.wdl" as Utils
import "ExtractVcfGqFilterProperties.wdl" as ExtractVcfGqFilterProperties
import "ExtractVcfScores.wdl" as ExtractVcfScores


workflow BenchmarkGqFilter {
    input {
        String data_label
        Array[ScoresDataSet] comparison_scores
        File benchmark_truth_json
        File? ped_file
        Array[String] benchmark_args = []
        String samtools_cloud_docker
        String sv_utils_docker
        # Overridable selection of which figures to make, and what order they'll appear in the final
        # PDF. By default make all of them, but you could pass a subset of this list to make fewer
        # (or change the order). Note on runtime: most of these plots take about one hour on a large
        # (4151 samples x 437931 records) data set with original scores, recalibrated,
        # cross-validated, and two additional comparison scores. However, the hardy-weinberg plot
        # takes about 2.5 hours. All plots run independently via cromwell scatter.
        Array[String] make_figures = [
            "precision-recall", "inheritance", "variants-per-sample", "violation-curve",
            "scores-histogram", "hardy-weinberg"
        ]
        Array[String] needed_benchmark_properties = [
            "ac", "svlen", "svtype", "is_autosome"
        ]
    }

    scatter(comparison_data in comparison_scores) {
        if(!defined(comparison_data.scores_parquet_tar)) {
            call ExtractVcfScores.ExtractVcfScores as ExtractComparisonScores {
                input:
                    vcf=comparison_data.vcf,
                    vcf_index=comparison_data.vcf_index,
                    wanted_properties=flatten(
                        [needed_benchmark_properties, [comparison_data.score_property]]
                    ),
                    compute_weights=false,
                    samtools_cloud_docker=samtools_cloud_docker,
                    sv_utils_docker=sv_utils_docker
            }
        }
        File comparison_scores_parquet_tar_ = select_first(
            [comparison_data.scores_parquet_tar, ExtractComparisonScores.properties_parquet_tar]
        )
        ExtractedScoresDataSet scores_data_sets = {
            "label": comparison_data.label,
            "scores_parquet_tar": comparison_scores_parquet_tar_,
            "score_property": comparison_data.score_property,
            "passing_score": comparison_data.passing_score
        }
    }
    Float scores_data_sets_gib = size(comparison_scores_parquet_tar_, "GiB")

    scatter(make_figure in make_figures) {
        call BenchmarkFilter {
            input:
                data_label=data_label,
                scores_data_sets=scores_data_sets,
                scores_data_sets_gib=scores_data_sets_gib,
                benchmark_truth_json=benchmark_truth_json,
                ped_file=ped_file,
                make_figure=make_figure,
                benchmark_figure_filename="benchmark-" + make_figure + ".pdf",
                benchmark_args=benchmark_args,
                sv_utils_docker=sv_utils_docker
        }
    }

    call ConcatBenchmarkPdfs {
        input:
            input_pdfs=BenchmarkFilter.benchmark_figure,
            concat_pdfs_docker=sv_utils_docker
    }

    output {
        File benchmark_figure = ConcatBenchmarkPdfs.benchmark_figure
        Array[File] comparison_scores_parquet_tar = comparison_scores_parquet_tar_
    }
}

# for when scores_parquet_tar may or may not be available, used to pass data to the workflow
struct ScoresDataSet {
    String label
    File vcf
    File vcf_index
    File? scores_parquet_tar
    String score_property
    Int passing_score
}

# for when the scores_parquet_tar is definitely available, only used from within the Workflow
struct ExtractedScoresDataSet {
    String label
    File scores_parquet_tar
    String score_property
    Int passing_score
}

task BenchmarkFilter {
    input {
        String data_label
        Array[ExtractedScoresDataSet] scores_data_sets
        Float scores_data_sets_gib  # size in GiB of data on disk, not in RAM
        File benchmark_truth_json
        File? ped_file
        String make_figure
        String benchmark_figure_filename = "quality-benchmark.pdf"
        Array[String] benchmark_args = []
        String sv_utils_docker
        Float mem_gb_baseline = 2.0
        Float mem_gb_per_cpu = 1.5
        Int num_cpu = 16
    }

    Float mem_gb = mem_gb_baseline + mem_gb_per_cpu * num_cpu
    # some files need to be extracted from tars, so they can double in size. Some processing also
    # uses temp space
    Float extract_file_scale = 4
    Float files_size = scores_data_sets_gib
                     + size([benchmark_truth_json, ped_file], "GiB")
    Int disk_gb = round(50 + files_size * extract_file_scale)
    String scores_data_json = "scores-data.json"
    String performance_report_html = "performance-report.html"

    runtime {
        docker: sv_utils_docker
        cpu: num_cpu
        preemptible: 3
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -euo pipefail

        # transform scores_data_sets into input format expected by benchmark_variant_filter
        python - <<'CODE'
import json
with open("~{write_json(scores_data_sets)}", 'r') as f_in:
  input_json = json.load(f_in)
output_json={
  scores_source["label"] : {
    key: value for key, value in scores_source.items()
    if key != "label"
  }
  for scores_source in input_json
}

with open("~{scores_data_json}", 'w') as f_out:
  json.dump(output_json, f_out, indent=2)
CODE

        # just for debugging:
        echo "~{scores_data_json}:"
        cat ~{scores_data_json}
        echo

        if ~{defined(ped_file)}; then
            PED_ARG="--ped-file ~{ped_file}"
        else
            PED_ARG=""
        fi

        gq-recalibrator benchmark-variant-filter \
            --truth-json "~{benchmark_truth_json}" \
            $PED_ARG \
            --scores-data-json "~{scores_data_json}" \
            --data-label "~{data_label}" \
            --figure-save-file "~{benchmark_figure_filename}" \
            --make-figure "~{make_figure}" \
            --performance-report "~{performance_report_html}" \
            --temp-dir "$(mktemp -d --tmpdir=.)" \
            ~{sep=' ' benchmark_args}
    >>>

    output {
        File benchmark_figure = benchmark_figure_filename
        File performance_report = performance_report_html
    }
}

task ConcatBenchmarkPdfs {
    input {
        Array[File] input_pdfs
        String concat_pdfs_docker
        String benchmark_figure_filename = "quality-benchmark.pdf"
    }

    Float mem_gb = 2.0
    Int disk_gb = round(100 + 2 * size(input_pdfs, "GiB"))

    runtime {
        docker: concat_pdfs_docker
        cpu: 1
        preemptible: 1
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite  -dCompatibilityLevel=1.4 -dPDFSETTINGS=/default -r150 \
            -sOutputFile="~{benchmark_figure_filename}" \
            ~{sep=' ' input_pdfs}
    >>>

    output {
        File benchmark_figure = benchmark_figure_filename
    }
}
