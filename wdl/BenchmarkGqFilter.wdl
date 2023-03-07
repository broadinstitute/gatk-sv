version 1.0

import "Utils.wdl" as Utils
import "ExtractVcfGqFilterProperties.wdl" as ExtractVcfGqFilterProperties


workflow BenchmarkGqFilter {
    input {
        String data_label
        ScoresDataSet original_scores
        File variant_properties_parquet_tar
        Array[ScoresDataSet] comparison_scores
        File benchmark_truth_json
        File? ped_file
        Array[String] benchmark_args = []
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
    }

    if(!defined(original_scores.scores_parquet_tar)) {
        call ExtractVcfScores as ExtractOriginalScores {
            input:
                vcf=original_scores.vcf,
                wanted_properties=[original_scores.score_property],
                sv_utils_docker=sv_utils_docker
        }
    }
    File original_scores_parquet_tar_ = select_first([original_scores.scores_parquet_tar,
                                                     ExtractOriginalScores.scores_file])

    ExtractedScoresDataSet original_scores_ = {
        "label": original_scores.label,
        "scores_parquet_tar": original_scores_parquet_tar_,
        "score_property": original_scores.property,
        "passing_score": original_scores.passing_score
    }

    scatter(comparison_data in comparison_scores) {
        if(!defined(comparison_data.scores_parquet_tar)) {
            call ExtractVcfScores as ExtractComparisonScores {
                input:
                    vcf=comparison_data.vcf,
                    wanted_properties=[comparison_data.score_property],
                    sv_utils_docker=sv_utils_docker
            }
        }
        File comparison_scores_parquet_tar_ = select_first([comparison_data.scores_parquet_tar,
                                                            ExtractComparisonScores.scores_file])
        PickledScoresDataSet comparison_scores_ = {
            "label": comparison_data.label,
            "scores_parquet_tar": comparison_scores_parquet_tar_,
            "score_property": comparison_data.property,
            "passing_score": comparison_data.passing_score
        }
    }

    Array[ExtractedScoresDataSet] scores_data_sets_ = flatten(
        [[original_scores_,], comparison_scores_]
    )

    scatter(make_figure in make_figures) {
        call BenchmarkFilter {
            input:
                data_label=data_label,
                variant_properties_parquet_tar=variant_properties_parquet_tar,
                scores_data_sets=scores_data_sets_,
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
        File original_scores_parquet_tar = original_scores_parquet_tar_
        Array[File] comparison_scores_parquet_tar = comparison_scores_parquet_tar_
    }
}

# for when pickled scores may or may not be available, used to pass data to the workflow
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
        File variant_properties_parquet_tar
        Array[ExtractedScoresDataSet] scores_data_sets
        File benchmark_truth_json
        File? ped_file
        String make_figure
        String benchmark_figure_filename = "quality-benchmark.pdf"
        Array[String] benchmark_args = []
        String sv_utils_docker
        Float mem_gb = 12
    }

    String pedigree_arg = if defined(ped_file)
        then "--pedigree " + ped_file
        else ""

    # some files need to be extracted from tars, so they can double in size. Some processing also
    # uses temp space
    Float extract_file_scale = 4
    Float files_size = size(scores_data_sets, "GiB")
                   + size([variant_properties_parquet_tar, benchmark_truth_json, ped_file], "GiB")
    Int disk_gb = round(50 + files_size * extract_file_scale)
    String scores_data_json = "scores_data.json"

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
        python - <<'CODE'
import json
with open("~{write_json(scores_data_sets)}", 'r') as f_in:
  input_json = json.load(f_in)
output_json={
  "~{data_label}" : {
    "variants_parquet" : "~{variant_properties_parquet_tar}",
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

        gq-recalibrator benchmark-variant-filter \
            --overlap-results ~{benchmark_truth_json} \
            ~{pedigree_arg} \
            --scores-data-json ~{scores_data_json} \
            --figure-save-file ~{benchmark_figure_filename} \
            --make-figure ~{make_figure} \
            ~{sep=' ' benchmark_args}
    >>>

    output {
        File benchmark_figure = benchmark_figure_filename
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

task ExtractVcfScores {
    input {
        File vcf
        Array[String] wanted_properties
        String sv_utils_docker
        Float mem_gb = 8
    }

    Int disk_gb = round(50 + size([vcf], "GiB"))

    # create output filename:
    #   strip .vcf.gz from end (if present)
    #   add wanted properties (separated by underscores)
    #   and add ".pickle.bz2"
    #   (sep doesn't work the way it's supposed to outside of command block, so move part of this logic to command)
    String vcf_basename = sub(sub(basename(vcf), ".gz$", ""), ".vcf$", "")
    String scores_file_name = vcf_basename + ".pq.tar"

    runtime {
        docker: sv_utils_docker
        cpu: 1
        preemptible: 1
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<

    python - <<'CODE'
from pathlib import Path
from sv_utils import genomics_io
from gq_recalibrator import tarred_properties_to_parquet

df = genomics_io.vcf_to_dask(
  "~{vcf}", wanted_properties=["~{sep='", "' wanted_properties}"]
)
tarred_properties_to_parquet.df_to_parquet(
        df=df, output_path=Path("~{scores_file_name}")
)
CODE
    >>>

    output {
        File scores_file = scores_file_name
    }
}