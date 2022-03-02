version development

import "GetTruthOverlap.wdl" as GetTruthOverlap

workflow BenchmarkGqFilter {
    input {
        File train_vcf
        String? train_vcf_label
        File filtered_vcf
        File cross_validated_filtered_vcf
        File truth_overlap_info
        File ped_file
        File? pickled_variant_properties
        File? pickled_original_scores
        File? pickled_filtered_scores
        File? pickled_cross_validated_scores
        Int new_pipeline_passing_score = 0
        Int old_pipeline_passing_score = 3
        Array[String] benchmark_args = []
        String module03_docker
        String sv_base_docker
    }

    call GetTruthOverlap.GetVcfSize as GetVcfSize {
        input:
            vcf=train_vcf,
            sv_base_docker=sv_base_docker
    }

    # extract properties from VCF and pickle them. It allows extraction on a cheaper (lower memory) machine in parallel
    # the pickled properties can be loaded very quickly. Some of them can be saved and (optionally) passed to the WDL
    # for faster production of benchmark figures when prototyping.
    if(!defined(pickled_variant_properties)) {
        call PickleVcfProperties as PickleVariantData {
            input:
                vcf=train_vcf,
                wanted_properties=["svtype", "svlen", "gt"],
                module03_docker=module03_docker,
                num_entries=GetVcfSize.num_entries
        }
    }
    File pickled_variant_properties_ = select_first([pickled_variant_properties, PickleVariantData.pickled_properties])

    if(!defined(pickled_original_scores)) {
        call PickleVcfProperties as PickleTrainScores {
            input:
                vcf=train_vcf,
                wanted_properties=["gq"],
                module03_docker=module03_docker,
                num_entries=GetVcfSize.num_entries
        }
    }
    File pickled_original_scores_ = select_first([pickled_original_scores, PickleTrainScores.pickled_properties])

    if(!defined(pickled_filtered_scores)) {
        call PickleVcfProperties as PickleFilteredScores {
            input:
                vcf=filtered_vcf,
                wanted_properties=["gq"],
                module03_docker=module03_docker,
                num_entries=GetVcfSize.num_entries
        }
    }
    File pickled_filtered_scores_ = select_first([pickled_filtered_scores, PickleFilteredScores.pickled_properties])

    if(!defined(pickled_cross_validated_scores)) {
        call PickleVcfProperties as PickleCrossValidatedScores {
            input:
                vcf=cross_validated_filtered_vcf,
                wanted_properties=["gq"],
                module03_docker=module03_docker,
                num_entries=GetVcfSize.num_entries
        }
    }
    File pickled_cross_validated_scores_ = select_first([pickled_cross_validated_scores,
                                                        PickleCrossValidatedScores.pickled_properties])
    call BenchmarkFilter {
        input:
            variant_properties=pickled_variant_properties_,
            original_scores=pickled_original_scores_,
            filtered_scores=pickled_filtered_scores_,
            cross_validated_scores=pickled_cross_validated_scores_,
            truth_overlap_info=truth_overlap_info,
            ped_file=ped_file,
            benchmark_args=benchmark_args,
            old_pipeline_passing_score=old_pipeline_passing_score,
            new_pipeline_passing_score=new_pipeline_passing_score,
            module03_docker=module03_docker,
            num_entries=GetVcfSize.num_entries
    }

    output {
        File benchmark_figure = BenchmarkFilter.benchmark_figure
        File variant_properties=pickled_variant_properties_
        File original_scores=pickled_original_scores_
        File filtered_scores=pickled_filtered_scores_
        File cross_validated_scores=pickled_cross_validated_scores_
    }
}

task PickleVcfProperties {
    input {
        File vcf
        Array[String] wanted_properties
        String module03_docker
        Int? num_entries
    }

    # High disk size for large throughput. A large proportion of run time is loading data from huge VCFs. Disk is cheap.
    Int disk_gb = round(
        1000 + size([vcf], "GiB")
    )
    Float mem_scale_vcf_size = 10.0
    Float mem_scale_vcf_entries = if length(wanted_properties) == 1 then "5e-8" else "1.4e-8"
    Float mem_gb_overhead = 2.0
    Float mem_gb = mem_gb_overhead + if defined(num_entries)
        then mem_scale_vcf_entries * select_first([num_entries])
        else mem_scale_vcf_size * size(vcf, "GiB")

    # create output filename:
    #   strip .vcf.gz from end (if present)
    #   add wanted properties (separated by underscores)
    #   and add ".pickle.bz2"
    String pickled_properties_name = sub(sub(basename(vcf), ".gz$", ""), ".vcf$", "")
                                   + "_" + sep("_", wanted_properties)
                                   + ".pickle.bz2"
#    This doesn't work outside command block: it adds an extra bracket, and seemingly quotes in the command block!
#       e.g. pickled_properties_name = hgdp_and_hgsv.cleaned_fixed_cross_validated_["gq"].pickle.bz2 inside
#            pickled_properties_name = hgdp_and_hgsv.cleaned_fixed_cross_validated_[gq].pickle.bz2 outside
#    String pickled_properties_name = sub(sub(basename(vcf), ".gz$", ""), ".vcf$", "")
#                                   + "_" + "~{sep='_' wanted_properties}"
#                                   + ".pickle.bz2"

    runtime {
        docker: module03_docker
        cpu: 1
        preemptible: 1
        max_retries: 0
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<

    python - <<'____EoF'
from sv_utils import genomics_io

variants = genomics_io.vcf_to_pandas(
  "~{vcf}", wanted_properties=["~{sep='", "' wanted_properties}"], drop_trivial_multi_index=True
)

variants.to_pickle("~{pickled_properties_name}")
____EoF
    >>>

    output {
        File pickled_properties = pickled_properties_name
    }
}


struct ScoreDataSet {
  File scores_file
  String label
}

task BenchmarkFilter {
    input {
        File variant_properties
        String? train_vcf_label
        File original_scores
        File filtered_scores
        File cross_validated_scores
        File truth_overlap_info
        File ped_file
        Array[String] benchmark_args = []
        Int new_pipeline_passing_score
        Int old_pipeline_passing_score
        String module03_docker
        Int? num_entries
    }

    # High disk size for large throughput. A large proportion of run time is loading data from huge VCFs. Disk is cheap.
    Int disk_gb = round(
        1000 + size([variant_properties, original_scores, filtered_scores, cross_validated_scores, ped_file,
                     truth_overlap_info], "GiB")
    )

    Boolean is_vcf = sub(sub(basename(variant_properties), ".gz$", ""), ".vcf$", "") != basename(variant_properties)
    Float mem_scale_vcf_size = if is_vcf then 75.0 else 1000.0
    Float mem_scale_vcf_entries = "2.5e-8"
    Float mem_gb_overhead = 2.0
    Float mem_gb = mem_gb_overhead + if defined(num_entries)
        then mem_scale_vcf_entries * select_first([num_entries])
        else mem_scale_vcf_size * size(variant_properties, "GiB")

    String data_label = select_first([train_vcf_label, "cohort"])
    String scores_data_json = "scores_data.json"
    String benchmark_figure_filename = "quality-benchmark.pdf"
    String args_str = if length(benchmark_args) > 0 then sep(" ", benchmark_args) else ""

    runtime {
        docker: module03_docker
        cpu: 1
        preemptible: 1
        max_retries: 0
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -eu -o pipefail

        cat << 'END_SCORES_DATA' > ~{scores_data_json}
{
    "~{data_label}": {
        "vcf": "~{variant_properties}",
        "scores_source": {
            "original": {
                "scores_file": "~{original_scores}",
                "score_property": "gq",
                "passing_score": ~{old_pipeline_passing_score}
            },
            "recalibrated": {
                "scores_file": "~{filtered_scores}",
                "score_property": "gq",
                "passing_score": ~{new_pipeline_passing_score}
            },
            "cross-validated": {
                "scores_file": "~{cross_validated_scores}",
                "score_property": "gq",
                "passing_score": ~{new_pipeline_passing_score}
            }
        }
    }
}
END_SCORES_DATA
        # just for debugging:
        echo "~{scores_data_json}:"
        cat ~{scores_data_json}

        module03 benchmark_variant_filter \
            --overlap-results ~{truth_overlap_info} \
            --ped-file ~{ped_file} \
            --scores-data-json ~{scores_data_json} \
            --figure-save-file ~{benchmark_figure_filename} \
            ~{args_str}
    >>>

    output {
        File benchmark_figure = benchmark_figure_filename
    }
}
