version development

import "GetTruthOverlap.wdl" as GetTruthOverlap

workflow TrainModule03VariantFilter {
    input {
        Array[File] train_vcfs
        Array[File] train_ped_files
        Array[File] train_metrics_files
        Array[File] test_vcfs
        Array[File] test_ped_files
        Array[File] test_metrics_files
        Array[File] test_legacy_scores_files
        Array[File] truth_vcfs
        Map[String, File]? vapor_files
        Array[File] genome_tracts
        File? optimal_overlap_cutoffs
        File? genome_tract_split_json
        File? variant_filter_net
        File? metrics_preprocessing
        String module03_docker
        String old_module03_docker = "us.gcr.io/broad-dsde-methods/tbrookin/sv_utils:8d05ff5"
        String less_old_module03_docker = "us.gcr.io/broad-dsde-methods/tbrookin/sv_utils:258b64d"
    }

    scatter(train_vcf in train_vcfs) {
        call GetGenomeTractOverlaps as TrainGenomeTractOverlaps {
            input:
                vcf=train_vcf,
                genome_tracts=genome_tracts,
                genome_tract_split_json=genome_tract_split_json,
                module03_docker=less_old_module03_docker
        }
    }

    call TrainVariantFilter {
        input:
            train_vcfs=train_vcfs,
            metrics_files=train_metrics_files,
            ped_files=train_ped_files,
            genome_tract_overlaps=TrainGenomeTractOverlaps.genome_tract_overlaps,
            variant_filter_net=variant_filter_net,
            metrics_preprocessing=metrics_preprocessing,
            module03_docker=module03_docker
    }

    scatter(test_vcf in test_vcfs) {
        call GetGenomeTractOverlaps as TestGenomeTractOverlaps {
            input:
                vcf=test_vcf,
                genome_tracts=genome_tracts,
                genome_tract_split_json=genome_tract_split_json,
                module03_docker=less_old_module03_docker
        }
    }

    call GetTruthOverlap.GetTruthOverlap as GetTruthOverlap {
        input:
            test_vcfs=test_vcfs,
            truth_vcfs=truth_vcfs,
            vapor_files=vapor_files,
            ped_files=test_ped_files,
            optimal_overlap_cutoffs=optimal_overlap_cutoffs,
            module03_docker=old_module03_docker
    }

    call CheckVariantFilterPrecisionRecall {
        input:
            truth_overlap_info=GetTruthOverlap.truth_overlap_info,
            optimal_overlap_cutoffs=GetTruthOverlap.output_optimal_overlap_cutoffs,
            metrics_files=test_metrics_files,
            legacy_scores_files=test_legacy_scores_files,
            genome_tract_overlaps=TestGenomeTractOverlaps.genome_tract_overlaps,
            net_save_file=TrainVariantFilter.output_variant_filter_net,
            metrics_preprocessing_save_file=TrainVariantFilter.output_metrics_preprocessing,
            module03_docker=module03_docker
    }

    output {
        Array[File] train_genome_tract_overlaps = TrainGenomeTractOverlaps.genome_tract_overlaps
        Array[File] test_genome_tract_overlaps = TestGenomeTractOverlaps.genome_tract_overlaps
        File truth_overlap_info = GetTruthOverlap.truth_overlap_info
        File output_optimal_overlap_cutoffs = GetTruthOverlap.output_optimal_overlap_cutoffs
        File variant_filter_net = TrainVariantFilter.output_variant_filter_net
        File metrics_preprocessing = TrainVariantFilter.output_metrics_preprocessing
        File precision_recall_figure = CheckVariantFilterPrecisionRecall.precision_recall_figure
    }
}


task GetGenomeTractOverlaps {
    input {
        File vcf
        Array[File] genome_tracts
        File? genome_tract_split_json
        String module03_docker
    }

    Int disk_gb = round(500 + 2 * size(vcf, "GiB") + 2 * size(genome_tracts, "GiB")
                        + size(genome_tract_split_json, "GiB"))
    String genome_tract_overlaps_filename = sub(sub(basename(vcf), ".gz", ""), ".vcf", "_overlaps.tsv")
    String tract_split_arg = if defined(genome_tract_split_json)
        then "--tract-split-json " + select_first([genome_tract_split_json])
        else ""

    runtime {
        docker: module03_docker
        cpu: 12
        preemptible: 3
        max_retries: 1
        memory: "64 GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -eu -o pipefail

        module03 get_genome_tract_overlaps \
            --vcf ~{vcf} \
            --tract ~{sep=" --tract " genome_tracts} \
            ~{tract_split_arg} \
            --output ~{genome_tract_overlaps_filename}
    >>>

    output {
        File genome_tract_overlaps = genome_tract_overlaps_filename
    }
}


task TrainVariantFilter {
    input {
        Array[File] train_vcfs
        Array[File] metrics_files
        Array[File] ped_files
        Array[File] genome_tract_overlaps
        File? variant_filter_net
        File? metrics_preprocessing
        String module03_docker
    }

    Int disk_gb = round(500 + size(train_vcfs, "GiB") + size(metrics_files, "GiB") + size(ped_files, "GiB")
                        + size(genome_tract_overlaps, "GiB"))
    String truth_overlap_info_filename = "truth_overlap.json"
    String net_save_file = if defined(variant_filter_net)
        then basename(select_first([variant_filter_net]))
        else "variant_filter.pickle"
    String metrics_preprocessing_save_file = if defined(metrics_preprocessing)
        then basename(select_first([metrics_preprocessing]))
        else "metrics_preprocessing.pickle"
    runtime {
        docker: module03_docker
        cpu: 8
        preemptible: 3
        max_retries: 0
        memory: "64 GiB"
        disks: "local-disk " + disk_gb + " HDD"
        gpuType: "nvidia-tesla-t4"
        gpuCount: 1
    }

    command <<<
        set -eu -o pipefail

        ~{if defined(variant_filter_net) then "cp " + select_first([variant_filter_net]) + " ." else ""}
        ~{if defined(metrics_preprocessing) then "cp " + select_first([metrics_preprocessing]) + " ." else ""}

        module03 train_variant_filter \
            --vcf ~{sep=" --vcf " train_vcfs} \
            --metrics-file ~{sep=" --metrics-file " metrics_files} \
            --ped-file ~{sep=" --ped-file " ped_files} \
            --tract-properties-file ~{sep=" --tract-properties-file " genome_tract_overlaps} \
            --net-save-file ~{net_save_file} \
            --metrics-preprocessing-save-file ~{metrics_preprocessing_save_file}
    >>>

    output {
        File output_variant_filter_net = net_save_file
        File output_metrics_preprocessing = metrics_preprocessing_save_file
    }
}

task CheckVariantFilterPrecisionRecall {
    input {
        File truth_overlap_info
        File optimal_overlap_cutoffs
        Array[File] metrics_files
        Array[File] legacy_scores_files
        Array[File] genome_tract_overlaps
        File net_save_file
        File metrics_preprocessing_save_file
        String module03_docker
    }

    Int disk_gb = round(
        500 + size(truth_overlap_info, "GiB") + size(optimal_overlap_cutoffs, "GiB") + size(metrics_files, "GiB")
        + size(legacy_scores_files, "GiB") + size(genome_tract_overlaps, "GiB") + size(net_save_file, "GiB")
        + size(metrics_preprocessing_save_file, "GiB")
    )
    String figure_save_file = "precision_recall.pdf"


    runtime {
        docker: module03_docker
        cpu: 12
        preemptible: 3
        max_retries: 1
        memory: "64 GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -eu -o pipefail

        module03 check_variant_filter_precision_recall \
            --overlap_results ~{truth_overlap_info} \
            --optimal-overlap-cutoffs-file ~{optimal_overlap_cutoffs} \
            --metrics-file ~{sep=" --metrics-file " metrics_files} \
            --legacy-scores-file ~{sep=" --legacy-scores-file " legacy_scores_files} \
            --tract-properties-file ~{sep=" --tract-properties-file" genome_tract_overlaps} \
            --net-save-file ~{net_save_file} \
            --metrics-preprocessing-save-file ~{metrics_preprocessing_save_file} \
            --figure-save-file ~{figure_save_file}
    >>>

    output {
        File precision_recall_figure = figure_save_file
    }
}
