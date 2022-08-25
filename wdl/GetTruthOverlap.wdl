version 1.0

import "Utils.wdl" as Utils

workflow GetTruthOverlap {
    input {
        Array[File] test_vcfs
        Array[File] test_vcf_indices
        Array[File] truth_vcfs
        Array[File] truth_vcf_indices
        Array[String]? vapor_sample_ids
        Array[File]? vapor_files
        Array[File] ped_files
        File? optimal_overlap_cutoffs
        Array[String] get_truth_overlap_args = []
        String sv_utils_docker
        String samtools_cloud_docker
    }

    scatter(test_index in range(length(test_vcfs))) {
        call Utils.GetVcfSize as GetTestVcfSize {
            input:
                vcf=test_vcfs[test_index],
                vcf_index=test_vcf_indices[test_index],
                samtools_cloud_docker=samtools_cloud_docker
        }
    }
    call Utils.MaxInts as MaxTestVcfNumRecords {
        input:
            ints=GetTestVcfSize.num_records
    }
    Int max_test_records = MaxTestVcfNumRecords.max_int

    scatter(truth_index in range(length(truth_vcfs))) {
        call Utils.GetVcfSize as GetTruthVcfSize {
            input:
                vcf=truth_vcfs[truth_index],
                vcf_index=truth_vcf_indices[truth_index],
                samtools_cloud_docker=samtools_cloud_docker
        }
    }

    call Utils.MaxInts as MaxTruthVcfNumRecords {
        input:
            ints=GetTruthVcfSize.num_records
    }
    Int max_truth_records = MaxTruthVcfNumRecords.max_int


    call GetTruthOverlapTask {
        input:
            test_vcfs=test_vcfs,
            truth_vcfs=truth_vcfs,
            vapor_sample_ids=vapor_sample_ids,
            vapor_files=vapor_files,
            ped_files=ped_files,
            optimal_overlap_cutoffs=optimal_overlap_cutoffs,
            get_truth_overlap_args=get_truth_overlap_args,
            sv_utils_docker=sv_utils_docker,
            max_test_records=max_test_records,
            max_truth_records=max_test_records
    }

    output {
        File truth_overlap_info = GetTruthOverlapTask.truth_overlap_info
        File output_optimal_overlap_cutoffs = GetTruthOverlapTask.output_optimal_overlap_cutoffs
    }
}


struct VaporDatum {
    String sample_id
    File vapor_file
}


task GetTruthOverlapTask {
    input {
        Array[File] test_vcfs
        Array[File] truth_vcfs
        Array[String]? vapor_sample_ids
        Array[File]? vapor_files
        Array[File] ped_files
        File? optimal_overlap_cutoffs
        Array[String] get_truth_overlap_args = []
        String sv_utils_docker
        Int? max_test_records
        Int? max_truth_records
    }

    # High disk size for large throughput. A large proportion of run time is loading data from huge VCFs. Disk is cheap.
    Int disk_gb = round(100 + size(test_vcfs, "GiB") + size(truth_vcfs, "GiB") + size(ped_files, "GiB")
                        + size(optimal_overlap_cutoffs, "GiB") + size(vapor_files, "GiB"))
    String optimal_overlap_cutoffs_filename = if defined(optimal_overlap_cutoffs)
        then basename(select_first([optimal_overlap_cutoffs]))
        else "optimal_overlap_cutoffs.pickle"

    String truth_overlap_info_filename = "truth_overlap.json"
    String args_str = if length(get_truth_overlap_args) > 0 then "~{sep=' ' get_truth_overlap_args}" else ""

    Float mem_baseline = 1.0
    Float mem_scale = "3.5e-6"
    Float mem_gb = if (defined(max_test_records) && defined(max_truth_records))
        then mem_baseline + mem_scale * (select_first([max_test_records]) + select_first([max_truth_records]))
        else 3.0
    String vapor_json = "vapor_data.json"

    runtime {
        docker: sv_utils_docker
        cpu: 12
        preemptible: 3
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -euo pipefail

        ~{if defined(optimal_overlap_cutoffs) then "cp " + select_first([optimal_overlap_cutoffs]) + " ." else ""}

        # construct a vapor JSON file if vapor data was passed
        ~{if defined(vapor_files) then
          "echo '{' ; paste -d: ~{write_lines(vapor_sample_ids)} ~{write_lines(vapor_files)} | sed -e 's/^/\"/' -e 's/:/\":\"/' -e 's/$/\",/'; echo '}'; > ~{vapor_json}"
          else ""}

        sv-utils get-truth-overlap \
            --test-vcf ~{sep=" --test-vcf " test_vcfs} \
            --truth-vcf ~{sep=" --truth-vcf " truth_vcfs} \
            ~{if defined(vapor_files) then "--vapor-json ~{vapor_json}" else ""} \
            --ped-file ~{sep=" --ped-file " ped_files} \
            --optimal-overlap-cutoffs-file ~{optimal_overlap_cutoffs_filename} \
            --output ~{truth_overlap_info_filename} \
            ~{args_str}
    >>>

    output {
        File output_optimal_overlap_cutoffs = optimal_overlap_cutoffs_filename
        File truth_overlap_info = truth_overlap_info_filename
    }
}


