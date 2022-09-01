version 1.0

workflow CombineTruthJsons {
    # Get truth JSON using trios from ped files to find "good" variants (that obey mendelian
    # inheritance) or "bad" variants that are de-novo. All members of a good / bad trio are marked
    # accordingly, because e.g. it's not clear if a de-novo is a false positive in the proband or
    # a false negative in one of the parents.
    input {
        Array[File] input_truth_jsons
        String combined_truth_json_file_name
        # combine_strategy must be one of: "prefer-first", "prefer-last", or "omit-conflicting"
        String combine_strategy = "omit-conflicting"
        Array[String] combine_truth_jsons_args = []
        String sv_utils_docker
    }

    call CombineTruthJsonsTask {
        input:
            input_truth_jsons=input_truth_jsons,
            combined_truth_json_file_name=combined_truth_json_file_name,
            combine_strategy=combine_strategy,
            combine_truth_jsons_args=combine_truth_jsons_args,
            sv_utils_docker=sv_utils_docker
    }

    output {
        File combined_truth_json = CombineTruthJsonsTask.combined_truth_json
    }
}


task CombineTruthJsonsTask {
    input {
        Array[File] input_truth_jsons
        String combined_truth_json_file_name
        # combine_strategy must be one of: "prefer-first", "prefer-last", or "omit-conflicting"
        String combine_strategy = "omit-conflicting"
        Array[String] combine_truth_jsons_args = []
        String sv_utils_docker
        Float mem_gb = 8.0
    }

    Int disk_gb = round(50 + 2 * size(input_truth_jsons, "GiB"))

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


        sv-utils combine-truth-jsons \
            --input-json ~{sep=" --input-json " input_truth_jsons} \
            --output-json ~{combined_truth_json_file_name} \
            --combine-strategy ~{combine_strategy} \
            ~{sep=' ' combine_truth_jsons_args}
    >>>

    output {
        File combined_truth_json = combined_truth_json_file_name
    }
}
