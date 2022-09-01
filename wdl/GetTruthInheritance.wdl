version 1.0

workflow GetTruthInheritance {
    # Get truth JSON using trios from ped files to find "good" variants (that obey mendelian
    # inheritance) or "bad" variants that are de-novo. All members of a good / bad trio are marked
    # accordingly, because e.g. it's not clear if a de-novo is a false positive in the proband or
    # a false negative in one of the parents.
    input {
        Array[File] vcfs
        Array[File]? vcf_indices
        Array[File] ped_files
        String truth_json_file_name
        Array[String] get_truth_inheritance_args = []
        String sv_utils_docker
    }

    call InheritanceToTruthJson {
        input:
            vcfs=vcfs,
            vcf_indices=vcf_indices,
            ped_files=ped_files,
            truth_json_file_name=truth_json_file_name,
            get_truth_inheritance_args=get_truth_inheritance_args,
            sv_utils_docker=sv_utils_docker
    }

    output {
        File inheritance_truth_json = InheritanceToTruthJson.inheritance_truth_json
    }
}


task InheritanceToTruthJson {
    input {
        Array[File] vcfs
        Array[File]? vcf_indices
        Array[File] ped_files
        String truth_json_file_name
        Array[String] get_truth_inheritance_args = []
        String sv_utils_docker
        Float mem_gb = 2.0
    }

    # Double input vcf size in case output json is somehow huge (it shouldn't be, but this should be safe)
    Int disk_gb = round(
        10 + 2 * size(vcfs, "GiB") + size(select_first([vcf_indices]), "GiB") + size(ped_files, "GiB")
    )

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


        sv-utils inheritance-to-truth-json \
            --vcf ~{sep=" --vcf " vcfs} \
            --ped-file ~{sep=" --ped-file " ped_files} \
            --output ~{truth_json_file_name} \
            ~{sep=' ' get_truth_inheritance_args}
    >>>

    output {
        File inheritance_truth_json = truth_json_file_name
    }
}
