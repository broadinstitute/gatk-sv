version development

workflow GetTruthOverlap {
    input {
        Array[File] test_vcfs
        Array[File]? test_vcf_indices
        Array[File] truth_vcfs
        Array[File]? truth_vcf_indices
        Map[String, File]? vapor_files
        Array[File] ped_files
        File? optimal_overlap_cutoffs
        Array[String] get_truth_overlap_args = []
        String module03_docker
        String sv_base_docker
    }

    if(defined(vapor_files)) {
        scatter(vapor_pair in as_pairs(select_first([vapor_files]))) {
            File vapor_file = vapor_pair.right
        }
        Float vapor_size_gib = size(vapor_file, "GiB")
    }


    scatter(test_index in range(length(test_vcfs))) {
        File test_vcf = test_vcfs[test_index]
        if(defined(test_vcf_indices)) {
            File test_vcf_index = select_first([test_vcf_indices])[test_index]
        }
        call GetVcfSize as GetTestVcfSize {
            input:
                vcf=test_vcf,
                vcf_index=test_vcf_index,
                sv_base_docker=sv_base_docker
        }
    }
    call MaxInts as MaxTestVcfNumRecords {
        input:
            ints=GetTestVcfSize.num_records
    }
    Int max_test_records = MaxTestVcfNumRecords.max_int

    scatter(truth_index in range(length(truth_vcfs))) {
        File truth_vcf = truth_vcfs[truth_index]
        if(defined(truth_vcf_indices)) {
            File truth_vcf_index = select_first([truth_vcf_indices])[truth_index]
        }
        call GetVcfSize as GetTruthVcfSize {
            input:
                vcf=truth_vcf,
                vcf_index=truth_vcf_index,
                sv_base_docker=sv_base_docker
        }
    }

    call MaxInts as MaxTruthVcfNumRecords {
        input:
            ints=GetTruthVcfSize.num_records
    }
    Int max_truth_records = MaxTruthVcfNumRecords.max_int


    call GetTruthOverlapTask {
        input:
            test_vcfs=test_vcfs,
            truth_vcfs=truth_vcfs,
            vapor_files=vapor_files,
            vapor_size_gib=select_first([vapor_size_gib, 0.0]),
            ped_files=ped_files,
            optimal_overlap_cutoffs=optimal_overlap_cutoffs,
            get_truth_overlap_args=get_truth_overlap_args,
            module03_docker=module03_docker,
            max_test_records=max_test_records,
            max_truth_records=max_test_records
    }

    output {
        File truth_overlap_info = GetTruthOverlapTask.truth_overlap_info
        File output_optimal_overlap_cutoffs = GetTruthOverlapTask.output_optimal_overlap_cutoffs
    }
}


task MaxInts {
    input {
        Array[Int] ints
    }

    command <<<
        awk 'BEGIN {z=0;} {z=(z>=$0?z:$0);} END {print z;}' "~{write_lines(ints)}"
    >>>

    output {
        Int max_int = read_int(stdout())
    }

    runtime {
        docker: "ubuntu:latest"
        cpu: 1
        preemptible: 3
        max_retries: 0
        memory: "1 GiB"
        disks: "local-disk 10 HDD"
    }
}


task GetVcfSize {
    input {
        File vcf
        File? vcf_index
        String sv_base_docker
    }

    parameter_meta {
        vcf: {
          localization_optional: true
        }
    }

    Int disk_gb = round(10 + size(vcf_index, "GiB"))
    String num_records_file = "num_records.txt"
    String num_samples_file = "num_samples.txt"
    File localized_vcf_index = select_first([vcf_index, vcf + ".tbi"])

    runtime {
        docker: sv_base_docker
        cpu: 1
        preemptible: 3
        max_retries: 0
        memory: "2 GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -eu -o pipefail

        ln -s ~{localized_vcf_index} .
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        bcftools query -l ~{vcf} | wc -w > ~{num_samples_file}
        # get num records from index. If index is wrong kind or somehow missing, build index on the fly.
        {
            bcftools index --nrecords ~{vcf} || {
                bcftools index -f -t -o "$(basename ~{vcf}).tbi" ~{vcf}
                bcftools index --nrecords ~{vcf}
            }
        } > ~{num_records_file}
    >>>

    output {
        Int num_records = read_int(num_records_file)
        Int num_samples = read_int(num_samples_file)
        Int num_entries = num_records * num_samples
    }
}


task GetTruthOverlapTask {
    input {
        Array[File] test_vcfs
        Array[File] truth_vcfs
        Map[String, File]? vapor_files
        Float vapor_size_gib
        Array[File] ped_files
        File? optimal_overlap_cutoffs
        Array[String] get_truth_overlap_args = []
        String module03_docker
        Int? max_test_records
        Int? max_truth_records
    }

    # High disk size for large throughput. A large proportion of run time is loading data from huge VCFs. Disk is cheap.
    Int disk_gb = round(4000 + size(test_vcfs, "GiB") + size(truth_vcfs, "GiB") + size(ped_files, "GiB")
                        + size(optimal_overlap_cutoffs, "GiB") + vapor_size_gib)
    String optimal_overlap_cutoffs_filename = if defined(optimal_overlap_cutoffs)
        then basename(select_first([optimal_overlap_cutoffs]))
        else "optimal_overlap_cutoffs.pickle"

    String truth_overlap_info_filename = "truth_overlap.json"
    String args_str = if length(get_truth_overlap_args) > 0 then sep(" ", get_truth_overlap_args) else ""

    Float mem_baseline = 1.0
    Float mem_scale = "3.5e-6"
    Int mem_gb = if (defined(max_test_records) && defined(max_truth_records))
        then round(mem_baseline + mem_scale * (select_first([max_test_records]) + select_first([max_truth_records])))
        else 3.0

    runtime {
        docker: module03_docker
        cpu: 12
        preemptible: 3
        max_retries: 0
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -eu -o pipefail

        ~{if defined(optimal_overlap_cutoffs) then "cp " + select_first([optimal_overlap_cutoffs]) + " ." else ""}

        # have to declare bash variable because cromwell doesn't localize json file if the arg string is added to it
        VAPOR_JSON=~{if defined(vapor_files) then write_json(select_first([vapor_files])) else ""}
        module03 get_truth_overlap \
            --test-vcf ~{sep=" --test-vcf " test_vcfs} \
            --truth-vcf ~{sep=" --truth-vcf " truth_vcfs} \
            ~{if defined(vapor_files) then "--vapor-json \"$VAPOR_JSON\"" else ""} \
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


