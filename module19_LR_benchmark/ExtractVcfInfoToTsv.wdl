version 1.0

import "Structs.wdl"

workflow ExtractVcfInfoToTsv {
    input {
        Array[File] vcf_list
        File extract_script          # extract_vcf_info_to_tsv.py
        File analyze_script          # analyze_vep_by_gene.py
        String output_suffix = ".info.tsv"
        String analysis_output_suffix = ".gene_consequence_counts.tsv"
        File vep_rank_file
        File gene_list_file
        Int vep_consequence_col = 9
        Int vep_gene_col = 10
        String predicted_cols_spec
        String python_docker
        RuntimeAttr? runtime_attr_override
    }

    scatter (vcf in vcf_list) {
        call ExtractInfoTask {
            input:
                vcf            = vcf,
                extract_script = extract_script,
                output_suffix  = output_suffix,
                docker         = python_docker,
                runtime_attr_override = runtime_attr_override
        }
        
        call AnalyzeVepByGeneTask {
            input:
                input_tsv      = ExtractInfoTask.out_tsv,
                analyze_script = analyze_script,
                vep_rank_file  = vep_rank_file,
                gene_list_file = gene_list_file,
                vep_consequence_col = vep_consequence_col,
                vep_gene_col   = vep_gene_col,
                predicted_cols = predicted_cols_spec,
                output_suffix  = analysis_output_suffix,
                docker         = python_docker,
                runtime_attr_override = runtime_attr_override
        }
    }

    call ConcatenateGeneCounts {
        input:
            count_files = AnalyzeVepByGeneTask.out_counts,
            python_docker = python_docker,
            runtime_attr_override = runtime_attr_override
    }

    output {
        Array[File] tsv_list = ExtractInfoTask.out_tsv
        Array[File] gene_consequence_counts = AnalyzeVepByGeneTask.out_counts
        File concatenated_gene_counts = ConcatenateGeneCounts.out_merged
    }
}

task ExtractInfoTask {
    input {
        File   vcf
        File   extract_script
        String output_suffix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size   = size(vcf, "GB")
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
        mem_gb:            4,
        disk_gb:           ceil(base_disk_gb + input_size * 2.0),
        cpu_cores:         1,
        preemptible_tries: 2,
        max_retries:       1,
        boot_disk_gb:      10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    String prefix = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail

        python3 ~{extract_script} \
            --input-vcf  ~{vcf} \
            --output-txt ~{prefix}~{output_suffix}
    >>>

    output {
        File out_tsv = "~{prefix}~{output_suffix}"
    }

    runtime {
        docker:         docker
        memory:         select_first([runtime_override.mem_gb,            runtime_default.mem_gb])           + " GB"
        disks:          "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb])   + " HDD"
        cpu:            select_first([runtime_override.cpu_cores,          runtime_default.cpu_cores])
        preemptible:    select_first([runtime_override.preemptible_tries,  runtime_default.preemptible_tries])
        maxRetries:     select_first([runtime_override.max_retries,        runtime_default.max_retries])
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb,       runtime_default.boot_disk_gb])
    }
}

task AnalyzeVepByGeneTask {
    input {
        File   input_tsv
        File   analyze_script
        File   vep_rank_file
        File   gene_list_file
        Int    vep_consequence_col
        Int    vep_gene_col
        String predicted_cols
        String output_suffix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size   = size(input_tsv, "GB")
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
        mem_gb:            4,
        disk_gb:           ceil(base_disk_gb + input_size * 2.0),
        cpu_cores:         1,
        preemptible_tries: 2,
        max_retries:       1,
        boot_disk_gb:      10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    String prefix = basename(input_tsv, ".info.tsv")

    command <<<
        set -euo pipefail

        python3 ~{analyze_script} \
            --input ~{input_tsv} \
            --vep-rank ~{vep_rank_file} \
            --gene-list ~{gene_list_file} \
            --vep-consequence-col ~{vep_consequence_col} \
            --vep-gene-col ~{vep_gene_col} \
            --predicted-cols "~{predicted_cols}" \
            --output ~{prefix}~{output_suffix}
    >>>

    output {
        File out_counts = "~{prefix}~{output_suffix}"
    }

    runtime {
        docker:         docker
        memory:         select_first([runtime_override.mem_gb,            runtime_default.mem_gb])           + " GB"
        disks:          "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb])   + " HDD"
        cpu:            select_first([runtime_override.cpu_cores,          runtime_default.cpu_cores])
        preemptible:    select_first([runtime_override.preemptible_tries,  runtime_default.preemptible_tries])
        maxRetries:     select_first([runtime_override.max_retries,        runtime_default.max_retries])
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb,       runtime_default.boot_disk_gb])
    }
}

task ConcatenateGeneCounts {
    input {
        Array[File] count_files
        String python_docker
        RuntimeAttr? runtime_attr_override
    }

    Float total_size  = size(count_files, "GB")
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
        mem_gb:            2,
        disk_gb:           ceil(base_disk_gb + total_size * 2.0),
        cpu_cores:         1,
        preemptible_tries: 2,
        max_retries:       1,
        boot_disk_gb:      10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    command <<<
        set -euo pipefail

        # Get header from first file
        head -1 ~{count_files[0]} > merged_gene_counts.tsv

        # Append data rows from all files (skipping headers)
        for file in ~{sep=' ' count_files}; do
            tail -n +2 "$file" >> merged_gene_counts.tsv
        done

        # Remove duplicate rows and sort by gene name
        (head -1 merged_gene_counts.tsv; tail -n +2 merged_gene_counts.tsv | sort -u) > temp.tsv
        mv temp.tsv merged_gene_counts.tsv
    >>>

    output {
        File out_merged = "merged_gene_counts.tsv"
    }

    runtime {
        docker:         python_docker
        memory:         select_first([runtime_override.mem_gb,            runtime_default.mem_gb])           + " GB"
        disks:          "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb])   + " HDD"
        cpu:            select_first([runtime_override.cpu_cores,          runtime_default.cpu_cores])
        preemptible:    select_first([runtime_override.preemptible_tries,  runtime_default.preemptible_tries])
        maxRetries:     select_first([runtime_override.max_retries,        runtime_default.max_retries])
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb,       runtime_default.boot_disk_gb])
    }
}
