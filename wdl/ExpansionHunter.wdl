##  This WDL implements workflow for ExpansionHunter.

version 1.0

import "Structs.wdl"

struct FilenamePostfixes {
    String locus
    String motif
    String profile
    String merged_profile
    Int profile_len
}

workflow ExpansionHunter {

    input {
        File bam_or_cram
        File? bam_or_cram_index
        File reference_fasta
        File? reference_fasta_index
        File variant_catalog
        File? output_prefix
        String expansion_hunter_docker
        String sv_base_mini_docker
        Int variant_catalog_batch_size
        RuntimeAttr? runtime_attr
        RuntimeAttr? runtime_override_concat
    }

    Boolean is_bam = basename(bam_or_cram, ".bam") + ".bam" == basename(bam_or_cram)
    File bam_or_cram_index_ =
        if defined(bam_or_cram_index) then
            select_first([bam_or_cram_index])
        else
            bam_or_cram + if is_bam then ".bai" else ".crai"

    File reference_fasta_index_ = select_first([
        reference_fasta_index,
        reference_fasta + ".fai"])

    String output_prefix_ =
        if defined(output_prefix) then
            select_first([output_prefix])
        else
            if is_bam then
                basename(bam_or_cram, ".bam")
            else
                basename(bam_or_cram, ".cram")

    call SplitVariantCatalog {
        input:
            variant_catalog = variant_catalog,
            batch_size = variant_catalog_batch_size
    }

    scatter(i in range(length(SplitVariantCatalog.catalogs_json))) {
        File variant_catalog_split = SplitVariantCatalog.catalogs_json[i]
        call RunExpansionHunter {
            input:
                bam_or_cram = bam_or_cram,
                bam_or_cram_index = bam_or_cram_index_,
                reference_fasta = reference_fasta,
                reference_fasta_index = reference_fasta_index_,
                variant_catalog = variant_catalog_split,
                output_prefix = output_prefix_,
                expansion_hunter_docker = expansion_hunter_docker,
                runtime_attr_override = runtime_attr,
        }
    }

    call ConcatEHOutputs {
        input:
            vcfs_gz = RunExpansionHunter.vcf_gz,
            jsons = RunExpansionHunter.json,
            overlapping_reads = RunExpansionHunter.overlapping_reads,
            timings = RunExpansionHunter.timing,
            output_prefix = output_prefix_,
            expansion_hunter_docker = expansion_hunter_docker
    }

    output {
        File json = ConcatEHOutputs.json
        File vcf_gz = ConcatEHOutputs.vcf_gz
        File overlapping_reads = ConcatEHOutputs.overlapping_reads
        File timing = ConcatEHOutputs.timing
    }
}

task RunExpansionHunter {
    input {
        File bam_or_cram
        File bam_or_cram_index
        File reference_fasta
        File reference_fasta_index
        File variant_catalog
        String output_prefix
        String expansion_hunter_docker
        RuntimeAttr? runtime_attr_override
    }

    output {
        File json = "${output_prefix}.json"
        File vcf_gz = "${output_prefix}.vcf.gz"
        File overlapping_reads = "${output_prefix}_realigned.bam"
        File timing = "${output_prefix}_timing.tsv"
    }

    command <<<
        set -euxo pipefail

        ExpansionHunter \
            --reads ~{bam_or_cram} \
            --reference ~{reference_fasta} \
            --variant-catalog ~{variant_catalog} \
            --output-prefix ~{output_prefix} \
            --cache-mates \
            --record-timing

        bgzip ~{output_prefix}.vcf
    >>>

    RuntimeAttr runtime_attr_str_profile_default = object {
        cpu_cores: 1,
        mem_gb: 4,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1,
        disk_gb: 10 + ceil(size([
            bam_or_cram,
            bam_or_cram_index,
            reference_fasta,
            reference_fasta_index], "GiB"))
    }
    RuntimeAttr runtime_attr = select_first([
        runtime_attr_override,
        runtime_attr_str_profile_default])

    runtime {
        docker: expansion_hunter_docker
        cpu: runtime_attr.cpu_cores
        memory: runtime_attr.mem_gb + " GiB"
        disks: "local-disk " + runtime_attr.disk_gb + " HDD"
        bootDiskSizeGb: runtime_attr.boot_disk_gb
        preemptible: runtime_attr.preemptible_tries
        maxRetries: runtime_attr.max_retries
    }
}

task SplitVariantCatalog {
    input {
        File variant_catalog
        Int batch_size
        String? split_variant_catalog_docker
        String? output_prefix
        RuntimeAttr? runtime_attr_override
    }

    output {
        Array[File] catalogs_json = glob("${output_prefix_}*.json")
    }

    String output_prefix_ =
        if defined(output_prefix) then
            select_first([output_prefix])
        else
            "output_"

    String split_variant_catalog_docker_ =
        if defined(split_variant_catalog_docker) then
            select_first([split_variant_catalog_docker])
        else
            "python:slim-buster"

    command <<<
        set -euxo pipefail

        python <<CODE
        import json
        import sys
        import os
        from pathlib import Path

        filename = "~{variant_catalog}"
        filename_without_ext = Path(filename).stem
        output_prefix = "~{output_prefix_}"
        i = 0
        subset_counter = 0

        def serialize():
            print(f"subset counter: {subset_counter}")
            with open(f"{output_prefix}{filename_without_ext}_{subset_counter}.json", "w") as f_out:
                json.dump(subset_catalogs, f_out, indent=4)

        with open(filename, "r") as f_in:
            catalogs = json.load(f_in)

            subset_catalogs = []
            for catalog in catalogs:
                subset_catalogs.append(catalog)
                i += 1
                if i == ~{batch_size}:
                    serialize()
                    i = 0
                    subset_counter += 1
                    subset_catalogs = []
            serialize()
        CODE
    >>>

    RuntimeAttr runtime_attr_str_profile_default = object {
        cpu_cores: 1,
        mem_gb: 4,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1,
        disk_gb: 10
    }
    RuntimeAttr runtime_attr = select_first([
        runtime_attr_override,
        runtime_attr_str_profile_default])

    runtime {
        docker: split_variant_catalog_docker_
        cpu: runtime_attr.cpu_cores
        memory: runtime_attr.mem_gb + " GiB"
        disks: "local-disk " + runtime_attr.disk_gb + " HDD"
        bootDiskSizeGb: runtime_attr.boot_disk_gb
        preemptible: runtime_attr.preemptible_tries
        maxRetries: runtime_attr.max_retries
    }
}

task ConcatEHOutputs {
    input {
        Array[File] vcfs_gz
        Array[File] jsons
        Array[File] overlapping_reads
        Array[File] timings
        String? output_prefix
        String expansion_hunter_docker
        RuntimeAttr? runtime_attr_override
    }

    output {
        File json = "${output_prefix}.json"
        File vcf_gz = "${output_prefix}.vcf.gz"
        File overlapping_reads = "${output_prefix}.bam"
        File timing = "${output_prefix}_timing.tsv"
    }

    command <<<
        set -euxo pipefail

        jq -s 'reduce .[] as $item ({}; . * $item)' ~{sep=" " jsons} > ~{output_prefix}.json

        VCFS="~{write_lines(vcfs_gz)}"
        bcftools concat --no-version --naive --naive-force --output-type z --file-list ${VCFS} --output "~{output_prefix}.vcf.gz"

        BAMS="~{write_lines(overlapping_reads)}"
        samtools merge ~{output_prefix}.bam -b ${BAMS}

        TIMINGS=(~{sep=" " timings})
        FIRST_TIMING=${TIMINGS[0]}
        head -1 ${FIRST_TIMING} > ~{output_prefix}_timing.tsv
        tail -n +2 -q ~{sep=" " timings} >> ~{output_prefix}_timing.tsv
    >>>

    RuntimeAttr runtime_attr_str_profile_default = object {
        cpu_cores: 1,
        mem_gb: 4,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1,
        disk_gb: 10
    }
    RuntimeAttr runtime_attr = select_first([
        runtime_attr_override,
        runtime_attr_str_profile_default])

    runtime {
        docker: expansion_hunter_docker
        cpu: runtime_attr.cpu_cores
        memory: runtime_attr.mem_gb + " GiB"
        disks: "local-disk " + runtime_attr.disk_gb + " HDD"
        bootDiskSizeGb: runtime_attr.boot_disk_gb
        preemptible: runtime_attr.preemptible_tries
        maxRetries: runtime_attr.max_retries
    }
}
