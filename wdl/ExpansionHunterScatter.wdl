version 1.0

import "Structs.wdl"
import "ExpansionHunter.wdl" as ExpansionHunter

workflow ExpansionHunterScatter {

    input {
        Array[File] bams_or_crams
        Array[File]? bams_or_crams_indexes
        Array[String]? sample_ids
        File reference_fasta
        File? reference_fasta_index
        File variant_catalog
        Int? variant_catalog_batch_size
        String expansion_hunter_docker
        RuntimeAttr? runtime_attr
    }

    String variant_catalog_batch_size_ =
        if defined(variant_catalog_batch_size) then
            select_first([variant_catalog_batch_size])
        else
            10000

    call SplitVariantCatalog {
        input:
            variant_catalog = variant_catalog,
            batch_size = variant_catalog_batch_size_
    }

    scatter (i in range(length(bams_or_crams))) {
        File bam_or_cram_ = bams_or_crams[i]
        Boolean is_bam =
            basename(bam_or_cram_, ".bam") + ".bam" == basename(bam_or_cram_)
        File bam_or_cram_index_ =
            if defined(bams_or_crams_indexes) then
                select_first([bams_or_crams_indexes])[i]
            else
                bam_or_cram_ + if is_bam then ".bai" else ".crai"
        File reference_fasta_index_ = select_first([
            reference_fasta_index, reference_fasta + ".fai"])

        String output_prefix =
            if defined(sample_ids) then
                select_first([sample_ids])[i]
            else
                if is_bam then
                    basename(bam_or_cram_, ".bam")
                else
                    basename(bam_or_cram_, ".cram")

        scatter (j in range(length(SplitVariantCatalog.catalogs_json))) {
            File variant_catalog_split = SplitVariantCatalog.catalogs_json[j]
            call ExpansionHunter.ExpansionHunter as expanionHunter {
                input:
                    bam_or_cram=bam_or_cram_,
                    bam_or_cram_index=bam_or_cram_index_,
                    reference_fasta=reference_fasta,
                    reference_fasta_index=reference_fasta_index_,
                    variant_catalog=variant_catalog_split,
                    output_prefix=output_prefix,
                    expansion_hunter_docker=expansion_hunter_docker
            }
        }
    }

    output {
        Array[Array[File]] jsons = expanionHunter.json
        Array[Array[File]] vcfs = expanionHunter.vcf
        Array[Array[File]] overlapping_reads = expanionHunter.overlapping_reads
        Array[Array[File]] timing = expanionHunter.timing
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
