version 1.0

import "Structs.wdl"
import "ExpansionHunter.wdl" as ExpansionHunter

workflow ExpansionHunterScatter {

    input {
        Array[File] bams_or_crams
        Array[File]? bams_or_crams_indexes
        Array[String] sample_ids
        File? ped_file
        File reference_fasta
        File? reference_fasta_index
        File variant_catalog_json
        Int? variant_catalog_batch_size
        String expansion_hunter_docker
        String python_docker
        RuntimeAttr? runtime_attr
    }

    parameter_meta {
        ped_file: "This file is used to extract the sex of the BAM/CRAM files."
        sample_ids: "One ID per sample, in the same order as the files in bams_or_crams. These IDs must match the ID given in the second column (`Individual ID` column) of the given PED file. These IDs will also be used as an output prefix."
    }

    String variant_catalog_batch_size_ = select_first([variant_catalog_batch_size, 1000])

    call SplitVariantCatalog as svc {
        input:
            variant_catalog = variant_catalog_json,
            batch_size = variant_catalog_batch_size_,
            output_prefix = basename(variant_catalog_json, ".json"),
            python_docker = python_docker
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

        String sample_id = sample_ids[i]

        call ExpansionHunter.ExpansionHunter as expanionHunter {
            input:
                bam_or_cram=bam_or_cram_,
                bam_or_cram_index=bam_or_cram_index_,
                reference_fasta=reference_fasta,
                reference_fasta_index=reference_fasta_index_,
                split_variant_catalogs=svc.catalogs_json,
                sample_id=sample_id,
                ped_file=ped_file,
                expansion_hunter_docker=expansion_hunter_docker,
                python_docker=python_docker
        }
    }

    output {
        Array[File] jsons = expanionHunter.json
        Array[File] vcfs_gz = expanionHunter.vcf_gz
        Array[File] overlapping_reads = expanionHunter.overlapping_reads
        Array[File] timing = expanionHunter.timing
    }
}

task SplitVariantCatalog {
    input {
        File variant_catalog
        Int batch_size
        String output_prefix
        String python_docker
        RuntimeAttr? runtime_attr_override
    }

    output {
        Array[File] catalogs_json = glob("~{output_prefix}*.json")
    }

    command <<<
        set -euxo pipefail

        python <<CODE
        import json
        from pathlib import Path

        filename = "~{variant_catalog}"
        filename_without_ext = Path(filename).stem
        output_prefix = "~{output_prefix}"
        i = 0
        subset_counter = 0


        def serialize():
            with open(f"{output_prefix}{filename_without_ext}_{subset_counter:06}.json", "w") as f_out:
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
            if len(subset_catalogs) > 0:
                serialize()
        CODE
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1,
        disk_gb: 10
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        docker: python_docker
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
