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
        Boolean? generate_realigned_bam
        Boolean? generate_vcf
        Boolean? seeking_analysis_mode
        Boolean? generate_reviewer_images
        Boolean? include_all_fields
        Int? thread_count
        String expansion_hunter_docker
        String python_docker
        RuntimeAttr? runtime_split_var_catalog
        RuntimeAttr? runtime_eh
        RuntimeAttr? runtime_concat
        RuntimeAttr? runtime_reviewer
    }

    parameter_meta {
        ped_file: "This file is used to extract the sex of the BAM/CRAM files."
        sample_ids: "One ID per sample, in the same order as the files in bams_or_crams. These IDs must match the ID given in the second column (`Individual ID` column) of the given PED file. These IDs will also be used as an output prefix."
    }

    String variant_catalog_batch_size_ = select_first([variant_catalog_batch_size, 1000])

    call SplitVariantCatalog {
        input:
            variant_catalog = variant_catalog_json,
            batch_size = variant_catalog_batch_size_,
            output_prefix = basename(variant_catalog_json, ".json"),
            python_docker = python_docker,
            runtime_override = runtime_split_var_catalog
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

        call ExpansionHunter.ExpansionHunter  {
            input:
                bam_or_cram = bam_or_cram_,
                bam_or_cram_index = bam_or_cram_index_,
                reference_fasta = reference_fasta,
                reference_fasta_index = reference_fasta_index_,
                split_variant_catalogs = SplitVariantCatalog.catalogs_json,
                sample_id = sample_id,
                ped_file = ped_file,
                generate_realigned_bam = generate_realigned_bam,
                generate_vcf = generate_vcf,
                seeking_analysis_mode = seeking_analysis_mode,
                generate_reviewer_images = generate_reviewer_images,
                include_all_fields = include_all_fields,
                thread_count = thread_count,
                expansion_hunter_docker = expansion_hunter_docker,
                python_docker = python_docker,
                runtime_eh = runtime_eh,
                runtime_concat = runtime_concat,
                runtime_reviewer = runtime_reviewer
        }
    }

    output {
        Array[File] variants_tsv_gz = ExpansionHunter.variants_tsv_gz
        Array[File] alleles_tsv_gz = ExpansionHunter.alleles_tsv_gz
        Array[File] vcfs_gz = ExpansionHunter.vcf_gz
        Array[File] realigned_bam = ExpansionHunter.realigned_bam
        Array[File] realigned_bam_index = ExpansionHunter.realigned_bam_index
        Array[Array[File]] jsons_gz = ExpansionHunter.jsons_gz
        Array[File] reviewer_metrics = ExpansionHunter.reviewer_metrics
        Array[File] reviewer_missing_metrics = ExpansionHunter.reviewer_missing_metrics
        Array[File] reviewer_images_gz = ExpansionHunter.reviewer_images_gz
    }
}

task SplitVariantCatalog {
    input {
        File variant_catalog
        Int batch_size
        String output_prefix
        String python_docker
        RuntimeAttr? runtime_override
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

    RuntimeAttr runtime_default = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1,
        disk_gb: 10 + (2 * ceil(size(variant_catalog, "GiB")))
    }
    RuntimeAttr runtime_attr = select_first([runtime_override, runtime_default])

    runtime {
        docker: python_docker
        cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, runtime_default.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    }
}
