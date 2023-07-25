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
        Array[File] split_variant_catalogs
        String sample_id
        Boolean? generate_realigned_bam
        Boolean? generate_vcf
        Boolean? seeking_analysis_mode
        Boolean? generate_reviewer_images
        Boolean? include_all_fields
        Int? thread_count
        File? ped_file
        String expansion_hunter_docker
        String python_docker
        RuntimeAttr? runtime_eh
        RuntimeAttr? runtime_concat
        RuntimeAttr? runtime_reviewer
    }

    parameter_meta {
        ped_file: "This file is used to extract the sex of the bam_or_cram file."
        sample_id: "The ped_file needs to be provided as well to determine sample sex. The ID must match the sample ID given in the second column (`Individual ID` column) of the given PED file. This ID will also be used as an output prefix."
    }

    Boolean generate_realigned_bam_ = select_first([generate_realigned_bam, false])
    Boolean generate_vcf_ = select_first([generate_vcf, false])
    Boolean generate_reviewer_images_ = select_first([generate_reviewer_images, false])
    Boolean include_all_fields_ = select_first([include_all_fields, false])

    Boolean is_bam = basename(bam_or_cram, ".bam") + ".bam" == basename(bam_or_cram)
    File bam_or_cram_index_ =
        if defined(bam_or_cram_index) then
            select_first([bam_or_cram_index])
        else
            bam_or_cram + if is_bam then ".bai" else ".crai"

    File reference_fasta_index_ = select_first([
        reference_fasta_index,
        reference_fasta + ".fai"])

    Int thread_count_ = select_first([thread_count, 2])
    String analysis_mode =
        if select_first([seeking_analysis_mode, true]) then
            "seeking"
        else
            "streaming"

    scatter (i in range(length(split_variant_catalogs))) {
        call RunExpansionHunter {
            input:
                bam_or_cram = bam_or_cram,
                bam_or_cram_index = bam_or_cram_index_,
                reference_fasta = reference_fasta,
                reference_fasta_index = reference_fasta_index_,
                variant_catalog = split_variant_catalogs[i],
                sample_id = sample_id,
                generate_realigned_bam = generate_realigned_bam_,
                generate_vcf = generate_vcf_,
                include_all_fields = include_all_fields_,
                analysis_mode = analysis_mode,
                thread_count = thread_count_,
                ped_file = ped_file,
                expansion_hunter_docker = expansion_hunter_docker,
                runtime_override = runtime_eh
        }

        if (generate_reviewer_images_ && generate_realigned_bam_ && generate_vcf_) {
            call RunReviewer {
                input:
                    sample_id = sample_id,
                    realigned_bam = RunExpansionHunter.realigned_bam,
                    realigned_bam_index = RunExpansionHunter.realigned_bam_index,
                    vcf_gz = RunExpansionHunter.vcf_gz,
                    reference_fasta = reference_fasta,
                    reference_fasta_index = reference_fasta_index_,
                    variant_catalog_json = split_variant_catalogs[i],
                    expansion_hunter_docker = expansion_hunter_docker,
                    runtime_override = runtime_reviewer
            }
        }
    }

    call ConcatEHOutputs {
        input:
            vcfs_gz = RunExpansionHunter.vcf_gz,
            variants_tsvs = RunExpansionHunter.variants_tsv,
            alleles_tsvs = RunExpansionHunter.alleles_tsv,
            realigned_bams = RunExpansionHunter.realigned_bam,
            realigned_bams_index = RunExpansionHunter.realigned_bam_index,
            generate_realigned_bam = generate_realigned_bam_,
            generate_reviewer_images = generate_reviewer_images_,
            generate_vcf = generate_vcf_,
            sample_id = sample_id,
            reviewer_images_gzs = RunReviewer.reviewer_images_gz,
            sample_metrics = RunReviewer.sample_metrics,
            output_prefix = sample_id,
            expansion_hunter_docker = expansion_hunter_docker,
            runtime_override = runtime_concat
    }

    output {
        File variants_tsv_gz = ConcatEHOutputs.variants_tsv_gz
        File alleles_tsv_gz = ConcatEHOutputs.alleles_tsv_gz
        File vcf_gz = ConcatEHOutputs.vcf_gz
        File realigned_bam = ConcatEHOutputs.realigned_bam
        File realigned_bam_index = ConcatEHOutputs.realigned_bam_index
        Array[File] jsons_gz = RunExpansionHunter.json_gz
        File reviewer_metrics = ConcatEHOutputs.reviewer_metrics
        File reviewer_missing_metrics = ConcatEHOutputs.reviewer_missing_metrics
        File reviewer_images_gz = ConcatEHOutputs.reviewer_images_gz
    }
}

task RunExpansionHunter {
    input {
        File bam_or_cram
        File bam_or_cram_index
        File reference_fasta
        File reference_fasta_index
        File variant_catalog
        String sample_id
        Boolean generate_realigned_bam
        Boolean generate_vcf
        Boolean include_all_fields
        String analysis_mode
        Int thread_count
        File? ped_file
        String expansion_hunter_docker
        RuntimeAttr? runtime_override
    }

    output {
        File variants_tsv = "${sample_id}_variants.tsv"
        File alleles_tsv = "${sample_id}_alleles.tsv"
        File vcf_gz = "${sample_id}.vcf.gz"
        File json_gz = "${sample_id}.json.gz"
        File realigned_bam = "${sample_id}_realigned.bam"
        File realigned_bam_index = "${sample_id}_realigned.bam.bai"
    }

    command <<<
        set -euxo pipefail

        BAM_OR_CRAM_DIR="$(dirname "~{bam_or_cram}")"
        BAM_OR_CRAM_INDEX_FILENAME="$(basename "~{bam_or_cram_index}")"
        DEST="$BAM_OR_CRAM_DIR/$BAM_OR_CRAM_INDEX_FILENAME"
        if [ $DEST != ~{bam_or_cram_index} ]; then
            mv ~{bam_or_cram_index} $DEST
        fi

        REF="$(basename "~{reference_fasta}")"
        ln -s ~{reference_fasta} $REF
        ln -s ~{reference_fasta_index} $REF.fai

        sex=""
        if ~{defined(ped_file)}; then
            sex=$(awk -F '\t' '{if ($2 == "~{sample_id}") {if ($5 == "1") {print "--sex male"; exit 0} else if ($5 == "2") {print "--sex female"; exit 0}}}' < "~{if defined(ped_file) then select_first([ped_file]) else "none"}" )
            if [ "$sex" = "" ]; then
                echo "The Sex of the sample defined in the PED file is other than male or female. ExpansionHunter only supports male or female samples."
                exit 1
            fi
        fi

        ExpansionHunter \
            --reads ~{bam_or_cram} \
            --reference $REF \
            --variant-catalog ~{variant_catalog} \
            --output-prefix ~{sample_id} \
            --analysis-mode ~{analysis_mode} \
            --threads ~{thread_count} \
            $sex

        if [ ~{generate_realigned_bam} = false ]; then
            rm ~{sample_id}_realigned.bam
            touch ~{sample_id}_realigned.bam
            touch ~{sample_id}_realigned.bam.bai
        else
            mv ~{sample_id}_realigned.bam ~{sample_id}_realigned_unsorted.bam
            samtools sort ~{sample_id}_realigned_unsorted.bam -o ~{sample_id}_realigned.bam
            samtools index ~{sample_id}_realigned.bam
        fi

        if ~{generate_vcf}; then
            bgzip ~{sample_id}.vcf
        else
            rm ~{sample_id}.vcf
            touch ~{sample_id}.vcf.gz
        fi

        python /opt/str/combine_expansion_hunter_json_to_tsv.py \
            ~{if include_all_fields then "--include-all-fields " else ""} \
            -o ~{sample_id} ~{sample_id}.json
        mv ~{sample_id}.*_json_files_alleles.tsv ~{sample_id}_alleles.tsv
        mv ~{sample_id}.*_json_files_variants.tsv ~{sample_id}_variants.tsv

        gzip ~{sample_id}.json
    >>>

    RuntimeAttr runtime_default = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1,
        disk_gb: 10 + (
            2 * ceil(size([
                bam_or_cram,
                bam_or_cram_index,
                reference_fasta,
                reference_fasta_index], "GiB")))
    }
    RuntimeAttr runtime_attr = select_first([runtime_override, runtime_default])

    runtime {
        docker: expansion_hunter_docker
        cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, runtime_default.disk_gb])  + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    }
}

task ConcatEHOutputs {
    input {
        Array[File] vcfs_gz
        Array[File] variants_tsvs
        Array[File] alleles_tsvs
        Array[File] realigned_bams
        Array[File] realigned_bams_index
        Array[File?] reviewer_images_gzs
        Array[File?] sample_metrics
        String sample_id
        Boolean generate_realigned_bam
        Boolean generate_reviewer_images
        Boolean generate_vcf
        String? output_prefix
        String expansion_hunter_docker
        RuntimeAttr? runtime_override
    }

    output {
        File variants_tsv_gz = "${output_prefix}_variants.tsv.gz"
        File alleles_tsv_gz = "${output_prefix}_alleles.tsv.gz"
        File vcf_gz = "${output_prefix}.vcf.gz"
        File realigned_bam = "${output_prefix}.bam"
        File realigned_bam_index = "${output_prefix}.bam.bai"
        File reviewer_metrics = "${output_prefix}_metrics.csv"
        File reviewer_missing_metrics = "${output_prefix}_missing_metrics.csv"
        File reviewer_images_gz = "${output_prefix}_reviewer_images.tar.gz"
    }

    Array[File] sample_metrics_ = select_all(sample_metrics)
    Array[File] reviewer_images_gz_ = select_all(reviewer_images_gzs)

    command <<<
        set -euxo pipefail

        if ~{generate_vcf}; then
            VCFS="~{write_lines(vcfs_gz)}"
            bcftools concat --no-version --naive-force --output-type z --file-list ${VCFS} --output "~{output_prefix}.vcf.gz"
        else
            touch ~{output_prefix}.vcf.gz
        fi

        if ~{generate_realigned_bam}; then
            BAMS="~{write_lines(realigned_bams)}"
            samtools merge ~{output_prefix}_unsorted.bam -b ${BAMS}
            samtools sort ~{output_prefix}_unsorted.bam -o ~{output_prefix}.bam
            samtools index ~{output_prefix}.bam
        else
            touch ~{output_prefix}.bam
            touch ~{output_prefix}.bam.bai
        fi

        function merge_tsv {
            INPUTS=$1
            OUTPUT_FILENAME=$2

            FIRST_TSV=$(head -n 1 $INPUTS)
            head -1 $FIRST_TSV > $OUTPUT_FILENAME
            while IFS= read -r line; do
                awk FNR!=1 $line >> $OUTPUT_FILENAME
            done < $INPUTS
        }

        merge_tsv "~{write_lines(alleles_tsvs)}" "~{output_prefix}_alleles.tsv"
        merge_tsv "~{write_lines(variants_tsvs)}" "~{output_prefix}_variants.tsv"

        gzip "~{output_prefix}_alleles.tsv"
        gzip "~{output_prefix}_variants.tsv"

        if ~{generate_reviewer_images}; then
            # This will output two files: prefix_metrics.csv & prefix_missing_metrics.csv
            python /opt/str/merge_csv_files.py \
                --input-filename ~{write_lines(sample_metrics_)} \
                --output-prefix ~{output_prefix}

            # Combine multiple archives into one archive,
            # by first unzipping all to a common directory,
            # the archiving all the contents of the directory into a single archive.
            TEMP_DIR_NAME="~{output_prefix}_reviewer_images"
            mkdir $TEMP_DIR_NAME
            while read -r archive_filename; do
                tar -xzf "$archive_filename" -C $TEMP_DIR_NAME
            done < ~{write_lines(reviewer_images_gz_)}
            tar -czf ~{output_prefix}_reviewer_images.tar.gz -C $TEMP_DIR_NAME .
        else
            touch ~{output_prefix}_metrics.csv
            touch ~{output_prefix}_missing_metrics.csv
            touch ~{output_prefix}_reviewer_images.tar.gz
        fi
    >>>

    RuntimeAttr runtime_default = object {
        cpu_cores: 1,
        mem_gb: 4,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1,
        disk_gb: 10 +
            (2 * ceil(
                size(vcfs_gz, "GiB") +
                size(variants_tsvs, "GiB") +
                size(alleles_tsvs, "GiB") +
                size(realigned_bams, "GiB")))
    }
    RuntimeAttr runtime_attr = select_first([runtime_override, runtime_default])

    runtime {
        docker: expansion_hunter_docker
        cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, runtime_default.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    }
}

task RunReviewer {
    input {
        String sample_id
        File realigned_bam
        File realigned_bam_index
        File vcf_gz
        File reference_fasta
        File? reference_fasta_index
        File variant_catalog_json
        String expansion_hunter_docker
        RuntimeAttr? runtime_override
    }

    output {
        File? reviewer_images_gz = "reviewer_images.tar.gz"
        File? sample_metrics = "${sample_id}_metrics.csv"
    }

    command <<<
        set -euxo pipefail

        REF="$(basename "~{reference_fasta}")"
        ln -s ~{reference_fasta} $REF
        ln -s ~{reference_fasta_index} $REF.fai

        gunzip -c ~{vcf_gz} > genotypes.vcf

        python /opt/str/get_locus_metrics.py \
            --sample-id ~{sample_id} \
            --realigned-bam ~{realigned_bam} \
            --reference $REF \
            --genotypes genotypes.vcf \
            --variant-catalog ~{variant_catalog_json} \
            --output ~{sample_id}_metrics.csv

        tar -czf reviewer_images.tar.gz ~{sample_id}_*.svg
    >>>

    RuntimeAttr runtime_default = object {
        cpu_cores: 1,
        mem_gb: 4,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1,
        disk_gb: 20 +
            (2 * ceil(
                size(realigned_bam, "GiB") +
                (size(vcf_gz, "GiB") * 10) +
                size(reference_fasta, "GiB") +
                size(variant_catalog_json, "GiB")))
    }
    RuntimeAttr runtime_attr = select_first([runtime_override, runtime_default])

    runtime {
        docker: expansion_hunter_docker
        cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, runtime_default.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    }
}
