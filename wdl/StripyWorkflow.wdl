version 1.0

import "Utils.wdl" as utils
import "Structs.wdl"

workflow StripyWorkflow {
    input {
        File bam_or_cram_file
        File? bam_or_cram_index
        File ped_file
        String? genome_build
        File reference_fasta
        String sample_name
        String? locus
        File? custom_catalog
        String? analysis
        File? config
        Boolean? verbose
        String stripy_docker
        String linux_docker
        RuntimeAttr? runtime_attr_override
    }

    Boolean is_bam_ = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)
    String index_ext_ = if is_bam_ then ".bai" else ".crai"
    File bam_or_cram_index_ = select_first([bam_or_cram_index, bam_or_cram_file + index_ext_])

    call utils.GetSampleSex {
        input:
            ped_file = ped_file,
            sample_id = sample_name,
            unknown_sex = "female",
            linux_docker = linux_docker
    }

    call RunStripy {
        input:
            bam_or_cram_file = bam_or_cram_file,
            bam_or_cram_index = bam_or_cram_index_,
            genome_build = genome_build,
            reference_fasta = reference_fasta,
            sample_name = sample_name,
            locus = locus,
            sex = GetSampleSex.out_string,
            custom_catalog = custom_catalog,
            analysis = analysis,
            config = config,
            verbose = verbose,
            stripy_docker = stripy_docker,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File json_output = RunStripy.json_output
        File tsv_output = RunStripy.tsv_output
        File html_output = RunStripy.html_output
        File? vcf_output = RunStripy.vcf_output
    }
}

task RunStripy {
    input {
        File bam_or_cram_file
        File bam_or_cram_index
        String genome_build = "hg38"
        File reference_fasta
        String sample_name
        String locus = "AFF2,AR,ARX_1,ARX_2,ATN1,ATXN1,ATXN10,ATXN2,ATXN3,ATXN7,ATXN8OS,BEAN1,C9ORF72,CACNA1A,CBL,CNBP,COMP,DAB1,DIP2B,DMD,DMPK,FGF14,FMR1,FOXL2,FXN,GIPC1,GLS,HOXA13_1,HOXA13_2,HOXA13_3,HOXD13,HTT,JPH3,LRP12,MARCHF6,NIPA1,NOP56,NOTCH2NLC,NUTM2B-AS1,PABPN1,PHOX2B,PPP2R2B,PRDM12,RAPGEF2,RFC1,RILPL1,RUNX2,SAMD12,SOX3,STARD7,TBP,TBX1,TCF4,TNRC6A,XYLT1,YEATS2,ZIC2,ZIC3"
        String? sex
        File? custom_catalog
        String analysis = "standard"
        File? config
        Boolean verbose = false
        String stripy_docker
        RuntimeAttr? runtime_attr_override
    }
    
    String output_dir = "STRipy_output"

    # Safety factor (you can adjust this as needed)
    Float safety_factor = 1.5
    Float input_size = size(bam_or_cram_file, "GiB") + size(bam_or_cram_index, "GiB") + size(reference_fasta, "GiB")
    Int provision_size_gb = 50 + ceil(input_size * safety_factor)
    RuntimeAttr default_attr = object {
                                   cpu_cores: 2,
                                   mem_gb: 8,
                                   disk_gb: provision_size_gb,
                                   boot_disk_gb: 10,
                                   preemptible_tries: 3,
                                   max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String bam_filename = basename(bam_or_cram_file)
    String bam_base_default = sub(bam_filename, "\\\.[^.]+$", "")

    String json_path = output_dir + "/" + sample_name + ".json"
    String tsv_path = output_dir + "/" + sample_name + ".tsv"
    String html_path = output_dir + "/" + sample_name + ".html"
    String vcf_path = output_dir + "/" + sample_name + ".vcf"

    command <<<
        # Run STRipy pipeline using our wrapper
        set -euxo pipefail
        mkdir -p ~{output_dir}
        stripy \
            --input ~{bam_or_cram_file} \
            --genome ~{genome_build} \
            --reference ~{reference_fasta} \
            --output ~{output_dir} \
            --analysis ~{analysis} \
            --output-json true \
            --output-tsv true \
            --output-html true \
            --output-vcf true \
            --verbose ~{verbose} \
            --num-threads $(nproc) \
            ~{if defined(config) then "--base-config " + config else ""} \
            ~{if defined(locus) then "--locus " + locus else ""} \
            ~{if defined(sex) then "--sex " + sex else ""} \
            ~{if defined(custom_catalog) then "--custom " + custom_catalog else ""}

        ACTUAL_FILENAME=$(basename "~{bam_or_cram_file}")
        ACTUAL_BASE=$(echo "${ACTUAL_FILENAME}" | sed 's/\.[^.]*$//')
        TARGET_BASE='~{sample_name}'
        echo "ACTUAL_FILENAME: ${ACTUAL_FILENAME}"
        echo "ACTUAL_BASE: ${ACTUAL_BASE}"
        echo "TARGET_BASE: ${TARGET_BASE}"
        ls ~{output_dir}
        if [ -f "~{output_dir}/${ACTUAL_FILENAME}.json" ]; then
            mv "~{output_dir}/${ACTUAL_FILENAME}.json" "~{output_dir}/${TARGET_BASE}.json"
        fi
        if [ -f "~{output_dir}/${ACTUAL_FILENAME}.tsv" ]; then
            mv "~{output_dir}/${ACTUAL_FILENAME}.tsv" "~{output_dir}/${TARGET_BASE}.tsv"
        fi
        if [ -f "~{output_dir}/${ACTUAL_FILENAME}.html" ]; then
            mv "~{output_dir}/${ACTUAL_FILENAME}.html" "~{output_dir}/${TARGET_BASE}.html"
        fi
        if [ -f "~{output_dir}/${ACTUAL_BASE}.vcf" ]; then
            mv "~{output_dir}/${ACTUAL_BASE}.vcf" "~{output_dir}/${TARGET_BASE}.vcf"
        fi
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: stripy_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        File json_output = json_path
        File tsv_output = tsv_path
        File html_output = html_path
        File vcf_output = vcf_path
    }
}
