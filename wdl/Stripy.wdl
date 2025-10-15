version 1.0

workflow STRipyPipeline {
    input {
        File input_bam
        File input_bam_index
        String genome_build = "hg38"
        File reference_fasta
        String? locus
        String? sex
        File? custom_catalog
        String analysis = "standard"
        File? config
        Boolean output_vcf = true
        Boolean verbose = false
        String docker_image
        Int memory_gb = 8
        Int cpu = 2
    }

    call RunSTRipy {
        input:
            input_bam = input_bam,
            input_bam_index = input_bam_index,
            genome_build = genome_build,
            reference_fasta = reference_fasta,
            locus = locus,
            sex = sex,
            custom_catalog = custom_catalog,
            analysis = analysis,
            config = config,
            output_vcf = output_vcf,
            verbose = verbose,
            docker_image = docker_image,
            memory_gb = memory_gb,
            cpu = cpu
    }

    output {
        Array[File] output_files = RunSTRipy.output_files
    }
}

task RunSTRipy {
    input {
        File input_bam
        File input_bam_index
        String genome_build
        File reference_fasta
        String? locus
        String? sex
        File? custom_catalog
        String analysis
        File? config
        Boolean output_vcf
        Boolean verbose
        String docker_image
        Int cpu
        Int memory_gb
    }
    
    String output_dir = "STRipy_output"

    # Safety factor (you can adjust this as needed)
    Float safety_factor = 1.5
    Float input_size = size(input_bam, "GiB") + size(input_bam_index, "GiB") + size(reference_fasta, "GiB")
    Int provision_size_gb = 10 + ceil(input_size * safety_factor)

    command {
        # Run STRipy pipeline using our wrapper
        mkdir -p ${output_dir}

        stripy \
            --input ${input_bam} \
            --genome ${genome_build} \
            --reference ${reference_fasta} \
            --output ${output_dir} \
            --analysis ${analysis} \
            --output-json true \
            --output-tsv true \
            --output-html true \
            --output-vcf ${output_vcf} \
            --verbose ${verbose} \
            --num-threads ${cpu} \
            ${if defined(config) then "--base-config " + config else ""} \
            ${if defined(locus) then "--locus " + locus else ""} \
            ${if defined(sex) then "--sex " + sex else ""} \
            ${if defined(custom_catalog) then "--custom " + custom_catalog else ""}
    }

    runtime {
        disks: "local-disk ${provision_size_gb} GiB SSD"
        docker: docker_image
        memory: "${memory_gb}G"
        cpu: cpu
        maxRetries: 2
    }

    output {
        Array[File] output_files = glob("${output_dir}/*")
    }
}
