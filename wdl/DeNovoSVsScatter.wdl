version 1.0
    
import "Structs.wdl"
import "Utils.wdl" as Utils

workflow DeNovoSVsScatter {

    input {
        File ped_input
        Array[File] vcf_files
        File disorder_input
        String chromosome
        File raw_proband
        File raw_parents
        File raw_depth_proband
        File raw_depth_parents
        File exclude_regions
        File sample_batches
        File batch_bincov_index
        File python_config
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_denovo
        RuntimeAttr? runtime_attr_vcf_to_bed
        RuntimeAttr? runtime_attr_merge_bed
    }
    
    Array[String] coverage_index_files = transpose(read_tsv(batch_bincov_index))[2]

    # Scatter genotyping over shards
    scatter (shard in vcf_files) {
        call Utils.VcfToBed as VcfToBed {
            input:
                vcf_file=shard,
                args="--info ALL --include-filters",
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_vcf_to_bed
        }

        call RunDeNovo {
            input:
                bed_input=VcfToBed.bed_output,
                ped_input=ped_input,
                vcf_input=shard,
                disorder_input=disorder_input,
                chromosome=chromosome,
                raw_proband=raw_proband,
                raw_parents=raw_parents,
                raw_depth_proband=raw_depth_proband,
                raw_depth_parents=raw_depth_parents,
                exclude_regions = exclude_regions,
                coverage_indices = coverage_index_files,
                sample_batches = sample_batches,
                batch_bincov_index = batch_bincov_index,
                python_config=python_config,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_denovo
        }   
    }

    call MergeBedFiles as MergeBedFilesAnnotated {
        input:
            bed_files = RunDeNovo.annotation_output,
            chromosome = chromosome,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_override = runtime_attr_merge_bed
    }

    call MergeBedFiles as MergeBedFilesFinal {
        input:
            bed_files = RunDeNovo.denovo_output,
            chromosome = chromosome,
            variant_interpretation_docker = variant_interpretation_docker,
            runtime_attr_override = runtime_attr_merge_bed
    }

    output {
        File per_chromosome_annotation_output_file = MergeBedFilesAnnotated.per_chromosome_denovo_output
        File per_chromosome_final_output_file = MergeBedFilesFinal.per_chromosome_denovo_output
    }
}

task RunDeNovo {

    input {
        File bed_input
        File ped_input
        File vcf_input
        File disorder_input
        String chromosome
        File raw_proband
        File raw_parents
        File raw_depth_proband
        File raw_depth_parents
        File exclude_regions
        Array[File] coverage_indices
        File batch_bincov_index
        File sample_batches
        File python_config
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float vcf_size = size(vcf_input, "GB")
    Float bed_size = size(bed_input, "GB")

    RuntimeAttr default_attr = object {
        mem_gb: 16,
        disk_gb: ceil(15 + vcf_size + bed_size * 1.5),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File denovo_output = "~{basename}.denovo.bed.gz"
        File annotation_output = "~{basename}.annotation.bed.gz"
    }

    String basename = basename(vcf_input, ".vcf.gz")
    command {
        set -exuo pipefail

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

        bcftools view ~{vcf_input} | grep -v ^## | bgzip -c > ~{basename}.noheader.vcf.gz
        python /src/denovo/denovo_svs.py \
            --bed ~{bed_input} \
            --ped ~{ped_input} \
            --vcf ~{basename}.noheader.vcf.gz \
            --disorder ~{disorder_input} \
            --out ~{basename}.annotation.bed \
            --out_de_novo ~{basename}.denovo.bed \
            --raw_proband ~{raw_proband} \
            --raw_parents ~{raw_parents} \
            --raw_depth_proband ~{raw_depth_proband} \
            --raw_depth_parents ~{raw_depth_parents} \
            --config ~{python_config} \
            --exclude_regions ~{exclude_regions} \
            --coverage ~{batch_bincov_index} \
            --sample_batches ~{sample_batches} \
            --verbose True

        bgzip ~{basename}.denovo.bed
        bgzip ~{basename}.annotation.bed
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task MergeBedFiles {

    input {
        Array[File] bed_files
        String chromosome
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float bed_files_size = size(bed_files, "GB")

    RuntimeAttr default_attr = object {
        mem_gb: 3.75,
        disk_gb: ceil(10 + (bed_files_size) * 2.0),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File per_chromosome_denovo_output = "~{chromosome}.denovo.merged.bed.gz"
    }

    command {
        set -exuo pipefail

        zcat ~{bed_files[0]} | awk 'NR==1' > ~{chromosome}.denovo.merged.bed
        zcat ~{sep=" " bed_files} | grep -v ^chrom >> ~{chromosome}.denovo.merged.bed
        bgzip ~{chromosome}.denovo.merged.bed
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}