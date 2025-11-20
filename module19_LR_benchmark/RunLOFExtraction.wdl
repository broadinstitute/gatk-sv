version 1.0

import "Structs.wdl"
import "RunLOFExtractionPerChrom.wdl" as RunLOFExtractionPerChrom

workflow RunLOFExtraction {
    input {
        Array[File] vcfs
        Array[File] vcf_indexes
        Array[String] prefix_list
        String output_prefix
        File extract_script    # the python script extract_SV_gene.lof.py
        Int chunk_size = 100000
        String sv_pipeline_base_docker
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_concat_bed
        RuntimeAttr? runtime_attr_concat_bed_2
        RuntimeAttr? runtime_attr_extract_lof
    }

    scatter (i in range(length(vcfs))) {
        call SplitVCF {
            input:
                vcf = vcfs[i],
                index = vcf_indexes[i],
                chunk_size = chunk_size,
                sv_pipeline_base_docker = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_split_vcf
        }

        call RunLOFExtractionPerChrom.RunLOFExtractionPerChrom as run_lof_extraction_per_chrom{
            input:
                vcfs = SplitVCF.split_vcfs,
                vcf_indexes = SplitVCF.split_indexes,
                extract_script = extract_script,
                prefix= prefix_list[i],
                sv_pipeline_base_docker = sv_pipeline_base_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_extract_lof = runtime_attr_extract_lof,
                runtime_attr_concat_bed = runtime_attr_concat_bed
        }
    }

    call RunLOFExtractionPerChrom.ConcatBeds as concat_beds{
        input:
            shard_bed_files = run_lof_extraction_per_chrom.bed_file,
            prefix = output_prefix,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat_bed_2
    }

    output {
        File output_bed = concat_beds.merged_bed_file
    }

}

################################################################################
# TASK 1 — Split VCF into ~100k-variant chunks + index each chunk
################################################################################

task SplitVCF {
    input {
        File vcf
        File index
        Int chunk_size
        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 10.0,
        disk_gb: ceil(10.0 + vcf * 3.0),
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10
    }   

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])


    command <<<
        set -euo pipefail

        mkdir splits
        cd splits

        # Count total variants
        total=$(bcftools index --nrecords ../${vcf})

        # Compute number of chunks
        chunks=$(( (total + ${chunk_size} - 1) / ${chunk_size} ))

        for i in $(seq 0 $((chunks - 1))); do
            start=$(( i * ${chunk_size} + 1 ))
            end=$(( (i + 1) * ${chunk_size} ))

            bcftools view -H ../${vcf} | \
                sed -n "${start},${end}p" | \
                bcftools view -O z -o chunk_${i}.vcf.gz

            bcftools index chunk_${i}.vcf.gz
        done
    >>>

    output {
        Array[File] split_vcfs = glob("splits/chunk_*.vcf.gz")
        Array[File] split_indexes = glob("splits/chunk_*.vcf.gz.tbi")
    }

    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_base_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

}

