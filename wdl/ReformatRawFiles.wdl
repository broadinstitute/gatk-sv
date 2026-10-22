version 1.0
    
import "Structs.wdl"
import "Utils.wdl" as Utils
import "TasksMakeCohortVcf.wdl" as TasksMakeCohortVcf

workflow ReformatRawFiles {

    input {
        Array[String] contigs
        File raw_files_list
        File ped_input
        String variant_interpretation_docker
        String sv_base_mini_docker
        Boolean depth
        RuntimeAttr? runtime_attr_vcf_to_bed
        RuntimeAttr? runtime_attr_merge_bed
        RuntimeAttr? runtime_attr_divide_by_chrom
        RuntimeAttr? runtime_attr_reformat_bed
    }

    Array[String] raw_files = transpose(read_tsv(raw_files_list))[1]

    scatter(raw_file in raw_files) {
        call Utils.VcfToBed as VcfToBed {
            input:
                vcf_file=raw_file,
                args="--info SVTYPE",
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_vcf_to_bed
        }
    }

    call TasksMakeCohortVcf.CatUncompressedFiles as MergeBedFiles {
        input:
            shards=VcfToBed.bed_output,
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_merge_bed
    }

    scatter (contig in contigs) {
        call RawDivideByChrom {
            input:
                bed_file=MergeBedFiles.outfile,
                chromosome=contig,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override=runtime_attr_divide_by_chrom
        }

        if (depth) {
            call RawReformatBedDepth {
                input:
                    per_chromosome_bed_file=RawDivideByChrom.per_chromosome_bed_output,
                    ped_input=ped_input,
                    chromosome=contig,
                    variant_interpretation_docker=variant_interpretation_docker,
                    runtime_attr_override=runtime_attr_reformat_bed
            }
        }

        if (!(depth)) {
            call RawReformatBed {
                input:
                    per_chromosome_bed_file=RawDivideByChrom.per_chromosome_bed_output,
                    ped_input=ped_input,
                    chromosome=contig,
                    variant_interpretation_docker=variant_interpretation_docker,
                    runtime_attr_override=runtime_attr_reformat_bed
            }
        }
        File reformatted_parents_output_ = select_first([RawReformatBed.reformatted_parents_output, RawReformatBedDepth.reformatted_parents_depth_output])
        File reformatted_proband_output_ = select_first([RawReformatBed.reformatted_proband_output, RawReformatBedDepth.reformatted_proband_depth_output])
    }


    output {
        Array[File] reformatted_parents_raw_files = reformatted_parents_output_
        Array[File] reformatted_proband_raw_files = reformatted_proband_output_
    }
}

task RawDivideByChrom {

    input {
        File bed_file
        String chromosome
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(bed_file, "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
        mem_gb: 3.75,
        disk_gb: ceil(10 + input_size * 1.3),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File per_chromosome_bed_output = "${chromosome}.bed.gz"
    }

    command {
        set -exuo pipefail
        
        zcat ~{bed_file} | \
        awk '$1 == "~{chromosome}"' | \
        bgzip -c > ~{chromosome}.bed.gz
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

task RawReformatBed {

    input {
        File per_chromosome_bed_file
        File ped_input
        String chromosome
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float bed_file_size = size(per_chromosome_bed_file, "GB")
    Float ped_size = size(ped_input, "GB")

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(10 + ped_size + bed_file_size * 2.0),
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
       File reformatted_proband_output = "${chromosome}.proband.reformatted.sorted.bed.gz"
       File reformatted_parents_output = "${chromosome}.parents.reformatted.sorted.bed.gz"
    }

    command {
        set -exuo pipefail

        # Reformat bed file
        Rscript /src/denovo/reformat_raw_bed.R ~{per_chromosome_bed_file} ~{ped_input} ~{chromosome}.proband.reformatted.bed ~{chromosome}.parents.reformatted.bed
        bedtools sort -i ~{chromosome}.proband.reformatted.bed | bgzip -c > ~{chromosome}.proband.reformatted.sorted.bed.gz
        bedtools sort -i ~{chromosome}.parents.reformatted.bed | bgzip -c > ~{chromosome}.parents.reformatted.sorted.bed.gz
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

task RawReformatBedDepth {
    input {
        File per_chromosome_bed_file
        File ped_input
        String chromosome
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float bed_file_size = size(per_chromosome_bed_file, "GB")
    Float ped_size = size(ped_input, "GB")

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(10 + ped_size + bed_file_size * 2.0),
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File reformatted_proband_depth_output = "${chromosome}.proband.depth.reformatted.sorted.bed.gz"
        File reformatted_parents_depth_output = "${chromosome}.parents.depth.reformatted.sorted.bed.gz"
    }

    command {
        set -exuo pipefail

        # Reformat bed file
        Rscript /src/denovo/reformat_raw_bed.R ~{per_chromosome_bed_file} ~{ped_input} ~{chromosome}.proband.reformatted.bed ~{chromosome}.parents.reformatted.bed
        bedtools sort -i ~{chromosome}.proband.reformatted.bed | bgzip > ~{chromosome}.proband.depth.reformatted.sorted.bed.gz
        bedtools sort -i ~{chromosome}.parents.reformatted.bed | bgzip > ~{chromosome}.parents.depth.reformatted.sorted.bed.gz
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
