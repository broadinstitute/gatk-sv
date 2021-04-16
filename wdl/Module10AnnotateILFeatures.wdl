version 1.0

import "Structs.wdl"
import "AnnotateILFeatures.wdl" as anno_il

workflow Module10AnnoILFeatures {
    input {
        File cleanVcf

        Array[String] prefix
        Array[String] samples
        Array[String] il_bams
        Array[String] il_bam_bais
        
        Array[File] pe_metrics
        Array[File] pe_indexes
        Array[File] sr_metrics
        Array[File] sr_indexes
        Array[File] rd_metrics
        Array[File] rd_indexes
        Array[File] raw_mantas
        Array[File] raw_whams
        Array[File] raw_melts

        File denovo_list
        File ref_SegDup
        File ref_SimpRep
        File ref_RepMask

        File ref_fasta
        File ref_fai
        File ref_dict
        File contig_list

        String rdpesr_benchmark_docker
        String vapor_docker
        String duphold_docker
        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_Vapor 
        RuntimeAttr? runtime_attr_duphold
        RuntimeAttr? runtime_attr_rdpesr
        RuntimeAttr? runtime_attr_bcf2vcf
        RuntimeAttr? runtime_attr_LocalizeCram
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_SplitVcf
        RuntimeAttr? runtime_attr_ConcatBeds
        RuntimeAttr? runtime_attr_ConcatVcfs
        RuntimeAttr? runtime_inte_anno
    }

    scatter(i in range(length(prefix))){
        call split_per_sample_vcf{
            input:
                vcf = cleanVcf,
                sample = samples[i],
                sv_pipeline_docker = sv_pipeline_docker
        }
        call anno_il.AnnoILFeatures as anno_il_features{
            input:
                prefix = samples[i],
                il_bam = il_bams[i],
                il_bam_bai = il_bam_bais[i],
                vcf_file = split_per_sample_vcf.vcf_file,
                pe_matrix = pe_metrics[i],
                pe_index = pe_indexes[i],
                sr_matrix = sr_metrics[i],
                sr_index = sr_indexes[i],
                rd_matrix = rd_metrics[i],
                rd_index = rd_indexes[i],
                ref_SegDup = ref_SegDup,
                ref_SimpRep = ref_SimpRep,
                ref_RepMask = ref_RepMask,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                contig_list = contig_list,
                raw_vcfs = [raw_mantas[i],raw_whams[i],raw_melts[i]],
                raw_algorithms = ["manta","wham","melt"],

                rdpesr_benchmark_docker = rdpesr_benchmark_docker,
                vapor_docker = vapor_docker,
                duphold_docker = duphold_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,

                runtime_attr_Vapor = runtime_attr_Vapor,
                runtime_attr_duphold = runtime_attr_duphold,
                runtime_attr_rdpesr = runtime_attr_rdpesr,
                runtime_attr_bcf2vcf = runtime_attr_bcf2vcf,
                runtime_attr_LocalizeCram = runtime_attr_LocalizeCram,
                runtime_attr_vcf2bed = runtime_attr_vcf2bed,
                runtime_attr_SplitVcf = runtime_attr_SplitVcf,
                runtime_attr_ConcatBeds = runtime_attr_ConcatBeds,
                runtime_attr_ConcatVcfs = runtime_attr_ConcatVcfs
        }

        call IntegrateAnno{
            input:
                prefix = prefix[i],
                sample = samples[i],
                gc_anno = anno_il_features.GCAnno,
                duphold_il = anno_il_features.duphold_vcf_il,
                duphold_il_le = anno_il_features.duphold_vcf_il_le,
                duphold_il_ri = anno_il_features.duphold_vcf_il_ri,
                pesr_anno = anno_il_features.PesrAnno,
                rd_anno = anno_il_features.RdAnno,
                rd_le_anno = anno_il_features.RdAnno_le,
                rd_ri_anno = anno_il_features.RdAnno_ri,
                gt_anno = anno_il_features.GTGQ,
                info_anno = anno_il_features.vcf_info,
                raw_manta = anno_il_features.vs_raw[0],
                raw_wham = anno_il_features.vs_raw[1],
                raw_melt = anno_il_features.vs_raw[2],
                de_novo = denovo_list,
                rdpesr_benchmark_docker = rdpesr_benchmark_docker,
                runtime_attr_override = runtime_inte_anno
        }
    }

    output{
        Array[File] integrated_file = IntegrateAnno.anno_file
    }
}


task split_per_sample_vcf{
    input{
        File vcf
        String sample
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
        }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 3.75, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    output {
        File vcf_file = "~{sample}.vcf.gz"
        File vcf_idx ="~{sample}.vcf.gz.tbi"
    }
    command <<<

        set -Eeuo pipefail
        
        bcftools view -s ~{sample} ~{vcf} | grep -v "0/0" | bgzip > ~{sample}.vcf.gz
        tabix -p vcf ~{sample}.vcf.gz

    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task IntegrateAnno{
    input{
        File gc_anno
        File duphold_il
        File duphold_il_le
        File duphold_il_ri
        File rd_anno
        File rd_le_anno
        File rd_ri_anno
        File pesr_anno
        File info_anno
        File gt_anno
        File raw_manta
        File raw_wham
        File raw_melt
        File de_novo
        String prefix
        String sample
        String rdpesr_benchmark_docker
        RuntimeAttr? runtime_attr_override
        }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 3.75, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File anno_file = "~{prefix}.anno.bed.gz"
    }
    
    command <<<
        zcat ~{rd_anno} | grep ~{sample}  > tmp.rd_anno
        zcat ~{pesr_anno} | grep ~{sample} > tmp.pesr_anno
        Rscript /src/integrate_annotations.R \
            --gc_anno ~{gc_anno} \
            --duphold_il ~{duphold_il} \
            --duphold_il_le ~{duphold_il_le} \
            --duphold_il_ri ~{duphold_il_ri} \
            --rd_le ~{rd_le_anno} \
            --rd_ri ~{rd_ri_anno} \
            --rd tmp.rd_anno \
            --pesr tmp.pesr_anno \
            --info ~{info_anno} \
            --gt ~{gt_anno} \
            --raw_manta ~{raw_manta} \
            --raw_wham ~{raw_wham} \
            --raw_melt ~{raw_melt} \
            --denovo ~{de_novo} \
            --output ~{prefix}.anno.bed
        bgzip ~{prefix}.anno.bed
    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: rdpesr_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
