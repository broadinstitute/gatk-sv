##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

## Copyright Broad Institute, 2020
## 
## This WDL pipeline implements Duphold 
##
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

version 1.0

import "Structs.wdl"
import "AnnotateILFeaturesPerSample.wdl" as anno_il
import "VaPoR.wdl" as vapor

workflow AnnotateILFeatures{
    input{
        File cleanVcf

        Array[String] prefix
        Array[String] samples
        #Array[String?] il_bams
        #Array[String?] il_bam_bais
        
        #Array[File?] pe_metrics
        #Array[File?] pe_indexes
        #Array[File?] sr_metrics
        #Array[File?] sr_indexes
        #Array[File?] rd_metrics
        #Array[File?] rd_indexes
        
        #Array[File?] ref_SegDups
        #Array[File?] ref_SimpReps
        #Array[File?] ref_RepMasks


        Array[File] raw_mantas
        Array[File] raw_whams
        Array[File] raw_melts

        Array[File] pacbio_seqs
        Array[File] pacbio_indexes

        File ref_fasta
        File ref_fai
        File ref_dict
        File contig_list

        Boolean requester_pays_crams = false
        Boolean run_genomic_context_anno = false
        Boolean run_extract_algo_evi = false
        Boolean run_duphold = false
        Boolean run_extract_gt_gq = true
        Boolean run_versus_raw_vcf = true
        Boolean run_rdpesr_anno = true
        Boolean run_vapor = true

        String rdpesr_benchmark_docker
        String vapor_docker
        String duphold_docker
        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_vapor 
        RuntimeAttr? runtime_attr_duphold
        RuntimeAttr? runtime_attr_rdpesr
        RuntimeAttr? runtime_attr_bcf2vcf
        RuntimeAttr? runtime_attr_LocalizeCram
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_SplitVcf
        RuntimeAttr? runtime_attr_ConcatBeds
        RuntimeAttr? runtime_attr_ConcatVcfs
        RuntimeAttr? runtime_inte_anno
        RuntimeAttr? runtime_attr_split_vcf
    }

    Array[String] contigs = transpose(read_tsv(contig_list))[0]

    scatter (i in range(length(prefix))){
        call split_per_sample_vcf{
            input:
                vcf = cleanVcf,
                sample = samples[i],
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_split_vcf
        }

        call anno_il.AnnoILFeaturesPerSample as anno_il_features{
            input:
                prefix = samples[i],
                vcf_file = split_per_sample_vcf.vcf_file,

                #il_bam = il_bams[i],
                #il_bam_bai = il_bam_bais[i],

                #pe_matrix = pe_metrics[i],
                #pe_index = pe_indexes[i],
                #sr_matrix = sr_metrics[i],
                #sr_index = sr_indexes[i],
                #rd_matrix = rd_metrics[i],
                #rd_index = rd_indexes[i],

                #ref_SegDup = ref_SegDups[i],
                #ref_SimpRep = ref_SimpReps[i],
                #ref_RepMask = ref_RepMasks[i],

                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                contig_list = contig_list,

                raw_vcfs = [raw_mantas[i],raw_whams[i],raw_melts[i]],
                raw_algorithms = ["manta","wham","melt"],

                rdpesr_benchmark_docker = rdpesr_benchmark_docker,
                duphold_docker = duphold_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                vapor_docker = vapor_docker,

                requester_pays_crams = requester_pays_crams,
                run_genomic_context_anno = run_genomic_context_anno,
                run_extract_algo_evi = run_extract_algo_evi,
                run_duphold = run_duphold,
                run_extract_gt_gq = run_extract_gt_gq,
                run_versus_raw_vcf = run_versus_raw_vcf,
                run_rdpesr_anno = run_rdpesr_anno,
                run_vapor = run_vapor,

                runtime_attr_duphold = runtime_attr_duphold,
                runtime_attr_rdpesr = runtime_attr_rdpesr,
                runtime_attr_bcf2vcf = runtime_attr_bcf2vcf,
                runtime_attr_vapor = runtime_attr_vapor,
                runtime_attr_LocalizeCram = runtime_attr_LocalizeCram,
                runtime_attr_vcf2bed = runtime_attr_vcf2bed,
                runtime_attr_SplitVcf = runtime_attr_SplitVcf,
                runtime_attr_ConcatBeds = runtime_attr_ConcatBeds,
                runtime_attr_ConcatVcfs = runtime_attr_ConcatVcfs
        }

        call vapor.VaPoR as vapor{
            input:
                prefix = samples[i],
                bam_or_cram_file = pacbio_seqs[i],
                bam_or_cram_index = pacbio_indexes[i],
                vcf_file = split_per_sample_vcf.vcf_file,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                contigs = contigs,
                vapor_docker = vapor_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_vapor = runtime_attr_vapor,
                runtime_attr_vcf2bed = runtime_attr_vcf2bed,
                runtime_attr_SplitVcf = runtime_attr_SplitVcf,
                runtime_attr_ConcatBeds = runtime_attr_ConcatBeds
        }

        call IntegrateAnno{
            input:
                prefix = prefix[i],
                sample = samples[i],
                bed           = anno_il_features.bed,
                gc_anno       = anno_il_features.GCAnno,
                duphold_il    = anno_il_features.duphold_vcf_il,
                duphold_il_le = anno_il_features.duphold_vcf_il_le,
                duphold_il_ri = anno_il_features.duphold_vcf_il_ri,
                pesr_anno     = anno_il_features.PesrAnno,
                rd_anno       = anno_il_features.RdAnno,
                rd_le_anno    = anno_il_features.RdAnno_le,
                rd_ri_anno    = anno_il_features.RdAnno_ri,
                gt_anno       = anno_il_features.GTGQ,
                info_anno     = anno_il_features.vcf_info,
                raw_manta     = anno_il_features.vs_raw[0],
                raw_wham      = anno_il_features.vs_raw[1],
                raw_melt      = anno_il_features.vs_raw[2],
                vapor_info    = vapor.vcf_out,
                rdpesr_benchmark_docker = rdpesr_benchmark_docker,
                runtime_attr_override = runtime_inte_anno
        }
    }
    output{
        Array[File] integrated_file = IntegrateAnno.anno_file
        Array[File] vcf_plots = vapor.vcf_plots
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
        File bed
        File? gc_anno
        File? duphold_il
        File? duphold_il_le
        File? duphold_il_ri
        File? rd_anno
        File? rd_le_anno
        File? rd_ri_anno
        File? pesr_anno
        File? info_anno
        File? gt_anno
        File raw_manta
        File raw_wham
        File raw_melt
        File? denovo
        File? vapor_info
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
        
        ~{if defined(rd_anno) then "zcat ~{rd_anno} | grep ~{sample}  > ~{sample}.rd_anno"  else ""}
        ~{if defined(pesr_anno) then "zcat ~{pesr_anno} | grep ~{sample} > ~{sample}.pesr_anno"  else ""}

        Rscript /src/integrate_annotations.R --bed ~{bed} --output ~{prefix}.anno.bed \
            --raw_manta ~{raw_manta} \
            --raw_wham ~{raw_wham} \
            --raw_melt ~{raw_melt} \
            ~{"--gc_anno " + gc_anno} \
            ~{"--duphold_il " + duphold_il} \
            ~{"--duphold_il_le " + duphold_il_le} \
            ~{"--duphold_il_ri " + duphold_il_ri} \
            ~{"--rd_le " + rd_le_anno} \
            ~{"--rd_ri " + rd_ri_anno} \
            #~{"--rd " + sample}.rd_anno \
            #~{"--pesr " + sample}.pesr_anno \
            ~{"--denovo " + denovo} \
            ~{"--gt " + gt_anno} 
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
