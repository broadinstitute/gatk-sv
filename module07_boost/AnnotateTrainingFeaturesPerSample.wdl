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
import "TasksBenchmark.wdl" as mini_tasks
import "VaPoR.wdl" as vapor

workflow AnnotateTrainingFeaturesPerSample{
    input{
        String sample

        File IL_anno
        
        File ref_fasta
        File ref_fai
        File ref_dict
        Int min_shard_size

        File? vapor_result
        File? pacbio_seq
        File? pacbio_seq_index
        File? ont_query
        File? hgsv_query
        File? array_query
        File? pacbio_query
        File? bionano_query

        Boolean requester_pays_crams = false
        Boolean run_genomic_context_anno = false
        Boolean run_extract_algo_evi = false
        Boolean run_duphold = false
        Boolean run_extract_gt_gq = true
        Boolean run_versus_raw_vcf = true
        Boolean run_rdpesr_anno = true
        Boolean run_vapor = true

        String vapor_docker
        String sv_base_mini_docker
        String sv_benchmark_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_vapor 
        RuntimeAttr? runtime_attr_bcf2vcf
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_SplitVcf
        RuntimeAttr? runtime_attr_ConcatBeds
        RuntimeAttr? runtime_attr_bed_vs_ont
        RuntimeAttr? runtime_attr_bed_vs_hgsv
        RuntimeAttr? runtime_attr_bed_vs_array
        RuntimeAttr? runtime_attr_bed_vs_pacbio
        RuntimeAttr? runtime_attr_bed_vs_bionano
        RuntimeAttr? runtime_attr_override_inte_anno
        RuntimeAttr? runtime_attr_override_format_ref 
    }

    if (defined(pacbio_seq)){
        File bam_or_cram_file = select_first([pacbio_seq])
        File bam_or_cram_index = select_first([pacbio_seq_index])

        call vapor.VaPoR as VaPoR{
            input:
                prefix = sample,
                sample = sample,
                bam_or_cram_file = bam_or_cram_file,
                bam_or_cram_index = bam_or_cram_index,
                bed_file = IL_anno,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                min_shard_size = min_shard_size,
                vapor_docker = vapor_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                sv_base_mini_docker = sv_base_mini_docker,

                runtime_attr_vapor = runtime_attr_vapor,
                runtime_attr_bcf2vcf = runtime_attr_bcf2vcf,
                runtime_attr_vcf2bed = runtime_attr_vcf2bed,
                runtime_attr_SplitVcf = runtime_attr_SplitVcf,
                runtime_attr_ConcatBeds = runtime_attr_ConcatBeds
        }
    }

    if (defined(vapor_result) || defined(pacbio_seq)) {
       File? vapor_param = select_first([vapor_result, VaPoR.bed_out])
    }

    call FormatRef{
        input:
            bed_file = IL_anno,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_override_format_ref
    }

    if (defined(hgsv_query)){
        call mini_tasks.BedComparison as BedComparison_vs_hgsv{
            input:
                query = hgsv_query,
                ref = FormatRef.ref,
                prefix = "${sample}.vs.hgsv",
                sv_benchmark_docker=sv_benchmark_docker,
                runtime_attr_override = runtime_attr_bed_vs_hgsv
        }
    }

    if (defined(array_query)){
        call mini_tasks.BedComparison as BedComparison_vs_array{
            input:
                query = array_query,
                ref = FormatRef.ref,
                prefix = "${sample}.vs.array",
                sv_benchmark_docker=sv_benchmark_docker,
                runtime_attr_override = runtime_attr_bed_vs_array
        }
    }

    if (defined(pacbio_query)){
        call mini_tasks.BedComparison as BedComparison_vs_pacbio{
            input:
                query = pacbio_query,
                ref = FormatRef.ref,
                prefix = "${sample}.vs.pacbio",
                sv_benchmark_docker=sv_benchmark_docker,
                runtime_attr_override = runtime_attr_bed_vs_pacbio
        }
    }

    if (defined(ont_query)){
        call mini_tasks.BedComparison as BedComparison_vs_ont{
            input:
                query = ont_query,
                ref = FormatRef.ref,
                prefix = "${sample}.vs.ont",
                sv_benchmark_docker=sv_benchmark_docker,
                runtime_attr_override = runtime_attr_bed_vs_ont
        }
    }

    if (defined(bionano_query)){
        call mini_tasks.BedComparison as BedComparison_vs_bionano{
            input:
                query = bionano_query,
                ref = FormatRef.ref,
                prefix = "${sample}.vs.bionano",
                sv_benchmark_docker=sv_benchmark_docker,
                runtime_attr_override = runtime_attr_bed_vs_bionano
        }
    }

    call mini_tasks.IntegrateAnno{
        input:
            bed = IL_anno,
            prefix = sample,
            sample = sample,
            vs_ont =  BedComparison_vs_ont.comparison,
            vs_array = BedComparison_vs_array.comparison,
            vs_hgsv = BedComparison_vs_hgsv.comparison,
            vs_pacbio = BedComparison_vs_pacbio.comparison,
            vs_bionano = BedComparison_vs_bionano.comparison,
            vapor = vapor_param,
            sv_benchmark_docker = sv_benchmark_docker,
            runtime_attr_override = runtime_attr_override_inte_anno
    }


    output{
        File annotated_file = IntegrateAnno.anno_file
    }
}



task FormatRef{
    input{
        File bed_file
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 1, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File ref = "~{prefix}.ref.gz"
    }

    String prefix = basename(bed_file,'bed.gz')
    command <<<
        zcat ~{bed_file} | cut -f1-4,7-10 | bgzip  > ~{prefix}.ref.gz
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



