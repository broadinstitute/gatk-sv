version 1.0

import "Structs.wdl"
import "AnnoRdPeSr.wdl" as anno_pesrrd

workflow Module10AnnotateRdPeSr {
    input {
        Array[String] prefixes
        Array[String] samples
        
        Array[File] beds
        Array[File] bed_le_flanks
        Array[File] bed_ri_flanks

        Array[File] pe_metrics
        Array[File] pe_indexes
        Array[File] sr_metrics
        Array[File] sr_indexes
        Array[File] rd_metrics
        Array[File] rd_indexes

        File contig_list

        String rdpesr_benchmark_docker
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

    scatter(i in range(length(prefixes))){
        call anno_pesrrd.AnnoRdPeSr as anno_rd_pe_sr{
            input:
                prefix = prefixes[i],
                sample = samples[i],
                bed = beds[i],
                bed_le_flank = bed_le_flanks[i],
                bed_ri_flank = bed_ri_flanks[i],
                pe_matrix = pe_metrics[i],
                pe_index = pe_indexes[i],
                sr_matrix = sr_metrics[i],
                sr_index = sr_indexes[i],
                rd_matrix = rd_metrics[i],
                rd_index = rd_indexes[i],
                contig_list = contig_list,
                rdpesr_benchmark_docker=rdpesr_benchmark_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker

        }
     }

    output{
        Array[File] rd_anno = anno_rd_pe_sr.RdAnno
        Array[File] rd_anno_le = anno_rd_pe_sr.RdAnno_le
        Array[File] rd_anno_ri = anno_rd_pe_sr.RdAnno_ri
        Array[File] pesr_anno = anno_rd_pe_sr.PesrAnno
    }
}


