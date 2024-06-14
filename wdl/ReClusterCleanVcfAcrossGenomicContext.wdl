version 1.0

# Author: Xuefang Zhao <XZHAO12@mgh.harvard.edu>

import "Structs.wdl"
import "AnnotateGenomicContext.wdl" as AnnotateGenomicContext
import "ReClusterCleanVcfSVsUnit.wdl" as ReClusterCleanVcfSVsUnit
import "TasksBenchmark.wdl" as tasks_benchmark
import "TasksMakeCohortVcf.wdl" as tasks_mcv

# Workflow to annotate vcf file with genomic context
workflow ReClusterCleanVcfAcrossGenomicContext {
    input {
        File vcf
        File vcf_index
        File Repeat_Masks
        File Simple_Repeats
        File Segmental_Duplicates
        File? annotated_SV_origin

        String prefix
        String vid_prefix
        String svtype
        Array[String] Genomic_Context_list
        Array[Array[Int]] size_list

        Int num_vids
        Int num_samples
        Int dist
        Int svsize
        Float frac
        Float sample_overlap

        String sv_base_mini_docker
        String sv_benchmark_docker
        String sv_pipeline_docker

        # overrides for MiniTasks
        RuntimeAttr? runtime_override_vcf_to_bed
        RuntimeAttr? runtime_attr_override_anno_gc
        RuntimeAttr? runtime_attr_override_inte_gc
        RuntimeAttr? runtime_override_extract_SV_sites
        RuntimeAttr? runtime_attr_override_svtk_vcfcluster
        RuntimeAttr? runtime_attr_override_extract_svid
        RuntimeAttr? runtime_attr_override_integrate_vcfs
        RuntimeAttr? runtime_attr_override_concat_vcfs
        RuntimeAttr? runtime_attr_override_sort_reclustered_vcf
        RuntimeAttr? runtime_attr_override_concat_clustered_SVID

    }

    if(!defined(annotated_SV_origin)){
        call AnnotateGenomicContext.AnnotateSVsWithGenomicContext as AnnotateSVsWithGenomicContext {
            input:
                vcf = vcf,
                vcf_index = vcf_index,
                Repeat_Masks = Repeat_Masks,
                Simple_Repeats = Simple_Repeats,
                Segmental_Duplicates = Segmental_Duplicates,

                sv_base_mini_docker  =sv_base_mini_docker,
                sv_benchmark_docker = sv_benchmark_docker,
                sv_pipeline_docker = sv_pipeline_docker,

                runtime_override_extract_SV_sites = runtime_override_extract_SV_sites,
                runtime_override_vcf_to_bed = runtime_override_vcf_to_bed,
                runtime_attr_override_anno_gc = runtime_attr_override_anno_gc,
                runtime_attr_override_inte_gc = runtime_attr_override_inte_gc
        }
    }

    File annotated_SVs = select_first([annotated_SV_origin, AnnotateSVsWithGenomicContext.annotated_SVs])

    scatter (i in range(length(Genomic_Context_list))) {
        call ReClusterCleanVcfSVsUnit.ReClusterCleanVcfSVs as ReClusterCleanVcfSVs{
            input:
                vcf = vcf,
                vcf_index = vcf_index,
                annotated_SV_origin = annotated_SVs,

                Repeat_Masks = Repeat_Masks,
                Simple_Repeats = Simple_Repeats,
                Segmental_Duplicates = Segmental_Duplicates,
                
                prefix = "~{prefix}_~{Genomic_Context_list[i]}",
                vid_prefix = "~{vid_prefix}_~{Genomic_Context_list[i]}",
                svtype = svtype,
                Genomic_Context = Genomic_Context_list[i],
                min_size = size_list[i][0],
                max_size = size_list[i][1],

                num_vids = num_vids,
                num_samples = num_samples,
                dist = dist,
                svsize = svsize,
                frac = frac,
                sample_overlap = sample_overlap,

                sv_base_mini_docker  =sv_base_mini_docker,
                sv_benchmark_docker = sv_benchmark_docker,
                sv_pipeline_docker = sv_pipeline_docker,

                runtime_attr_override_extract_svid = runtime_attr_override_extract_svid,
                runtime_override_vcf_to_bed = runtime_override_vcf_to_bed,

                runtime_override_vcf_to_bed = runtime_override_vcf_to_bed,
                runtime_attr_override_anno_gc = runtime_attr_override_anno_gc,
                runtime_attr_override_inte_gc = runtime_attr_override_inte_gc,
                runtime_override_extract_SV_sites = runtime_override_extract_SV_sites,
                runtime_attr_override_svtk_vcfcluster = runtime_attr_override_svtk_vcfcluster,
                runtime_attr_override_extract_svid = runtime_attr_override_extract_svid,
                runtime_attr_override_integrate_vcfs = runtime_attr_override_integrate_vcfs,
                runtime_attr_override_sort_reclustered_vcf= runtime_attr_override_sort_reclustered_vcf
        }
    }

    call tasks_mcv.ConcatVcfs as ConcatVcfs{
        input:
            vcfs = ReClusterCleanVcfSVs.reclustered_SV,
            vcfs_idx = ReClusterCleanVcfSVs.reclustered_SV_idx,
            outfile_prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_override_concat_vcfs
    }

    output{
        File reclustered_SV_genomic_context = ConcatVcfs.concat_vcf
        File reclustered_SV_genomic_context_idx = ConcatVcfs.concat_vcf_idx
    }
}


