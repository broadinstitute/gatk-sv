version 1.0

# Author: Xuefang Zhao <XZHAO12@mgh.harvard.edu>

import "Structs.wdl"
import "AnnotateGenomicContext.wdl" as AnnotateGenomicContext
import "ReClusterCleanVcfAcrossGenomicContext.wdl" as ReClusterCleanVcfAcrossGenomicContext
import "TasksBenchmark.wdl" as mini_tasks
# Workflow to annotate vcf file with genomic context
workflow ReClusterCleanVcfAcrossSVTYPE {
    input {
        File vcf
        File vcf_index
        File Repeat_Masks
        File Simple_Repeats
        File Segmental_Duplicates
        File? annotated_SV_origin
        File empty_file


        String contig
        String cohort_name
        String prefix
        Array[String] svtype_list
        Array[Array[String]] Genomic_Context_list_list
        Array[Array[Array[Int]]] size_list_list

        Int num_samples
        Int dist
        Int svsize
        Float frac
        Float sample_overlap

        String sv_base_mini_docker
        String sv_benchmark_docker
        String sv_pipeline_docker
        String sv_pipeline_hail_docker

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
        RuntimeAttr? runtime_override_concat_sharded_cluster
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

    scatter (i in range(length(svtype_list))) {
        call ReClusterCleanVcfAcrossGenomicContext.ReClusterCleanVcfAcrossGenomicContext as ReClusterCleanVcfAcrossGenomicContext{
            input:
                vcf = vcf,
                vcf_index = vcf_index,
                annotated_SV_origin = annotated_SVs,

                Repeat_Masks = Repeat_Masks,
                Simple_Repeats = Simple_Repeats,
                Segmental_Duplicates = Segmental_Duplicates,

                empty_file = empty_file,
                cohort_name = "~{cohort_name}.~{svtype_list[i]}",
                contig = contig,
                
                prefix = "~{prefix}_~{svtype_list[i]}",
                svtype = svtype_list[i],
                Genomic_Context_list = Genomic_Context_list_list[i],
                size_list = size_list_list[i],

                num_samples = num_samples,
                dist = dist,
                svsize = svsize,
                frac = frac,
                sample_overlap = sample_overlap,

                sv_pipeline_docker = sv_pipeline_docker,
                sv_base_mini_docker  =sv_base_mini_docker,
                sv_benchmark_docker = sv_benchmark_docker,
                sv_pipeline_hail_docker = sv_pipeline_hail_docker,

                runtime_attr_override_extract_svid = runtime_attr_override_extract_svid,
                runtime_override_vcf_to_bed = runtime_override_vcf_to_bed,

                runtime_override_vcf_to_bed = runtime_override_vcf_to_bed,
                runtime_attr_override_anno_gc = runtime_attr_override_anno_gc,
                runtime_attr_override_inte_gc = runtime_attr_override_inte_gc,
                runtime_override_extract_SV_sites = runtime_override_extract_SV_sites,
                runtime_attr_override_svtk_vcfcluster = runtime_attr_override_svtk_vcfcluster,
                runtime_attr_override_extract_svid = runtime_attr_override_extract_svid,
                runtime_attr_override_integrate_vcfs = runtime_attr_override_integrate_vcfs,
                runtime_attr_override_sort_reclustered_vcf= runtime_attr_override_sort_reclustered_vcf,
                runtime_override_concat_sharded_cluster = runtime_override_concat_sharded_cluster

        }
    }

 
    call mini_tasks.ConcatVcfs as ConcatVcfs{
        input:
            vcfs = ReClusterCleanVcfAcrossGenomicContext.reclustered_SV_genomic_context,
            vcfs_idx = ReClusterCleanVcfAcrossGenomicContext.reclustered_SV_genomic_context_idx,
            outfile_prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_override_concat_vcfs
    }

    call mini_tasks.IntegrateReClusterdVcfs{
        input:
            vcf_all = vcf,
            vcf_all_idx = vcf_index,
            vcf_recluster = ConcatVcfs.concat_vcf,
            vcf_recluster_idx = ConcatVcfs.concat_vcf_idx,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_override_integrate_vcfs
    }

    call mini_tasks.SortReClusterdVcfs{
        input:
            vcf_1 = IntegrateReClusterdVcfs.reclustered_Part1,
            vcf_2 = IntegrateReClusterdVcfs.reclustered_Part2,
            vcf_1_idx = IntegrateReClusterdVcfs.reclustered_Part1_idx,
            vcf_2_idx = IntegrateReClusterdVcfs.reclustered_Part2_idx,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_override_sort_reclustered_vcf
    }



    output{
        File svid_annotation = IntegrateReClusterdVcfs.SVID_anno
        File reclustered_SV_svtype = SortReClusterdVcfs.sorted_vcf
        File reclustered_SV_svtype_idx = SortReClusterdVcfs.sorted_vcf_idx
    }
}


