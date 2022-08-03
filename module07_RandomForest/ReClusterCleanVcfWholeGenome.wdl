version 1.0

# Author: Xuefang Zhao <XZHAO12@mgh.harvard.edu>

import "Structs.wdl"
import "AnnotateGenomicContext.wdl" as AnnotateGenomicContext
import "ShardSVsByGenomicContext.wdl" as ShardSVsByGenomicContext
import "ReClusterCleanVcfSVs.wdl" as ReClusterCleanVcfSVs

# Workflow to annotate vcf file with genomic context

workflow ReClusterCleanVcfWholeGenome {
    input {
        Array[File] vcfs
        Array[File] vcf_indexes
        Array[String] prefixes
        Array[String] vid_prefixes

        File Repeat_Masks
        File Simple_Repeats
        File Segmental_Duplicates

        String svtype
        String Genomic_Context

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
    }

    scatter(i in range(length(vcfs))){
        call ReClusterCleanVcfSVs.ReClusterCleanVcfSVs as ReClusterCleanVcfSVs{
            input:
                vcf = vcfs[i],
                vcf_index = vcf_indexes[i],
                prefix = prefixes[i],
                vid_prefix = vid_prefixes[i],

                Repeat_Masks  = Repeat_Masks,
                Simple_Repeats = Simple_Repeats,
                Segmental_Duplicates  = Segmental_Duplicates,
                svtype = svtype,
                Genomic_Context = Genomic_Context,
                num_vids  = num_vids,
                num_samples = num_samples,
                dist  = dist,
                svsize = svsize,
                frac  = frac,
                sample_overlap = sample_overlap,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_benchmark_docker = sv_benchmark_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_override_vcf_to_bed = runtime_override_vcf_to_bed,
                runtime_attr_override_anno_gc = runtime_attr_override_anno_gc,
                runtime_attr_override_inte_gc = runtime_attr_override_inte_gc,
                runtime_override_extract_SV_sites = runtime_override_extract_SV_sites,
                runtime_attr_override_svtk_vcfcluster = runtime_attr_override_svtk_vcfcluster,
                runtime_attr_override_extract_svid = runtime_attr_override_extract_svid,
                runtime_attr_override_integrate_vcfs  = runtime_attr_override_integrate_vcfs
        }
    }

    output{
        Array[File] resolve_vcfs = ReClusterCleanVcfSVs.resolve_vcf
        Array[File] svid_annotation = ReClusterCleanVcfSVs.svid_annotation
        Array[File] out_vcfs = ReClusterCleanVcfSVs.out_vcf
        Array[File] out_vcf_idxes = ReClusterCleanVcfSVs.out_vcf_idx
    }
}



