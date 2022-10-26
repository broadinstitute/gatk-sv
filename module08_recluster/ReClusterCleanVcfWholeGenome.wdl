version 1.0

# Author: Xuefang Zhao <XZHAO12@mgh.harvard.edu>

import "Structs.wdl"
import "ReClusterCleanVcfAcrossSVTYPE.wdl" as ReClusterCleanVcfAcrossSVTYPE
# Workflow to annotate vcf file with genomic context

workflow ReClusterCleanVcfWholeGenome {
    input {
        Array[File] vcf_list
        Array[File] vcf_index_list
        Array[String] chromosome_name_list
        String prefix

        File Repeat_Masks
        File Simple_Repeats
        File Segmental_Duplicates
        File empty_file

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
        RuntimeAttr? runtime_attr_override_concat_clustered_SVID
        RuntimeAttr? runtime_attr_override_sort_reclustered_vcf
        RuntimeAttr? runtime_attr_override_concat_vcfs
        RuntimeAttr? runtime_override_concat_sharded_cluster

    }

    scatter(i in range(length(vcf_list))){
        call ReClusterCleanVcfAcrossSVTYPE.ReClusterCleanVcfAcrossSVTYPE as ReClusterCleanVcfAcrossSVTYPE{
            input:
                vcf = vcf_list[i],
                vcf_index = vcf_index_list[i],
                prefix = "~{prefix}.~{chromosome_name_list[i]}",
                cohort_name = "Redun_Removed",
                contig = chromosome_name_list[i],

                empty_file = empty_file,
                Repeat_Masks  = Repeat_Masks,
                Simple_Repeats = Simple_Repeats,
                Segmental_Duplicates  = Segmental_Duplicates,

                svtype_list = svtype_list,
                Genomic_Context_list_list = Genomic_Context_list_list,
                size_list_list = size_list_list,

                num_samples = num_samples,
                dist  = dist,
                svsize = svsize,
                frac  = frac,
                sample_overlap = sample_overlap,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_benchmark_docker = sv_benchmark_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                sv_pipeline_hail_docker = sv_pipeline_hail_docker,

                runtime_override_vcf_to_bed = runtime_override_vcf_to_bed,
                runtime_attr_override_anno_gc = runtime_attr_override_anno_gc,
                runtime_attr_override_inte_gc = runtime_attr_override_inte_gc,
                runtime_override_extract_SV_sites = runtime_override_extract_SV_sites,
                runtime_attr_override_svtk_vcfcluster = runtime_attr_override_svtk_vcfcluster,
                runtime_attr_override_extract_svid = runtime_attr_override_extract_svid,
                runtime_attr_override_integrate_vcfs  = runtime_attr_override_integrate_vcfs,
                runtime_attr_override_concat_vcfs = runtime_attr_override_concat_vcfs,
                runtime_attr_override_sort_reclustered_vcf = runtime_attr_override_sort_reclustered_vcf,
                runtime_attr_override_concat_clustered_SVID = runtime_attr_override_concat_clustered_SVID,
                runtime_override_concat_sharded_cluster = runtime_override_concat_sharded_cluster
       }
    }

    output{
        Array[File] svid_annotation = ReClusterCleanVcfAcrossSVTYPE.svid_annotation
        Array[File] out_vcfs = ReClusterCleanVcfAcrossSVTYPE.reclustered_SV_svtype
        Array[File] out_vcf_idxes = ReClusterCleanVcfAcrossSVTYPE.reclustered_SV_svtype_idx
    }
}



