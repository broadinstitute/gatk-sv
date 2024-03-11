#this script is developed to fix false postiive inversions in gnomad v3 vcf that are >1Mb in size but have PE_GT=0
#INVs with no passing samples have filter column revised from "PASS" to "UNRESOLVED"

#developing workdir on erisone: /data/talkowski/xuefang/data/gnomad_V3/module08/step9_sm_depyh_only_dup_fix

version 1.0

import "Structs.wdl"
import "RerunAnnoPerContig.wdl" as rerun_anno_per_contig

workflow RerunAnno {

    input {
        File CNV_to_Fix #1 column file including the two CNVs that were revised from DEL / DUP and have GT format other than ./.
        File anno_id
        File ped_file
        File par_bed
        File? allosomes_list
        File sample_pop_assignments
        File two_subsets_af_header

        Array[File] vcf_list
        Array[File] vcf_idx_list
        Array[File] contig_list
        Array[File] prefix_list
        Array[Boolean] fix_mCNV_formats_list
        Array[File] two_subsets_af_info_list

        Int  sv_per_shard

        File protein_coding_gtf
        File? noncoding_bed
        Int? promoter_window
        Int? max_breakend_as_cnv_length
        String? svannotate_additional_args

        Boolean use_hail
        String? gcs_project

        String gatk_docker
        String sv_base_mini_docker
        String sv_pipeline_docker
        String sv_pipeline_base_docker
        String sv_pipeline_hail_docker

        RuntimeAttr? runtime_attr_shard_vcf
        RuntimeAttr? runtime_attr_compute_AFs
        RuntimeAttr? runtime_attr_combine_vcfs
        RuntimeAttr? runtime_attr_concat_vcfs
        RuntimeAttr? runtime_attr_svannotate
        RuntimeAttr? runtime_attr_scatter_vcf
        RuntimeAttr? runtime_attr_revise_cpx_func_anno
        RuntimeAttr? runtime_attr_concat_sharded_cluster
    }

    scatter (i in range(length(vcf_list))) {
        call rerun_anno_per_contig.RerunAnnoPerContig{
            input:
                CNV_to_Fix = CNV_to_Fix,
                anno_id =  anno_id,
                ped_file = ped_file,
                par_bed =  par_bed,
                allosomes_list = allosomes_list,
                sample_pop_assignments = sample_pop_assignments,
                two_subsets_af_header = two_subsets_af_header,

                vcf =  vcf_list[i],
                vcf_idx = vcf_idx_list[i],
                contig = contig_list[i],
                prefix = prefix_list[i],
                fix_mCNV_formats = fix_mCNV_formats_list[i],
                two_subsets_af_info = two_subsets_af_info_list[i],

                sv_per_shard = sv_per_shard,

                protein_coding_gtf = protein_coding_gtf,
                noncoding_bed = noncoding_bed,
                promoter_window = promoter_window,
                max_breakend_as_cnv_length = max_breakend_as_cnv_length,
                svannotate_additional_args = svannotate_additional_args,

                use_hail  = use_hail,
                gcs_project =  gcs_project,

                gatk_docker =  gatk_docker,
                sv_base_mini_docker =  sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                sv_pipeline_base_docker =  sv_pipeline_base_docker,
                sv_pipeline_hail_docker = sv_pipeline_hail_docker,

                runtime_attr_shard_vcf = runtime_attr_shard_vcf,
                runtime_attr_compute_AFs  = runtime_attr_compute_AFs,
                runtime_attr_combine_vcfs = runtime_attr_combine_vcfs,
                runtime_attr_concat_vcfs  = runtime_attr_concat_vcfs,
                runtime_attr_svannotate =  runtime_attr_svannotate,
                runtime_attr_scatter_vcf  = runtime_attr_scatter_vcf,
                runtime_attr_revise_cpx_func_anno = runtime_attr_revise_cpx_func_anno,
                runtime_attr_concat_sharded_cluster = runtime_attr_concat_sharded_cluster
        }
    }

    output {
        Array[File] re_annotated_vcf = RerunAnnoPerContig.output_vcf
        Array[File] re_annotated_vcf_idx= RerunAnnoPerContig.output_vcf_idx

        Array[File] re_annotated_sites = RerunAnnoPerContig.output_site
        Array[File] re_annotated_sites_idx = RerunAnnoPerContig.output_site_idx

        Array[File] re_annotated_bed = RerunAnnoPerContig.output_bed

        Array[File] re_annotated_with_2_subset_af_vcf = RerunAnnoPerContig.with_two_subset_af_vcf
        Array[File] re_annotated_with_2_subset_af_vcf_idx= RerunAnnoPerContig.with_two_subset_af_vcf_idx

        Array[File] re_annotated_with_2_subset_af_sites = RerunAnnoPerContig.with_two_subset_af_site
        Array[File] re_annotated_with_2_subset_af_sites_idx = RerunAnnoPerContig.with_two_subset_af_site_idx

    }
  }



