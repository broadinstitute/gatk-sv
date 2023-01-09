#### Copyright (C) 2022 Ryan Collins & the Talkowski Laboratory
####
#### Workflow to compute and filter no-call GT rates from a GATK-SV VCF


version 1.0

import "Structs.wdl"
import "Utils.wdl" as Utils
import "EnforceMinNoCallRate.wdl" as enforce_min_no_call_rate

workflow EnforceMinNoCallRateVcfList {
  input {
    Array[File] vcf_list
    Array[File] vcf_idx_list
    Int records_per_shard

    Boolean always_shard_vcf = false
    Boolean? reannotate_ncrs_in_vcf
    Boolean exclude_CTX = true

    Array[String]? sample_subset_prefixes
    Array[File]? sample_subset_lists

    Array[String]? info_col_to_remove

    Array[String]? svtype_list
    Array[String]? ncr_filter_field
    Array[Float]? NCR_cff_list
    Float? global_max_ncr
    String? global_ncr_filter_field = "NCR"

    String sv_base_mini_docker
    String sv_pipeline_base_docker

    RuntimeAttr? runtime_attr_index_vcf
    RuntimeAttr? runtime_attr_shard_vcf
    RuntimeAttr? runtime_attr_split_vcf_by_type
    RuntimeAttr? runtime_attr_apply_ncr_filter
    RuntimeAttr? runtime_attr_concat_filtered_vcfs
    RuntimeAttr? runtime_attr_exclude_type_from_vcf


  }
  
  scatter(i in range(length(vcf_list))){
    call enforce_min_no_call_rate.EnforceMinNoCallRate as EnforceMinNoCallRate{
      input:
        vcf = vcf_list[i],
        vcf_idx = vcf_idx_list[i],
        records_per_shard = records_per_shard,
        always_shard_vcf = always_shard_vcf,
        reannotate_ncrs_in_vcf = reannotate_ncrs_in_vcf,

        sample_subset_prefixes = sample_subset_prefixes,
        sample_subset_lists = sample_subset_lists,
 
        info_col_to_remove = info_col_to_remove,
        svtype_list  = svtype_list,
        ncr_filter_field = ncr_filter_field,
        NCR_cff_list = NCR_cff_list,
        global_max_ncr = global_max_ncr,
        global_ncr_filter_field = global_ncr_filter_field,
        exclude_CTX = exclude_CTX,

        sv_pipeline_base_docker = sv_pipeline_base_docker,
        sv_base_mini_docker = sv_base_mini_docker,

        runtime_attr_index_vcf = runtime_attr_index_vcf,
        runtime_attr_shard_vcf = runtime_attr_shard_vcf,
        runtime_attr_split_vcf_by_type = runtime_attr_split_vcf_by_type,
        runtime_attr_apply_ncr_filter = runtime_attr_apply_ncr_filter,
        runtime_attr_concat_filtered_vcfs = runtime_attr_concat_filtered_vcfs, 
        runtime_attr_exclude_type_from_vcf = runtime_attr_exclude_type_from_vcf


    }
  }

  output{
    Array[File] ncr_annotated_vcf_list = EnforceMinNoCallRate.ncr_annotated_vcf
    Array[File] ncr_annotated_vcf_idx_list = EnforceMinNoCallRate.ncr_annotated_vcf_idx
    Array[File] ncr_table_list = EnforceMinNoCallRate.ncr_table
    Array[File?] ncr_filtered_vcf_list = EnforceMinNoCallRate.ncr_filtered_vcf
    Array[File?] ncr_filtered_vcf_idx_list = EnforceMinNoCallRate.ncr_filtered_vcf_idx
  }
}

