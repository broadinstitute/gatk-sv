version 1.0

import "CleanVcfChromosome.wdl" as CleanVcfChromosome
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "HailMerge.wdl" as HailMerge

workflow CleanVcf {
  input {
    String cohort_name

    Array[File] complex_genotype_vcfs
    Array[File] complex_resolve_bothside_pass_lists
    Array[File] complex_resolve_background_fail_lists
    File merged_ped_file

    File contig_list
    File allosome_fai
    Int max_shards_per_chrom_step1
    Int min_records_per_shard_step1
    Int samples_per_step2_shard
    Int? max_samples_per_shard_step3
    Int clean_vcf1b_records_per_shard
    Int clean_vcf5_records_per_shard

    String chr_x
    String chr_y

    File? outlier_samples_list

    Boolean use_hail = false
    String? gcs_project

    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_hail_docker
    String sv_pipeline_updates_docker

    # overrides for mini tasks
    RuntimeAttr? runtime_override_preconcat_clean_final
    RuntimeAttr? runtime_override_hail_merge_clean_final
    RuntimeAttr? runtime_override_fix_header_clean_final
    RuntimeAttr? runtime_override_concat_cleaned_vcfs
    RuntimeAttr? runtime_override_fix_bad_ends

    # overrides for CleanVcfContig
    RuntimeAttr? runtime_override_clean_vcf_1a
    RuntimeAttr? runtime_override_clean_vcf_2
    RuntimeAttr? runtime_override_clean_vcf_3
    RuntimeAttr? runtime_override_clean_vcf_4
    RuntimeAttr? runtime_override_clean_vcf_5_scatter
    RuntimeAttr? runtime_override_clean_vcf_5_make_cleangq
    RuntimeAttr? runtime_override_clean_vcf_5_find_redundant_multiallelics
    RuntimeAttr? runtime_override_clean_vcf_5_polish
    RuntimeAttr? runtime_override_stitch_fragmented_cnvs
    RuntimeAttr? runtime_override_final_cleanup

    # Clean vcf 1b
    RuntimeAttr? runtime_attr_override_subset_large_cnvs_1b
    RuntimeAttr? runtime_attr_override_sort_bed_1b
    RuntimeAttr? runtime_attr_override_intersect_bed_1b
    RuntimeAttr? runtime_attr_override_build_dict_1b
    RuntimeAttr? runtime_attr_override_scatter_1b
    RuntimeAttr? runtime_attr_override_filter_vcf_1b
    RuntimeAttr? runtime_override_concat_vcfs_1b
    RuntimeAttr? runtime_override_cat_multi_cnvs_1b

    RuntimeAttr? runtime_override_preconcat_step1
    RuntimeAttr? runtime_override_hail_merge_step1
    RuntimeAttr? runtime_override_fix_header_step1

    RuntimeAttr? runtime_override_preconcat_drc
    RuntimeAttr? runtime_override_hail_merge_drc
    RuntimeAttr? runtime_override_fix_header_drc

    RuntimeAttr? runtime_override_split_vcf_to_clean
    RuntimeAttr? runtime_override_combine_step_1_sex_chr_revisions
    RuntimeAttr? runtime_override_split_include_list
    RuntimeAttr? runtime_override_combine_clean_vcf_2
    RuntimeAttr? runtime_override_combine_revised_4
    RuntimeAttr? runtime_override_combine_multi_ids_4
    RuntimeAttr? runtime_override_drop_redundant_cnvs
    RuntimeAttr? runtime_override_combine_step_1_vcfs
    RuntimeAttr? runtime_override_sort_drop_redundant_cnvs
  }

  #Scatter per chromosome
  Array[String] contigs = transpose(read_tsv(contig_list))[0]
  scatter ( i in range(length(contigs)) ) {
    String contig = contigs[i]

    call CleanVcfChromosome.CleanVcfChromosome {
      input:
        vcf=complex_genotype_vcfs[i],
        contig=contig,
        background_list=complex_resolve_background_fail_lists[i],
        ped_file=merged_ped_file,
        bothsides_pass_list=complex_resolve_bothside_pass_lists[i],
        allosome_fai=allosome_fai,
        prefix="~{cohort_name}.~{contig}",
        max_shards_per_chrom_step1=max_shards_per_chrom_step1,
        min_records_per_shard_step1=min_records_per_shard_step1,
        samples_per_step2_shard=samples_per_step2_shard,
        max_samples_per_shard_step3=max_samples_per_shard_step3,
        outlier_samples_list=outlier_samples_list,
        use_hail=use_hail,
        gcs_project=gcs_project,
        clean_vcf1b_records_per_shard=clean_vcf1b_records_per_shard,
        clean_vcf5_records_per_shard=clean_vcf5_records_per_shard,
        chr_x=chr_x,
        chr_y=chr_y,
        linux_docker=linux_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_updates_docker=sv_pipeline_updates_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_hail_docker=sv_pipeline_hail_docker,
        runtime_override_clean_vcf_1a=runtime_override_clean_vcf_1a,
        runtime_override_clean_vcf_2=runtime_override_clean_vcf_2,
        runtime_override_clean_vcf_3=runtime_override_clean_vcf_3,
        runtime_override_clean_vcf_4=runtime_override_clean_vcf_4,
        runtime_override_clean_vcf_5_scatter=runtime_override_clean_vcf_5_scatter,
        runtime_override_clean_vcf_5_make_cleangq=runtime_override_clean_vcf_5_make_cleangq,
        runtime_override_clean_vcf_5_find_redundant_multiallelics=runtime_override_clean_vcf_5_find_redundant_multiallelics,
        runtime_override_clean_vcf_5_polish=runtime_override_clean_vcf_5_polish,
        runtime_override_stitch_fragmented_cnvs=runtime_override_stitch_fragmented_cnvs,
        runtime_override_final_cleanup=runtime_override_final_cleanup,
        runtime_override_split_vcf_to_clean=runtime_override_split_vcf_to_clean,
        runtime_override_combine_step_1_sex_chr_revisions=runtime_override_combine_step_1_sex_chr_revisions,
        runtime_override_split_include_list=runtime_override_split_include_list,
        runtime_override_combine_clean_vcf_2=runtime_override_combine_clean_vcf_2,
        runtime_override_combine_revised_4=runtime_override_combine_revised_4,
        runtime_override_combine_multi_ids_4=runtime_override_combine_multi_ids_4,
        runtime_override_preconcat_step1=runtime_override_preconcat_step1,
        runtime_override_hail_merge_step1=runtime_override_hail_merge_step1,
        runtime_override_fix_header_step1=runtime_override_fix_header_step1,
        runtime_override_preconcat_drc=runtime_override_preconcat_drc,
        runtime_override_hail_merge_drc=runtime_override_hail_merge_drc,
        runtime_override_fix_header_drc=runtime_override_fix_header_drc,
        runtime_override_drop_redundant_cnvs=runtime_override_drop_redundant_cnvs,
        runtime_attr_override_subset_large_cnvs_1b=runtime_attr_override_subset_large_cnvs_1b,
        runtime_attr_override_sort_bed_1b=runtime_attr_override_sort_bed_1b,
        runtime_attr_override_intersect_bed_1b=runtime_attr_override_intersect_bed_1b,
        runtime_attr_override_build_dict_1b=runtime_attr_override_build_dict_1b,
        runtime_attr_override_scatter_1b=runtime_attr_override_scatter_1b,
        runtime_attr_override_filter_vcf_1b=runtime_attr_override_filter_vcf_1b,
        runtime_override_concat_vcfs_1b=runtime_override_concat_vcfs_1b,
        runtime_override_cat_multi_cnvs_1b=runtime_override_cat_multi_cnvs_1b,
        runtime_override_fix_bad_ends=runtime_override_fix_bad_ends
    }
  }

  if (use_hail) {
    call HailMerge.HailMerge as ConcatVcfsHail {
      input:
        vcfs=CleanVcfChromosome.out,
        prefix="~{cohort_name}.cleaned",
        gcs_project=gcs_project,
        reset_cnv_gts=true,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_hail_docker=sv_pipeline_hail_docker,
        runtime_override_preconcat=runtime_override_preconcat_clean_final,
        runtime_override_hail_merge=runtime_override_hail_merge_clean_final,
        runtime_override_fix_header=runtime_override_fix_header_clean_final
    }
  }
  if (!use_hail) {
    call MiniTasks.ConcatVcfs as ConcatCleanedVcfs {
      input:
        vcfs=CleanVcfChromosome.out,
        vcfs_idx=CleanVcfChromosome.out_idx,
        allow_overlaps=true,
        outfile_prefix="~{cohort_name}.cleaned",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_concat_cleaned_vcfs
    }
  }

  output {
    File cleaned_vcf = select_first([ConcatCleanedVcfs.concat_vcf, ConcatVcfsHail.merged_vcf])
    File cleaned_vcf_index = select_first([ConcatCleanedVcfs.concat_vcf_idx, ConcatVcfsHail.merged_vcf_index])
  }
}