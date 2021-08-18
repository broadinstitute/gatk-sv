version 1.0

import "CleanVcfChromosome.wdl" as CleanVcfChromosome
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow CleanVcf {
  input {
    String cohort_name

    Array[File] complex_genotype_vcfs
    Array[File] complex_resolve_bothside_pass_lists
    Array[File] complex_resolve_background_fail_lists
    File merged_ped_file

    File contig_list
    File allosome_fai
    Int max_shards_per_chrom_clean_vcf_step1
    Int min_records_per_shard_clean_vcf_step1
    Int samples_per_clean_vcf_step2_shard

    File? outlier_samples_list

    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    # overrides for mini tasks
    RuntimeAttr? runtime_override_concat_cleaned_vcfs

    # overrides for CleanVcfContig
    RuntimeAttr? runtime_override_clean_vcf_1a
    RuntimeAttr? runtime_override_clean_vcf_1b
    RuntimeAttr? runtime_override_clean_vcf_2
    RuntimeAttr? runtime_override_clean_vcf_3
    RuntimeAttr? runtime_override_clean_vcf_4
    RuntimeAttr? runtime_override_clean_vcf_5
    RuntimeAttr? runtime_override_stitch_fragmented_cnvs
    RuntimeAttr? runtime_override_final_cleanup
    RuntimeAttr? runtime_override_split_vcf_to_clean
    RuntimeAttr? runtime_override_combine_step_1_vcfs
    RuntimeAttr? runtime_override_combine_step_1_sex_chr_revisions
    RuntimeAttr? runtime_override_split_include_list
    RuntimeAttr? runtime_override_combine_clean_vcf_2
    RuntimeAttr? runtime_override_combine_revised_4
    RuntimeAttr? runtime_override_combine_multi_ids_4
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
        prefix=cohort_name,
        max_shards_per_chrom_step1=max_shards_per_chrom_clean_vcf_step1,
        min_records_per_shard_step1=min_records_per_shard_clean_vcf_step1,
        samples_per_step2_shard=samples_per_clean_vcf_step2_shard,
        outlier_samples_list=outlier_samples_list,
        linux_docker=linux_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_override_clean_vcf_1a=runtime_override_clean_vcf_1a,
        runtime_override_clean_vcf_1b=runtime_override_clean_vcf_1b,
        runtime_override_clean_vcf_2=runtime_override_clean_vcf_2,
        runtime_override_clean_vcf_3=runtime_override_clean_vcf_3,
        runtime_override_clean_vcf_4=runtime_override_clean_vcf_4,
        runtime_override_clean_vcf_5=runtime_override_clean_vcf_5,
        runtime_override_stitch_fragmented_cnvs=runtime_override_stitch_fragmented_cnvs,
        runtime_override_final_cleanup=runtime_override_final_cleanup,
        runtime_override_split_vcf_to_clean=runtime_override_split_vcf_to_clean,
        runtime_override_combine_step_1_vcfs=runtime_override_combine_step_1_vcfs,
        runtime_override_combine_step_1_sex_chr_revisions=runtime_override_combine_step_1_sex_chr_revisions,
        runtime_override_split_include_list=runtime_override_split_include_list,
        runtime_override_combine_clean_vcf_2=runtime_override_combine_clean_vcf_2,
        runtime_override_combine_revised_4=runtime_override_combine_revised_4,
        runtime_override_combine_multi_ids_4=runtime_override_combine_multi_ids_4
    }
  }

  call MiniTasks.ConcatVcfs as ConcatCleanedVcfs {
    input:
      vcfs=CleanVcfChromosome.out,
      vcfs_idx=CleanVcfChromosome.out_idx,
      naive=true,
      outfile_prefix="~{cohort_name}.cleaned",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_cleaned_vcfs
  }

  output {
    File cleaned_vcf = ConcatCleanedVcfs.concat_vcf
    File cleaned_vcf_index = ConcatCleanedVcfs.concat_vcf_idx
  }
}