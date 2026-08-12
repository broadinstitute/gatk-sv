version 1.0

import "CleanVcfChromosome.wdl" as CleanVcfChromosome
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "MakeCohortVcfMetrics.wdl" as metrics

workflow CleanVcf {
  input {
    File complex_genotype_vcf
    File complex_resolve_bothside_pass_list
    File complex_resolve_background_fail_list
    File? outlier_samples_list
    File ped_file

    String chr_x
    String chr_y
    String cohort_name
    Int format_vcf_records_per_shard = 5000
    Int preprocess_records_per_shard = 5000
    Int postprocess_records_per_shard = 5000

    String contig
    File contig_list
    File allosome_fai

    File HERVK_reference
    File LINE1_reference
    File intron_reference

    # Module metrics parameters
    # Run module metrics workflow at the end - off by default to avoid resource errors
    Boolean? run_module_metrics
    File? primary_contigs_list  # required if run_module_metrics = true
    File? baseline_cluster_vcf  # baseline files are optional for metrics workflow
    File? baseline_complex_resolve_vcf
    File? baseline_complex_genotype_vcf
    File? baseline_cleaned_vcf
    File? combine_batches_merged_vcf  # intermediate merged VCFs are optional for metrics workflow
    File? resolve_complex_merged_vcf
    File? genotype_complex_merged_vcf

    String linux_docker
    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_create_ploidy
    RuntimeAttr? runtime_attr_format_to_clean_create_ploidy
    RuntimeAttr? runtime_attr_format_to_clean_scatter
    RuntimeAttr? runtime_attr_format_to_clean_format
    RuntimeAttr? runtime_attr_format_to_clean_concat
    RuntimeAttr? runtime_attr_scatter_preprocess
    RuntimeAttr? runtime_attr_preprocess
    RuntimeAttr? runtime_attr_concat_preprocess
    RuntimeAttr? runtime_attr_revise_overlapping_cnvs
    RuntimeAttr? runtime_attr_revise_large_cnvs
    RuntimeAttr? runtime_attr_revise_multiallelics
    RuntimeAttr? runtime_attr_scatter_postprocess
    RuntimeAttr? runtime_attr_postprocess
    RuntimeAttr? runtime_attr_concat_postprocess
    RuntimeAttr? runtime_override_drop_redundant_cnvs
    RuntimeAttr? runtime_override_sort_drop_redundant_cnvs
    RuntimeAttr? runtime_override_stitch_fragmented_cnvs
    RuntimeAttr? runtime_override_rescue_me_dels
    RuntimeAttr? runtime_attr_add_high_fp_rate_filters
    RuntimeAttr? runtime_attr_add_retro_del_filters
    RuntimeAttr? runtime_override_final_cleanup
    RuntimeAttr? runtime_attr_format_to_output_create_ploidy
    RuntimeAttr? runtime_attr_format_to_output_scatter
    RuntimeAttr? runtime_attr_format_to_output_format
    RuntimeAttr? runtime_attr_format_to_output_concat
    RuntimeAttr? runtime_override_concat_cleaned_vcfs
  }

  call CleanVcfChromosome.CleanVcfChromosome {
    input:
      vcf=complex_genotype_vcf,
      contig=contig,
      background_list=complex_resolve_background_fail_list,
      bothsides_pass_list=complex_resolve_bothside_pass_list,
      outlier_samples_list=outlier_samples_list,
      ped_file=ped_file,
      chr_x=chr_x,
      chr_y=chr_y,
      prefix="~{cohort_name}.~{contig}",
      format_vcf_records_per_shard=format_vcf_records_per_shard,
      preprocess_records_per_shard=preprocess_records_per_shard,
      postprocess_records_per_shard=postprocess_records_per_shard,
      contig_list=contig_list,
      allosome_fai=allosome_fai,
      HERVK_reference=HERVK_reference,
      LINE1_reference=LINE1_reference,
      intron_reference=intron_reference,
      gatk_docker=gatk_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_format_to_clean_create_ploidy=runtime_attr_format_to_clean_create_ploidy,
      runtime_attr_format_to_clean_scatter=runtime_attr_format_to_clean_scatter,
      runtime_attr_format_to_clean_format=runtime_attr_format_to_clean_format,
      runtime_attr_format_to_clean_concat=runtime_attr_format_to_clean_concat,
      runtime_attr_scatter_preprocess=runtime_attr_scatter_preprocess,
      runtime_attr_preprocess=runtime_attr_preprocess,
      runtime_attr_concat_preprocess=runtime_attr_concat_preprocess,
      runtime_attr_revise_overlapping_cnvs=runtime_attr_revise_overlapping_cnvs,
      runtime_attr_revise_large_cnvs=runtime_attr_revise_large_cnvs,
      runtime_attr_revise_multiallelics=runtime_attr_revise_multiallelics,
      runtime_attr_scatter_postprocess=runtime_attr_scatter_postprocess,
      runtime_attr_postprocess=runtime_attr_postprocess,
      runtime_attr_concat_postprocess=runtime_attr_concat_postprocess,
      runtime_override_drop_redundant_cnvs=runtime_override_drop_redundant_cnvs,
      runtime_override_sort_drop_redundant_cnvs=runtime_override_sort_drop_redundant_cnvs,
      runtime_override_stitch_fragmented_cnvs=runtime_override_stitch_fragmented_cnvs,
      runtime_override_rescue_me_dels=runtime_override_rescue_me_dels,
      runtime_attr_add_high_fp_rate_filters=runtime_attr_add_high_fp_rate_filters,
      runtime_attr_add_retro_del_filters=runtime_attr_add_retro_del_filters,
      runtime_override_final_cleanup=runtime_override_final_cleanup,
      runtime_attr_format_to_output_create_ploidy=runtime_attr_format_to_output_create_ploidy,
      runtime_attr_format_to_output_scatter=runtime_attr_format_to_output_scatter,
      runtime_attr_format_to_output_format=runtime_attr_format_to_output_format,
      runtime_attr_format_to_output_concat=runtime_attr_format_to_output_concat
  }

  output {
    File cleaned_vcf = CleanVcfChromosome.out
    File cleaned_vcf_index = CleanVcfChromosome.out_idx
  }
}
