version 1.0

import "CleanVcfChromosome.wdl" as CleanVcfChromosome
import "TasksClusterBatch.wdl" as TasksCluster
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "MakeCohortVcfMetrics.wdl" as metrics

workflow CleanVcf {
  input {
    String cohort_name

    Array[File] complex_genotype_vcfs
    File complex_resolve_bothside_pass_list
    File complex_resolve_background_fail_list
    File ped_file

    File contig_list
    File allosome_fai

    File HERVK_reference
    File LINE1_reference
    File intron_reference

    String chr_x
    String chr_y

    File? outlier_samples_list

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

    String gatk_docker
    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_create_ploidy
    RuntimeAttr? runtime_attr_format_to_clean
		RuntimeAttr? runtime_attr_preprocess
		RuntimeAttr? runtime_attr_revise_overlapping_cnvs
		RuntimeAttr? runtime_attr_revise_large_cnvs
		RuntimeAttr? runtime_attr_revise_abnormal_allosomes
		RuntimeAttr? runtime_attr_revise_multiallelics
		RuntimeAttr? runtime_attr_postprocess
		RuntimeAttr? runtime_override_drop_redundant_cnvs
		RuntimeAttr? runtime_override_sort_drop_redundant_cnvs
		RuntimeAttr? runtime_override_stitch_fragmented_cnvs
		RuntimeAttr? runtime_override_rescue_me_dels
		RuntimeAttr? runtime_attr_add_high_fp_rate_filters
		RuntimeAttr? runtime_attr_add_retro_del_filters
		RuntimeAttr? runtime_override_final_cleanup
		RuntimeAttr? runtime_attr_format_to_output
    RuntimeAttr? runtime_override_concat_cleaned_vcfs
  }

  call TasksCluster.CreatePloidyTableFromPed {
    input:
      ped_file=ped_file,
      contig_list=contig_list,
      retain_female_chr_y=false,
      chr_x=chr_x,
      chr_y=chr_y,
      output_prefix="~{cohort_name}.ploidy",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_create_ploidy
  }

  #Scatter per chromosome
  Array[String] contigs = transpose(read_tsv(contig_list))[0]
  scatter ( i in range(length(contigs)) ) {
    String contig = contigs[i]
    call CleanVcfChromosome.CleanVcfChromosome {
      input:
        vcf=complex_genotype_vcfs[i],
        contig=contig,
        background_list=complex_resolve_background_fail_list,
        bothsides_pass_list=complex_resolve_bothside_pass_list,
        ped_file=ped_file,
        allosome_fai=allosome_fai,
        prefix="~{cohort_name}.~{contig}",
        outlier_samples_list=outlier_samples_list,
        ploidy_table=CreatePloidyTableFromPed.out,
        HERVK_reference=HERVK_reference,
        LINE1_reference=LINE1_reference,
        intron_reference=intron_reference,
        chr_x=chr_x,
        chr_y=chr_y,
        gatk_docker=gatk_docker,
        linux_docker=linux_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_format_to_clean=runtime_attr_format_to_clean,
        runtime_attr_preprocess=runtime_attr_preprocess,
        runtime_attr_revise_overlapping_cnvs=runtime_attr_revise_overlapping_cnvs,
        runtime_attr_revise_large_cnvs=runtime_attr_revise_large_cnvs,
        runtime_attr_revise_abnormal_allosomes=runtime_attr_revise_abnormal_allosomes,
        runtime_attr_revise_multiallelics=runtime_attr_revise_multiallelics,
        runtime_attr_postprocess=runtime_attr_postprocess,
        runtime_override_drop_redundant_cnvs=runtime_override_drop_redundant_cnvs,
        runtime_override_sort_drop_redundant_cnvs=runtime_override_sort_drop_redundant_cnvs,
        runtime_override_stitch_fragmented_cnvs=runtime_override_stitch_fragmented_cnvs,
        runtime_override_rescue_me_dels=runtime_override_rescue_me_dels,
        runtime_attr_add_high_fp_rate_filters=runtime_attr_add_high_fp_rate_filters,
        runtime_attr_add_retro_del_filters=runtime_attr_add_retro_del_filters,
        runtime_override_final_cleanup=runtime_override_final_cleanup,
        runtime_attr_format_to_output=runtime_attr_format_to_output
    }
  }

  call MiniTasks.ConcatVcfs as ConcatCleanedVcfs {
    input:
      vcfs=CleanVcfChromosome.out,
      vcfs_idx=CleanVcfChromosome.out_idx,
      allow_overlaps=true,
      outfile_prefix="~{cohort_name}.cleaned",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_cleaned_vcfs
  }

  Boolean run_module_metrics_ = if defined(run_module_metrics) then select_first([run_module_metrics]) else false
  if (run_module_metrics_) {
    call metrics.MakeCohortVcfMetrics {
      input:
        name = cohort_name,
        cluster_vcf = combine_batches_merged_vcf,
        complex_resolve_vcf = resolve_complex_merged_vcf,
        complex_genotype_vcf = genotype_complex_merged_vcf,
        cleaned_vcf = ConcatCleanedVcfs.concat_vcf,
        baseline_cluster_vcf = baseline_cluster_vcf,
        baseline_complex_resolve_vcf = baseline_complex_resolve_vcf,
        baseline_complex_genotype_vcf = baseline_complex_genotype_vcf,
        baseline_cleaned_vcf = baseline_cleaned_vcf,
        contig_list = select_first([primary_contigs_list]),
        linux_docker = linux_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_base_mini_docker = sv_base_mini_docker
    }
  }

  output {
    File cleaned_vcf = ConcatCleanedVcfs.concat_vcf
    File cleaned_vcf_index = ConcatCleanedVcfs.concat_vcf_idx
    File? metrics_file_makecohortvcf = MakeCohortVcfMetrics.metrics_file
  }
}