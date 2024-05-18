version 1.0

import "PESRClustering.wdl" as pesr
import "TasksClusterBatch.wdl" as tasks
import "Utils.wdl" as util

workflow ClusterBatch {
  input {
    String batch

    File? vcf_tar
    String caller
    # Reference
    File contig_list
    File reference_fasta
    File reference_fasta_fai
    File reference_dict
    String? chr_x
    String? chr_y

    # For testing
    File? contig_subset_list  # subsets contig scatters
    File? ploidy_table_script
    File? cnv_bed_to_gatk_vcf_script
    File? svtk_to_gatk_script
    File? gatk_to_svtk_script


    # PESR-based variant clustering
    Int? pesr_min_size
    File pesr_exclude_intervals
    Float pesr_interval_overlap
    Int pesr_breakend_window
    String? pesr_clustering_algorithm

    # PlotSVCountsPerSample
    Int? N_IQR_cutoff_plotting

    # Module metrics parameters
    # Run module metrics workflow at the end - on by default
    Boolean? run_module_metrics
    String? linux_docker  # required if run_module_metrics = true
    File? baseline_depth_vcf  # baseline files are optional for metrics workflow
    File? baseline_manta_vcf
    File? baseline_wham_vcf
    File? baseline_melt_vcf
    File? baseline_scramble_vcf

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_ids_from_vcf_list
    RuntimeAttr? runtime_attr_create_ploidy
    RuntimeAttr? runtime_attr_prepare_pesr_vcfs
    RuntimeAttr? runtime_attr_svcluster_manta
    RuntimeAttr? runtime_attr_svcluster_melt
    RuntimeAttr? runtime_attr_svcluster_scramble
    RuntimeAttr? runtime_attr_svcluster_wham
    RuntimeAttr? runtime_override_concat_vcfs_pesr
    RuntimeAttr? runtime_attr_gatk_to_svtk_vcf_pesr
    RuntimeAttr? runtime_attr_scatter_bed
    RuntimeAttr? runtime_attr_cnv_bed_to_gatk_vcf
    RuntimeAttr? runtime_attr_exclude_intervals_depth
    RuntimeAttr? runtime_attr_svcluster_depth
    RuntimeAttr? runtime_attr_gatk_to_svtk_vcf_depth
    RuntimeAttr? runtime_override_concat_vcfs_depth
    RuntimeAttr? runtime_attr_exclude_intervals_pesr
    RuntimeAttr? runtime_attr_count_svs
    RuntimeAttr? runtime_attr_plot_svcounts
    RuntimeAttr? runtime_attr_cat_outliers_preview
  }

  call util.GetSampleIdsFromVcfTar {
    input:
      vcf_tar=vcf_tar,
      prefix="~{batch}.samples",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_ids_from_vcf_list
  }

  # TODO : properly set allosome ploidy, which creates problems in RDTest for allosomes at the moment
  call tasks.CreatePloidyTableFromPed {
    input:
      ped_file=ped_file,
      script=ploidy_table_script,
      contig_list=contig_list,
      retain_female_chr_y=true,
      chr_x=chr_x,
      chr_y=chr_y,
      output_prefix="~{batch}.ploidy",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_create_ploidy
  }

  if (defined(vcf_tar)) {
    call pesr.ClusterPESR as ClusterPESR_vcf {
      input:
        vcf_tar=select_first([vcf_tar]),
        ploidy_table=CreatePloidyTableFromPed.out,
        batch=batch,
        caller=caller,
        min_size=select_first([pesr_min_size, 50]),
        svtk_to_gatk_script=svtk_to_gatk_script,
        gatk_to_svtk_script=gatk_to_svtk_script,
        exclude_intervals=pesr_exclude_intervals,
        contig_list=contig_list,
        contig_subset_list=contig_subset_list,
        pesr_interval_overlap=pesr_interval_overlap,
        pesr_breakend_window=pesr_breakend_window,
        clustering_algorithm=pesr_clustering_algorithm,
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        gatk_docker=gatk_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_prepare_pesr_vcfs=runtime_attr_prepare_pesr_vcfs,
        runtime_attr_svcluster=runtime_attr_svcluster_manta,
        runtime_override_concat_vcfs_pesr=runtime_override_concat_vcfs_pesr,
        runtime_attr_gatk_to_svtk_vcf=runtime_attr_gatk_to_svtk_vcf_pesr,
        runtime_attr_exclude_intervals_pesr=runtime_attr_exclude_intervals_pesr
    }
  }



  output {
    File? clustered_vcf = ClusterPESR_vcf.clustered_vcf
    File? clustered_vcf_index = ClusterPESR_vcf.clustered_vcf_index
  }


}
