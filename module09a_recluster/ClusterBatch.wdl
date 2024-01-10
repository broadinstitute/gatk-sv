version 1.0

import "PESRClustering.wdl" as pesr
import "DepthClustering.wdl" as depth
import "ClusterBatchMetrics.wdl" as metrics
import "TasksClusterBatch.wdl" as tasks
import "Utils.wdl" as util
import "PlotSVCountsPerSample.wdl" as sv_counts

workflow ClusterBatch {
  input {
    String batch

    File? manta_vcf_tar
    File? wham_vcf_tar
    File? melt_vcf_tar
    File? scramble_vcf_tar

    File del_bed
    File dup_bed
    File ped_file

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

    # Depth-based variant clustering
    Int? depth_records_per_bed_shard
    File depth_exclude_intervals
    Float depth_exclude_overlap_fraction
    Float depth_interval_overlap
    String? depth_clustering_algorithm

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
    String? sv_pipeline_base_docker  # required if run_module_metrics = true
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
      vcf_tar=select_first([manta_vcf_tar, melt_vcf_tar, wham_vcf_tar]),
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

  if (defined(manta_vcf_tar)) {
    call pesr.ClusterPESR as ClusterPESR_manta {
      input:
        vcf_tar=select_first([manta_vcf_tar]),
        ploidy_table=CreatePloidyTableFromPed.out,
        batch=batch,
        caller="manta",
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

  if (defined(wham_vcf_tar)) {
    call pesr.ClusterPESR as ClusterPESR_wham {
      input:
        vcf_tar=select_first([wham_vcf_tar]),
        ploidy_table=CreatePloidyTableFromPed.out,
        batch=batch,
        caller="wham",
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
        runtime_attr_svcluster=runtime_attr_svcluster_wham,
        runtime_override_concat_vcfs_pesr=runtime_override_concat_vcfs_pesr,
        runtime_attr_gatk_to_svtk_vcf=runtime_attr_gatk_to_svtk_vcf_pesr,
        runtime_attr_exclude_intervals_pesr=runtime_attr_exclude_intervals_pesr
    }
  }

  if (defined(melt_vcf_tar)) {
    call pesr.ClusterPESR as ClusterPESR_melt {
      input:
        vcf_tar=select_first([melt_vcf_tar]),
        ploidy_table=CreatePloidyTableFromPed.out,
        batch=batch,
        caller="melt",
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
        runtime_attr_svcluster=runtime_attr_svcluster_melt,
        runtime_override_concat_vcfs_pesr=runtime_override_concat_vcfs_pesr,
        runtime_attr_gatk_to_svtk_vcf=runtime_attr_gatk_to_svtk_vcf_pesr,
        runtime_attr_exclude_intervals_pesr=runtime_attr_exclude_intervals_pesr
    }
  }

  if (defined(scramble_vcf_tar)) {
    call pesr.ClusterPESR as ClusterPESR_scramble {
      input:
        vcf_tar=select_first([scramble_vcf_tar]),
        ploidy_table=CreatePloidyTableFromPed.out,
        batch=batch,
        caller="scramble",
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
        runtime_attr_svcluster=runtime_attr_svcluster_scramble,
        runtime_override_concat_vcfs_pesr=runtime_override_concat_vcfs_pesr,
        runtime_attr_gatk_to_svtk_vcf=runtime_attr_gatk_to_svtk_vcf_pesr,
        runtime_attr_exclude_intervals_pesr=runtime_attr_exclude_intervals_pesr
    }
  }

  call depth.ClusterDepth {
  	input:
      del_bed=del_bed,
      dup_bed=dup_bed,
      batch=batch,
      ploidy_table=CreatePloidyTableFromPed.out,
      contig_list=contig_list,
      contig_subset_list=contig_subset_list,
      sample_list=GetSampleIdsFromVcfTar.out_file,
      records_per_bed_shard=select_first([depth_records_per_bed_shard, 1000000]),
      exclude_intervals=depth_exclude_intervals,
      exclude_overlap_fraction=depth_exclude_overlap_fraction,
      clustering_algorithm=depth_clustering_algorithm,
      depth_interval_overlap=depth_interval_overlap,
      gatk_to_svtk_script=gatk_to_svtk_script,
      cnv_bed_to_gatk_vcf_script=cnv_bed_to_gatk_vcf_script,
      java_mem_fraction=java_mem_fraction,
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_fasta_fai,
      reference_dict=reference_dict,
      gatk_docker=gatk_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_cnv_bed_to_gatk_vcf=runtime_attr_cnv_bed_to_gatk_vcf,
      runtime_attr_scatter_bed=runtime_attr_scatter_bed,
      runtime_attr_exclude_intervals_depth=runtime_attr_exclude_intervals_depth,
      runtime_attr_svcluster=runtime_attr_svcluster_depth,
      runtime_attr_gatk_to_svtk_vcf=runtime_attr_gatk_to_svtk_vcf_depth,
      runtime_override_concat_vcfs=runtime_override_concat_vcfs_depth
  }

  Boolean run_module_metrics_ = if defined(run_module_metrics) then select_first([run_module_metrics]) else true
  if (run_module_metrics_) {
    call metrics.ClusterBatchMetrics {
      input:
        name = batch,
        depth_vcf = ClusterDepth.clustered_vcf,
        manta_vcf = ClusterPESR_manta.clustered_vcf,
        wham_vcf = ClusterPESR_wham.clustered_vcf,
        melt_vcf = ClusterPESR_melt.clustered_vcf,
        scramble_vcf = ClusterPESR_scramble.clustered_vcf,
        baseline_depth_vcf = baseline_depth_vcf,
        baseline_manta_vcf = baseline_manta_vcf,
        baseline_wham_vcf = baseline_wham_vcf,
        baseline_scramble_vcf = baseline_scramble_vcf,
        baseline_melt_vcf = baseline_melt_vcf,
        contig_list = select_first([contig_list]),
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_base_docker = select_first([sv_pipeline_base_docker]),
        linux_docker = select_first([linux_docker])
    }
  }

  if (defined(N_IQR_cutoff_plotting)){
    call sv_counts.PlotSVCountsPerSample {
    input:
      prefix = batch,
      vcfs = [ClusterDepth.clustered_vcf, ClusterPESR_manta.clustered_vcf, ClusterPESR_wham.clustered_vcf, ClusterPESR_melt.clustered_vcf, ClusterPESR_scramble.clustered_vcf],
      N_IQR_cutoff = select_first([N_IQR_cutoff_plotting]),
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_count_svs = runtime_attr_count_svs,
      runtime_attr_plot_svcounts = runtime_attr_plot_svcounts,
      runtime_attr_cat_outliers_preview = runtime_attr_cat_outliers_preview
    }
  }

  output {
    File clustered_depth_vcf = ClusterDepth.clustered_vcf
    File clustered_depth_vcf_index = ClusterDepth.clustered_vcf_index
    File? clustered_manta_vcf = ClusterPESR_manta.clustered_vcf
    File? clustered_manta_vcf_index = ClusterPESR_manta.clustered_vcf_index
    File? clustered_wham_vcf = ClusterPESR_wham.clustered_vcf
    File? clustered_wham_vcf_index = ClusterPESR_wham.clustered_vcf_index
    File? clustered_melt_vcf = ClusterPESR_melt.clustered_vcf
    File? clustered_melt_vcf_index = ClusterPESR_melt.clustered_vcf_index
    File? clustered_scramble_vcf = ClusterPESR_scramble.clustered_vcf
    File? clustered_scramble_vcf_index = ClusterPESR_scramble.clustered_vcf_index
    Array[File]? clustered_sv_counts = PlotSVCountsPerSample.sv_counts
    Array[File]? clustered_sv_count_plots = PlotSVCountsPerSample.sv_count_plots
    File? clustered_outlier_samples_preview = PlotSVCountsPerSample.outlier_samples_preview
    File? clustered_outlier_samples_with_reason = PlotSVCountsPerSample.outlier_samples_with_reason
    Int? clustered_num_outlier_samples = PlotSVCountsPerSample.num_outlier_samples
    File? metrics_file_clusterbatch = ClusterBatchMetrics.metrics_file
  }


}