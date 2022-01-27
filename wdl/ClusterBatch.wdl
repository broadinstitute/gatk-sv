version 1.0

import "PESRClustering.wdl" as pesr
import "DepthClustering.wdl" as depth
import "ClusterBatchMetrics.wdl" as metrics

workflow ClusterBatch {
  input {
    Array[File]? manta_vcfs
    Array[File]? delly_vcfs
    Array[File]? wham_vcfs
    Array[File]? melt_vcfs
    File del_bed
    File dup_bed
    String batch
    Int pesr_svsize
    Float pesr_frac
    String pesr_flags
    Int pesr_distance
    File pesr_exclude_list

    File? depth_exclude_list
    Float? depth_exclude_list_frac_max

    String depth_flags
    Float depth_frac
    File contigs

    # Module metrics parameters
    # Run module metrics workflow at the end - on by default
    Boolean? run_module_metrics
    String? linux_docker  # required if run_module_metrics = true
    String? sv_pipeline_base_docker  # required if run_module_metrics = true
    File? primary_contigs_list  # required if run_module_metrics = true
    File? baseline_depth_vcf  # baseline files are optional for metrics workflow
    File? baseline_delly_vcf
    File? baseline_manta_vcf
    File? baseline_wham_vcf
    File? baseline_melt_vcf

    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_pesr_cluster
    RuntimeAttr? runtime_attr_pesr_concat
    RuntimeAttr? runtime_attr_depth_cluster
    RuntimeAttr? runtime_attr_depth_concat
    RuntimeAttr? runtime_attr_depth_vcf
    RuntimeAttr? runtime_attr_rdtest_bed
  }

  if (defined(manta_vcfs) && (length(select_first([manta_vcfs])) > 0)) {
    call pesr.ClusterPESR as ClusterPESR_manta {
      input:
        algorithm = "manta",
        vcfs = select_first([manta_vcfs]),
        svsize = pesr_svsize,
        frac = pesr_frac,
        svtypes = "DEL,DUP,INV,BND,INS",
        flags = pesr_flags,
        batch = batch,
        dist = pesr_distance,
        exclude_list = pesr_exclude_list,
        contigs = contigs,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_cluster = runtime_attr_pesr_cluster,
        runtime_attr_concat = runtime_attr_pesr_concat
    }
  }

  if (defined(delly_vcfs) && (length(select_first([delly_vcfs])) > 0)) {
    call pesr.ClusterPESR as ClusterPESR_delly {
      input:
        algorithm = "delly",
        vcfs = select_first([delly_vcfs]),
        svsize = pesr_svsize,
        frac = pesr_frac,
        svtypes = "DEL,DUP,INV,BND,INS",
        flags = pesr_flags,
        batch = batch,
        dist = pesr_distance,
        exclude_list = pesr_exclude_list,
        contigs = contigs,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_cluster = runtime_attr_pesr_cluster,
        runtime_attr_concat = runtime_attr_pesr_concat
    }
  }

  if (defined(wham_vcfs) && (length(select_first([wham_vcfs])) > 0)) {
    call pesr.ClusterPESR as ClusterPESR_wham {
      input:
        algorithm = "wham",
        vcfs = select_first([wham_vcfs]),
        svsize = pesr_svsize,
        frac = pesr_frac,
        svtypes = "DEL,DUP,INV,BND,INS",
        flags = pesr_flags,
        batch = batch,
        dist = pesr_distance,
        exclude_list = pesr_exclude_list,
        contigs = contigs,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_cluster = runtime_attr_pesr_cluster,
        runtime_attr_concat = runtime_attr_pesr_concat
    }
  }

  if (defined(melt_vcfs) && (length(select_first([melt_vcfs])) > 0)) {
    call pesr.ClusterPESR as ClusterPESR_melt {
      input:
        algorithm = "melt",
        vcfs = select_first([melt_vcfs]),
        svsize = pesr_svsize,
        frac = pesr_frac,
        svtypes = "DEL,DUP,INV,BND,INS",
        flags = pesr_flags,
        batch = batch,
        dist = pesr_distance,
        exclude_list = pesr_exclude_list,
        contigs = contigs,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_cluster = runtime_attr_pesr_cluster,
        runtime_attr_concat = runtime_attr_pesr_concat
    }
  }

  call depth.ClusterDepth as ClusterDepth {
  	input: 
  	  del_bed = del_bed,
  	  dup_bed = dup_bed,
  	  batch = batch,
  	  contigs = contigs,
  	  frac = depth_frac,
  	  exclude_list = depth_exclude_list,
      exclude_list_frac_max = depth_exclude_list_frac_max,
  	  flags = depth_flags,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_pipeline_docker = sv_pipeline_docker,
  	  runtime_attr_bed_cluster = runtime_attr_depth_cluster,
  	  runtime_attr_concat = runtime_attr_depth_concat,
  	  runtime_attr_depth_vcf = runtime_attr_depth_vcf,
  	  runtime_attr_rdtest_bed = runtime_attr_rdtest_bed
  }

  Boolean run_module_metrics_ = if defined(run_module_metrics) then select_first([run_module_metrics]) else true
  if (run_module_metrics_) {
    call metrics.ClusterBatchMetrics {
      input:
        name = batch,
        depth_vcf = ClusterDepth.clustered_vcf,
        delly_vcf = ClusterPESR_delly.clustered_vcf,
        manta_vcf = ClusterPESR_manta.clustered_vcf,
        wham_vcf = ClusterPESR_wham.clustered_vcf,
        melt_vcf = ClusterPESR_melt.clustered_vcf,
        baseline_depth_vcf = baseline_depth_vcf,
        baseline_delly_vcf = baseline_delly_vcf,
        baseline_manta_vcf = baseline_manta_vcf,
        baseline_wham_vcf = baseline_wham_vcf,
        baseline_melt_vcf = baseline_melt_vcf,
        contig_list = select_first([primary_contigs_list]),
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_base_docker = select_first([sv_pipeline_base_docker]),
        linux_docker = select_first([linux_docker])
    }
  }

  output {
    File depth_vcf = ClusterDepth.clustered_vcf
    File depth_vcf_index = ClusterDepth.clustered_vcf_index
    File? manta_vcf = ClusterPESR_manta.clustered_vcf
    File? manta_vcf_index = ClusterPESR_manta.clustered_vcf_index
    File? delly_vcf = ClusterPESR_delly.clustered_vcf
    File? delly_vcf_index = ClusterPESR_delly.clustered_vcf_index
    File? wham_vcf = ClusterPESR_wham.clustered_vcf
    File? wham_vcf_index = ClusterPESR_wham.clustered_vcf_index
    File? melt_vcf = ClusterPESR_melt.clustered_vcf
    File? melt_vcf_index = ClusterPESR_melt.clustered_vcf_index

    File? metrics_file_clusterbatch = ClusterBatchMetrics.metrics_file
  }
}
