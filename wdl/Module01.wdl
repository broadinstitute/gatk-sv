##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/01_pesr_clustering_MMDLW/13/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

version 1.0

import "PESRClustering.wdl" as pesr
import "DepthClustering.wdl" as depth

workflow Module01 {
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
    File pesr_blacklist

    File? depth_blacklist
    Float? depth_blacklist_frac_max

    String depth_flags
    Float depth_frac
    File contigs

    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_pesr_cluster
    RuntimeAttr? runtime_attr_pesr_concat
    RuntimeAttr? runtime_attr_depth_cluster
    RuntimeAttr? runtime_attr_depth_concat
    RuntimeAttr? runtime_attr_depth_vcf
    RuntimeAttr? runtime_attr_rdtest_bed
  }

  if (defined(manta_vcfs)) {
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
        blacklist = pesr_blacklist,
        contigs = contigs,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_cluster = runtime_attr_pesr_cluster,
        runtime_attr_concat = runtime_attr_pesr_concat
    }
  }

  if (defined(delly_vcfs)) {
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
        blacklist = pesr_blacklist,
        contigs = contigs,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_cluster = runtime_attr_pesr_cluster,
        runtime_attr_concat = runtime_attr_pesr_concat
    }
  }

  if (defined(wham_vcfs)) {
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
        blacklist = pesr_blacklist,
        contigs = contigs,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_cluster = runtime_attr_pesr_cluster,
        runtime_attr_concat = runtime_attr_pesr_concat
    }
  }

  if (defined(melt_vcfs)) {
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
        blacklist = pesr_blacklist,
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
  	  blacklist = depth_blacklist,
      blacklist_frac_max = depth_blacklist_frac_max,
  	  flags = depth_flags,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_pipeline_docker = sv_pipeline_docker,
  	  runtime_attr_bed_cluster = runtime_attr_depth_cluster,
  	  runtime_attr_concat = runtime_attr_depth_concat,
  	  runtime_attr_depth_vcf = runtime_attr_depth_vcf,
  	  runtime_attr_rdtest_bed = runtime_attr_rdtest_bed
  }

  output {
    File depth_vcf = ClusterDepth.clustered_vcf
    File? manta_vcf = ClusterPESR_manta.clustered_vcf
    File? delly_vcf = ClusterPESR_delly.clustered_vcf
    File? wham_vcf = ClusterPESR_wham.clustered_vcf
    File? melt_vcf = ClusterPESR_melt.clustered_vcf
  }
}
