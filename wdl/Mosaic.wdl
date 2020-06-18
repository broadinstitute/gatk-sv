# To obtaina list of likely mosaic variants that failed RF due to separation only 

version 1.0

import "Structs.wdl"
import "DepthMosaic.wdl" as depth_mosaic
import "PreRFCohort.wdl" as preRF 
import "PesrMosaicPart1.wdl" as mosaic_pesr_part1
import "PesrMosaicPart2.wdl" as mosaic_pesr_part2

workflow MosaicManualCheck{
  input{
    File fam_file
    Int rare_cutoff
    File outlier
    Array[File] per_batch_clustered_pesr_vcf_list # preRF 
    Array[File] clustered_depth_vcfs
    Array[File] coverage_files
    Array[File] coverage_file_idxs
    Array[File] median_files
    
    Array[File] agg_metrics
    Array[File] RF_cutoffs
    
    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker
  }
  scatter (i in range(length(per_batch_clustered_pesr_vcf_list))) {
    call mosaic_pesr_part1.Mosaic as pesr1{
      input:
        name=basename(clustered_depth_vcfs[i]),
        pesr_vcfs=read_lines(per_batch_clustered_pesr_vcf_list[i]),
        metrics=agg_metrics[i],
        cutoffs=RF_cutoffs[i],
        coverage_file=coverage_files[i],
        coverage_file_idx=coverage_file_idxs[i],
        fam_file=fam_file,
        median_file=median_files[i],
        sv_pipeline_docker=sv_pipeline_docker
        
    }
  }
  scatter (i in range(length(clustered_depth_vcfs))) {
    call depth_mosaic.Mosaic as depth{
      input:
        name=basename(clustered_depth_vcfs[i]),
        metrics=agg_metrics[i],
        cutoffs=RF_cutoffs[i],
        rare_cutoff=rare_cutoff,
        depth_vcf=clustered_depth_vcfs[i],
        lookup=LookupGen.depthlookup,
        coverage_file=coverage_files[i],
        coverage_file_idx=coverage_file_idxs[i],
        fam_file=fam_file,
        median_file=median_files[i],
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker
        
    }
  }
  call preRF.make_cohort_VCFs as LookupGen {
    input:
      pesr_vcfs = pesr1.merged_pesr,
      depth_vcfs = clustered_depth_vcfs,
      sv_pipeline_docker=sv_pipeline_docker
  }
  scatter (i in range(length(pesr1.common_potential))) {
    call mosaic_pesr_part2.Mosaic as pesr2{
      input:
        name=basename(pesr1.common_potential[i]),
        outlier=outlier,
        rare_cutoff=rare_cutoff,
        lookup=LookupGen.pesrlookup,
        potential=pesr1.common_potential[i],
        coverage_file=coverage_files[i],
        coverage_file_idx=coverage_file_idxs[i],
        fam_file=fam_file,
        median_file=median_files[i],
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker
    }
  }
  output{
    Array[File] depthplots = depth.igvplots
    Array[File] pesrplots =  pesr2.igvplots
    Array[File] raredepthpotential = depth.rare_potential
    Array[File ] rarepesrpotential = pesr2.potentialmosaic
  }
}
