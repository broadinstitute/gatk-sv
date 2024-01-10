##########################
## EXPERIMENTAL WORKFLOW
##########################

# To obtains list of likely mosaic variants that failed RF due to separation only

version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "PreRFCohort.wdl" as preRF 
import "MosaicDepth.wdl" as depth_mosaic
import "MosaicPesrPart1.wdl" as mosaic_pesr_part1
import "MosaicPesrPart2.wdl" as mosaic_pesr_part2

workflow MosaicManualCheck{
  input{
    File fam_file
    Int rare_cutoff
    File outlier
    String prefix

    Array[File] per_batch_clustered_pesr_vcf_list # preRF 
    Array[File] clustered_depth_vcfs
    Array[File] coverage_files
    Array[File] coverage_file_idxs
    Array[File] median_files
    
    Array[File] agg_metrics
    Array[File] RF_cutoffs
    
    String sv_pipeline_docker
    String sv_base_mini_docker
    String sv_pipeline_rdtest_docker

    RuntimeAttr? runtime_attr_concat_depth_bed
    RuntimeAttr? runtime_attr_concat_pesr_bed
    RuntimeAttr? runtime_attr_concat_depth_plot
    RuntimeAttr? runtime_attr_concat_pesr_plot
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

  call MiniTasks.ConcatBeds as concat_depth_bed{
    input:
      shard_bed_files = depth.rare_potential,
      prefix = "~{prefix}.depth",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_depth_bed
  }

  call MiniTasks.ConcatBeds as concat_pesr_bed{
    input:
      shard_bed_files = pesr2.potentialmosaic,
      prefix = "~{prefix}.pesr",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_pesr_bed
  }

  call ConcatPlots as concat_depth_plots{
    input:
      shard_plots = depth.igvplots,
      prefix = "~{prefix}.depth",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_depth_plot
  }

  call ConcatPlots as concat_pesr_plots{
    input:
      shard_plots = pesr2.igvplots,
      prefix = "~{prefix}.pesr",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_pesr_plot
  }

  output{
    File depth_bed = concat_depth_bed.merged_bed_file
    File pesr_bed = concat_pesr_bed.merged_bed_file
    File depth_plots = concat_depth_plots.merged_plots
    File pesr_plots = concat_pesr_plots.merged_plots
  }
}


# Merge plots from each shard:
task ConcatPlots {
  input{
    Array[File] shard_plots
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_file="~{prefix}.plot.tar.gz"

  Float input_size = size(shard_plots, "GB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + compression_factor)),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }

  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu

    # note head -n1 stops reading early and sends SIGPIPE to zcat,
    # so setting pipefail here would result in early termination
    mkdir output_folder/

    # no more early stopping
    set -o pipefail

    while read SPLIT; do
      tar zxvf $SPLIT
      mv plots/* output_folder/
    done < ~{write_lines(shard_plots)} \

    tar zcvf ~{prefix}.plots.tar.gz output_folder/
  >>>

  output {
    File merged_plots = "~{prefix}.plots.tar.gz"
  }
}

