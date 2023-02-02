version 1.0

import "CollectQcVcfWide.wdl" as vcfwideqc
import "CollectQcPerSample.wdl" as persample
import "CollectSiteLevelBenchmarking.wdl" as sitebench
import "CollectPerSampleBenchmarking.wdl" as samplebench
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "Utils.wdl" as Utils

# Main workflow to perform comprehensive quality control (QC) on
# an SV VCF output by GATK-SV
workflow MainVcfQc {
  input {
    Array[File] vcfs  # Option to provide a single GATK-SV VCF or an array of position-sharded SV VCFs. Must be indexed
    Boolean vcf_format_has_cn = true
    String? bcftools_preprocessing_options
    File? ped_file
    File? list_of_samples_to_include
    File? sample_renaming_tsv # File with mapping to rename per-sample benchmark sample IDs for compatibility with cohort
    Int max_trios = 1000
    String prefix
    Int sv_per_shard
    Int samples_per_shard
    Array[Array[String]]? site_level_comparison_datasets    # Array of two-element arrays, one per dataset, each of format [prefix, gs:// path to directory with one BED per population]
    Array[Array[String]]? sample_level_comparison_datasets  # Array of two-element arrays, one per dataset, each of format [prefix, gs:// path to per-sample tarballs]
    File primary_contigs_fai
    Int? random_seed
    Int? max_gq  # Max GQ for plotting. Default = 99, ie. GQ is on a scale of [0,99]. Prior to CleanVcf, use 999
    Int? downsample_qc_per_sample  # Number of samples to use for per-sample QC. Default: 1000

    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_qc_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_plot_qc_vcf_wide
    RuntimeAttr? runtime_override_site_level_benchmark_plot
    RuntimeAttr? runtime_override_plot_qc_per_sample
    RuntimeAttr? runtime_override_plot_qc_per_family
    RuntimeAttr? runtime_override_per_sample_benchmark_plot
    RuntimeAttr? runtime_override_sanitize_outputs

    # overrides for MiniTasks or Utils
    RuntimeAttr? runtime_override_subset_vcf
    RuntimeAttr? runtime_override_rename_vcf_samples
    RuntimeAttr? runtime_override_merge_vcfwide_stat_shards
    RuntimeAttr? runtime_override_merge_vcf_2_bed

    # overrides for CollectQcVcfWide
    RuntimeAttr? runtime_override_preprocess_vcf
    RuntimeAttr? runtime_override_collect_sharded_vcf_stats
    RuntimeAttr? runtime_override_svtk_vcf_2_bed
    RuntimeAttr? runtime_override_scatter_vcf
    RuntimeAttr? runtime_override_merge_subvcf_stat_shards

    # overrides for CollectSiteLevelBenchmarking
    RuntimeAttr? runtime_override_site_level_benchmark
    RuntimeAttr? runtime_override_merge_site_level_benchmark

    # overrides for CollectQcPerSample
    RuntimeAttr? runtime_override_collect_vids_per_sample
    RuntimeAttr? runtime_override_split_samples_list
    RuntimeAttr? runtime_override_tar_shard_vid_lists
    RuntimeAttr? runtime_override_merge_sharded_per_sample_vid_lists

    # overrides for CollectPerSampleBenchmarking
    RuntimeAttr? runtime_override_benchmark_samples
    RuntimeAttr? runtime_override_split_shuffled_list
    RuntimeAttr? runtime_override_merge_and_tar_shard_benchmarks
  }

  Array[String] contigs = transpose(read_tsv(primary_contigs_fai))[0]

  # Restrict to a subset of all samples, if optioned. This can be useful to 
  # exclude outlier samples, or restrict to males/females on X/Y (for example)

  if (defined(list_of_samples_to_include)) {
    scatter ( vcf in vcfs ) {
      call Utils.SubsetVcfBySamplesList {
        input:
          vcf=vcf,
          vcf_idx=vcf + ".tbi",
          list_of_samples=select_first([list_of_samples_to_include]),
          sv_base_mini_docker=sv_base_mini_docker,
          runtime_attr_override=runtime_override_subset_vcf
      }
    }
  }

  Array[File] vcfs_for_qc = select_first([SubsetVcfBySamplesList.vcf_subset, vcfs])

  # Scatter raw variant data collection per chromosome
  scatter ( contig in contigs ) {
    # Collect VCF-wide summary stats
    call vcfwideqc.CollectQcVcfWide {
      input:
        vcfs=vcfs_for_qc,
        contig=contig,
        sv_per_shard=sv_per_shard,
        bcftools_preprocessing_options=bcftools_preprocessing_options,
        prefix="~{prefix}.~{contig}",
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_override_preprocess_vcf=runtime_override_preprocess_vcf,
        runtime_override_collect_sharded_vcf_stats=runtime_override_collect_sharded_vcf_stats,
        runtime_override_svtk_vcf_2_bed=runtime_override_svtk_vcf_2_bed,
        runtime_override_scatter_vcf=runtime_override_scatter_vcf,
        runtime_override_merge_subvcf_stat_shards=runtime_override_merge_subvcf_stat_shards,
        runtime_override_merge_svtk_vcf_2_bed=runtime_override_merge_vcf_2_bed
    }
  }

  # Merge shards into single VCF stats file
  call MiniTasks.ConcatBeds as MergeVcfWideStats {
    input:
      shard_bed_files=CollectQcVcfWide.vcf_stats,
      prefix="~{prefix}.VCF_sites.stats",
      index_output=true,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_merge_vcfwide_stat_shards
  }

  # Merge vcf2bed output
  call MiniTasks.ConcatBeds as MergeVcf2Bed {
    input:
      shard_bed_files=CollectQcVcfWide.vcf2bed_out,
      prefix="~{prefix}.vcf2bed",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_merge_vcf_2_bed
  }

  # Plot VCF-wide summary stats
  call PlotQcVcfWide {
    input:
      vcf_stats=MergeVcfWideStats.merged_bed_file,
      samples_list=CollectQcVcfWide.samples_list[0],
      prefix=prefix,
      sv_pipeline_qc_docker=sv_pipeline_qc_docker,
      runtime_attr_override=runtime_override_plot_qc_vcf_wide
  }

  # Collect and plot site-level benchmarking vs. external datasets
  if (defined(site_level_comparison_datasets)) {
    scatter ( comparison_dataset_info in select_first([site_level_comparison_datasets, 
                                                       [[], []]]) ) {
      # Collect site-level external benchmarking data
      call sitebench.CollectSiteLevelBenchmarking {
        input:
          vcf_stats=MergeVcfWideStats.merged_bed_file,
          prefix=prefix,
          contigs=contigs,
          benchmark_url=comparison_dataset_info[1],
          benchmark_name=comparison_dataset_info[0],
          sv_pipeline_qc_docker=sv_pipeline_qc_docker,
          sv_base_mini_docker=sv_base_mini_docker,
          runtime_override_site_level_benchmark=runtime_override_site_level_benchmark,
          runtime_override_merge_site_level_benchmark=runtime_override_merge_site_level_benchmark
      }

      # Plot site-level benchmarking results
      call PlotSiteLevelBenchmarking {
        input:
          benchmarking_tarball=CollectSiteLevelBenchmarking.benchmarking_results_tarball,
          tarball_dir_name=CollectSiteLevelBenchmarking.tarball_dir_name,
          benchmark_name=comparison_dataset_info[0],
          sv_pipeline_qc_docker=sv_pipeline_qc_docker,
          runtime_attr_override=runtime_override_site_level_benchmark_plot
      }
    }
  }

  # Shard sample list
  call MiniTasks.SplitUncompressed as SplitSamplesList {
    input:
      whole_file=CollectQcVcfWide.samples_list[0],
      lines_per_shard=samples_per_shard,
      shard_prefix="~{prefix}.list_shard.",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_split_samples_list
  }

  # Collect per-sample VID lists for each sample shard
  scatter ( shard in SplitSamplesList.shards ) {
    call persample.CollectQcPerSample {
      input:
        vcfs=vcfs_for_qc,
        vcf_format_has_cn=vcf_format_has_cn,
        samples_list=shard,
        prefix=prefix,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_override_collect_vids_per_sample=runtime_override_collect_vids_per_sample,
        runtime_override_merge_sharded_per_sample_vid_lists=runtime_override_merge_sharded_per_sample_vid_lists
    }
  }
  
  # Merge all VID lists into single output directory and tar it
  call TarShardVidLists {
    input:
      in_tarballs=CollectQcPerSample.vid_lists,
      folder_name="~{prefix}_perSample_VIDs_merged",
      tarball_prefix="~{prefix}_perSample_VIDs",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_tar_shard_vid_lists
  }
  
  Int max_gq_ = select_first([max_gq, 99])
  # Plot per-sample stats
  call PlotQcPerSample {
    input:
      vcf_stats=MergeVcfWideStats.merged_bed_file,
      samples_list=CollectQcVcfWide.samples_list[0],
      per_sample_tarball=TarShardVidLists.vid_lists,
      prefix=prefix,
      max_gq=max_gq_,
      downsample_qc_per_sample=downsample_qc_per_sample,
      sv_pipeline_qc_docker=sv_pipeline_qc_docker,
      runtime_attr_override=runtime_override_plot_qc_per_sample
  }

  # Plot per-family stats if .ped file provided as input
  if (defined(ped_file)) {
    call PlotQcPerFamily {
      input:
        vcf_stats=MergeVcfWideStats.merged_bed_file,
        samples_list=CollectQcVcfWide.samples_list[0],
        ped_file=select_first([ped_file]),
        max_trios=max_trios,
        per_sample_tarball=TarShardVidLists.vid_lists,
        prefix=prefix,
        max_gq=max_gq_,
        sv_pipeline_qc_docker=sv_pipeline_qc_docker,
        runtime_attr_override=runtime_override_plot_qc_per_family
    }
  }

  # Collect and plot per-sample benchmarking vs. external callsets
  if (defined(sample_level_comparison_datasets)) {
    scatter ( comparison_dataset_info in select_first([sample_level_comparison_datasets, 
                                                       [[], []]]) ) {
      # Collect per-sample external benchmarking data
      call samplebench.CollectPerSampleBenchmarking {
        input:
          vcf_stats=MergeVcfWideStats.merged_bed_file,
          samples_list=CollectQcVcfWide.samples_list[0],
          per_sample_tarball=TarShardVidLists.vid_lists,
          comparison_tarball=comparison_dataset_info[1],
          sample_renaming_tsv=sample_renaming_tsv,
          prefix=prefix,
          contigs=contigs,
          comparison_set_name=comparison_dataset_info[0],
          samples_per_shard=samples_per_shard,
          random_seed=random_seed,
          sv_base_mini_docker=sv_base_mini_docker,
          sv_pipeline_docker=sv_pipeline_docker,
          sv_pipeline_qc_docker=sv_pipeline_qc_docker,
          runtime_override_benchmark_samples=runtime_override_benchmark_samples,
          runtime_override_split_shuffled_list=runtime_override_split_shuffled_list,
          runtime_override_merge_and_tar_shard_benchmarks=runtime_override_merge_and_tar_shard_benchmarks
      }

      # Plot per-sample benchmarking results
      call PlotPerSampleBenchmarking {
        input:
          per_sample_benchmarking_tarball=CollectPerSampleBenchmarking.benchmarking_results_tarball,
          samples_list=CollectQcVcfWide.samples_list[0],
          comparison_set_name=comparison_dataset_info[0],
          prefix=prefix,
          sv_pipeline_qc_docker=sv_pipeline_qc_docker,
          runtime_attr_override=runtime_override_per_sample_benchmark_plot
      }
    }
  }

  # Sanitize all outputs
  call SanitizeOutputs {
    input:
      prefix=prefix,
      samples_list=CollectQcVcfWide.samples_list[0],
      vcf_stats=MergeVcfWideStats.merged_bed_file,
      vcf_stats_idx=MergeVcfWideStats.merged_bed_idx,
      plot_qc_vcfwide_tarball=PlotQcVcfWide.plots_tarball,
      plot_qc_site_level_external_benchmarking_tarballs=PlotSiteLevelBenchmarking.tarball_wPlots,
      plot_qc_per_sample_tarball=PlotQcPerSample.perSample_plots_tarball,
      plot_qc_per_family_tarball=PlotQcPerFamily.perFamily_plots_tarball,
      cleaned_fam_file=PlotQcPerFamily.cleaned_fam_file,
      plot_qc_per_sample_external_benchmarking_tarballs=PlotPerSampleBenchmarking.perSample_plots_tarball,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_sanitize_outputs
  }

  # Final output
  output {
    File sv_vcf_qc_output = SanitizeOutputs.vcf_qc_tarball
    File vcf2bed_output = MergeVcf2Bed.merged_bed_file
  }
}


# Plot VCF-wide QC stats
task PlotQcVcfWide {
  input {
    File vcf_stats
    File samples_list
    String prefix
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: 20,
    cpu_cores: 1,
    preemptible_tries: 1,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_qc_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail
    
    # Plot VCF-wide distributions
    /opt/sv-pipeline/scripts/vcf_qc/plot_sv_vcf_distribs.R \
      -N $( cat ~{samples_list} | sort | uniq | wc -l ) \
      -S /opt/sv-pipeline/scripts/vcf_qc/SV_colors.txt \
      ~{vcf_stats} \
      plotQC_vcfwide_output/

    # Prep outputs
    tar -czvf ~{prefix}.plotQC_vcfwide_output.tar.gz \
      plotQC_vcfwide_output
  >>>

  output {
    File plots_tarball = "~{prefix}.plotQC_vcfwide_output.tar.gz"
  }
}


# Task to merge VID lists across shards
task TarShardVidLists {
  input {
    Array[File] in_tarballs
    String? folder_name
    String? tarball_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String tar_folder_name = select_first([folder_name, "merged"])
  String outfile_name = select_first([tarball_prefix, tar_folder_name]) + ".tar.gz"

  # Since the input files are often/always compressed themselves, assume compression factor for tarring is 1.0
  Float input_size = size(in_tarballs, "GB")
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(10.0 + input_size * 2.0),
    cpu_cores: 1,
    preemptible_tries: 1,
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
    # Create final output directory
    mkdir "~{tar_folder_name}"

    while read tarball_path; do
      tar -xzvf "$tarball_path" --directory ~{tar_folder_name}/
    done < ~{write_lines(in_tarballs)}

    # Compress final output directory
    tar -czvf "~{outfile_name}" "~{tar_folder_name}"
  >>>

  output {
    File vid_lists = outfile_name
  }
}


# Plot external benchmarking results
task PlotSiteLevelBenchmarking {
  input {
    File benchmarking_tarball
    String tarball_dir_name
    String benchmark_name
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: 20,
    cpu_cores: 1,
    preemptible_tries: 1,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_qc_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<  
    set -eu -o pipefail
      
    # Plot benchmarking stats
    /opt/sv-pipeline/scripts/vcf_qc/plotQC.external_benchmarking.helper.sh \
      ~{benchmarking_tarball} \
      ~{benchmark_name} \
      ~{tarball_dir_name}

    # Prep outputs
    tar -czvf ~{tarball_dir_name}.wPlots.tar.gz \
      ~{tarball_dir_name}
  >>>

  output {
    File tarball_wPlots = "~{tarball_dir_name}.wPlots.tar.gz"
  }
}


# Plot per-sample stats
task PlotQcPerSample {
  input {
    File vcf_stats
    File samples_list
    File per_sample_tarball
    String prefix
    Int max_gq
    Int? downsample_qc_per_sample
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr runtime_default = object {
    mem_gb: 7.75,
    disk_gb: 50,
    cpu_cores: 1,
    preemptible_tries: 1,
    max_retries: 1,
    boot_disk_gb: 10
  }

  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_qc_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail
    
    # Make per-sample directory
    mkdir ~{prefix}_perSample/

    # Untar per-sample VID lists
    mkdir tmp_untar/
    tar -xvzf ~{per_sample_tarball} \
      --directory tmp_untar/
    find tmp_untar/ -name "*.VIDs_genotypes.txt.gz" | while read FILE; do
      mv $FILE ~{prefix}_perSample/
    done

    # Plot per-sample distributions
    /opt/sv-pipeline/scripts/vcf_qc/plot_sv_perSample_distribs.R \
      -S /opt/sv-pipeline/scripts/vcf_qc/SV_colors.txt \
      ~{vcf_stats} \
      ~{samples_list} \
      ~{prefix}_perSample/ \
      ~{prefix}_perSample_plots/ \
      --maxgq ~{max_gq} \
      ~{"--downsample " + downsample_qc_per_sample}

    # Prepare output
    tar -czvf ~{prefix}.plotQC_perSample.tar.gz \
      ~{prefix}_perSample_plots
  >>>

  output {
    File perSample_plots_tarball = "~{prefix}.plotQC_perSample.tar.gz"
  }
}


# Plot per-family stats
task PlotQcPerFamily {
  input {
    File vcf_stats
    File samples_list
    File ped_file
    File per_sample_tarball
    Int max_trios
    Int? random_seed
    String prefix
    Int max_gq
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr_override
  }
  Int random_seed_ = select_first([random_seed, 0])
  RuntimeAttr runtime_default = object {
    mem_gb: 7.75,
    disk_gb: 50,
    cpu_cores: 1,
    preemptible_tries: 1,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_qc_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail
    
    # Clean fam file
    /opt/sv-pipeline/scripts/vcf_qc/cleanFamFile.sh \
      ~{samples_list} \
      ~{ped_file} \
      ~{prefix}.cleaned.fam
    rm ~{ped_file} ~{samples_list}

    # Only run if any families remain after cleaning
    n_fams=$( grep -Ev "^#" ~{prefix}.cleaned.fam | wc -l ) || true
    echo -e "DETECTED $n_fams FAMILIES"
    if [ $n_fams -gt 0 ]; then

      # Make per-sample directory
      mkdir ~{prefix}_perSample/

      # Untar per-sample VID lists
      mkdir tmp_untar/
      tar -xvzf ~{per_sample_tarball} \
        --directory tmp_untar/
      for FILE in $( find tmp_untar/ -name "*.VIDs_genotypes.txt.gz" ); do
        mv -v $FILE ~{prefix}_perSample/
      done

      # Subset fam file, if optioned
      n_trios=$( grep -Ev "^#" ~{prefix}.cleaned.fam \
                 | awk '{ if ($2 != "0" && $2 != "." && \
                              $3 != "0" && $3 != "." && \
                              $4 != "0" && $4 != ".") print $0 }' \
                 | wc -l )
      echo -e "DETECTED $n_trios COMPLETE TRIOS"
      if [ $n_trios -gt ~{max_trios} ]; then
        grep -E '^#' ~{prefix}.cleaned.fam > fam_header.txt
        grep -Ev "^#" ~{prefix}.cleaned.fam \
        | awk '{ if ($2 != "0" && $2 != "." && \
                     $3 != "0" && $3 != "." && \
                     $4 != "0" && $4 != ".") print $0 }' \
        | sort -R --random-source <( yes ~{random_seed_} ) \
        > cleaned.shuffled.fam
        awk -v max_trios="~{max_trios}" 'NR <= max_trios' cleaned.shuffled.fam \
        | cat fam_header.txt - \
        > cleaned.subset.fam
        echo -e "SUBSETTED TO $( cat cleaned.subset.fam | wc -l | awk '{ print $1-1 }' ) RANDOM FAMILIES"
      else
        echo -e "NUMBER OF TRIOS DETECTED ( $n_trios ) LESS THAN MAX_TRIOS ( ~{max_trios} ); PROCEEDING WITHOUT DOWNSAMPLING"
        cp ~{prefix}.cleaned.fam cleaned.subset.fam
      fi
      
      # Run family analysis
      echo -e "STARTING FAMILY-BASED ANALYSIS"
      /opt/sv-pipeline/scripts/vcf_qc/analyze_fams.R \
        -S /opt/sv-pipeline/scripts/vcf_qc/SV_colors.txt \
        ~{vcf_stats} \
        cleaned.subset.fam \
        ~{prefix}_perSample/ \
        ~{prefix}_perFamily_plots/ \
        --maxgq ~{max_gq}

    else

      mkdir ~{prefix}_perFamily_plots/

    fi

    # Prepare output
    echo -e "COMPRESSING RESULTS AS A TARBALL"
    tar -czvf ~{prefix}.plotQC_perFamily.tar.gz \
      ~{prefix}_perFamily_plots
  >>>

  output {
    File perFamily_plots_tarball = "~{prefix}.plotQC_perFamily.tar.gz"
    File cleaned_fam_file = "~{prefix}.cleaned.fam"
  }
}


# Plot per-sample benchmarking
task PlotPerSampleBenchmarking {
  input {
    File per_sample_benchmarking_tarball
    File samples_list
    String comparison_set_name
    String prefix
    String sv_pipeline_qc_docker
    Int? max_samples = 3000
    Int? random_seed
    RuntimeAttr? runtime_attr_override
  }
  Int random_seed_ = select_first([random_seed, 0])
  RuntimeAttr runtime_default = object {
    mem_gb: 7.75,
    disk_gb: 50,
    cpu_cores: 1,
    preemptible_tries: 1,
    max_retries: 1,
    boot_disk_gb: 10
  }

  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_qc_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail
    
    # Untar per-sample benchmarking results
    mkdir tmp_untar/
    tar -xvzf ~{per_sample_benchmarking_tarball} \
      --directory tmp_untar/
    mkdir results/

    # Subset to max_samples
    find tmp_untar/ -name "*.sensitivity.bed.gz" \
    | xargs -I {} basename {} | sed 's/\.sensitivity\.bed\.gz//g' \
    | sort -V \
    > all_samples.list
    n_samples_all=$( cat all_samples.list | wc -l )
    echo -e "IDENTIFIED $n_samples_all TOTAL SAMPLES"
    if [ $n_samples_all -gt ~{max_samples} ]; then
      echo -e "SUBSETTING TO ~{max_samples} SAMPLES"
      cat all_samples.list \
      | sort -R --random-source <( yes ~{random_seed_} ) \
      | awk -v max_samples=~{max_samples} '{ if (NR<=max_samples) print }' \
      > ~{prefix}.plotted_samples.list
    else
      cp all_samples.list ~{prefix}.plotted_samples.list
    fi
    
    while read ID; do
      find tmp_untar -name "$ID.*.bed.gz" | xargs -I {} mv {} results/
    done < ~{prefix}.plotted_samples.list
    
    # Plot per-sample benchmarking
    /opt/sv-pipeline/scripts/vcf_qc/plot_perSample_benchmarking.R \
      -c ~{comparison_set_name} \
      results/ \
      ~{samples_list} \
      /opt/sv-pipeline/scripts/vcf_qc/SV_colors.txt \
      ~{prefix}.~{comparison_set_name}_perSample_benchmarking_plots/

    # Prepare output
    tar -czvf ~{prefix}.~{comparison_set_name}_perSample_benchmarking_plots.tar.gz \
      ~{prefix}.~{comparison_set_name}_perSample_benchmarking_plots
  >>>

  output {
    File perSample_plots_tarball = "~{prefix}.~{comparison_set_name}_perSample_benchmarking_plots.tar.gz"
    File samples_plotted = "~{prefix}.plotted_samples.list"
  }
}


# Sanitize final output
task SanitizeOutputs {
  input {
    String prefix
    File samples_list
    File vcf_stats
    File vcf_stats_idx
    File plot_qc_vcfwide_tarball
    Array[File]? plot_qc_site_level_external_benchmarking_tarballs
    File plot_qc_per_sample_tarball
    File? plot_qc_per_family_tarball
    File? cleaned_fam_file
    Array[File]? plot_qc_per_sample_external_benchmarking_tarballs
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  # simple compress + tar workflow
  Float isize_1 = size([samples_list, vcf_stats, vcf_stats_idx, plot_qc_vcfwide_tarball, plot_qc_per_sample_tarball], "GiB")
  Float isize_2 = size(select_first([plot_qc_site_level_external_benchmarking_tarballs, []]), "GiB")
  Float isize_3 = size(select_first([plot_qc_per_family_tarball, []]), "GiB")
  Float isize_4 = size(select_first([plot_qc_per_sample_external_benchmarking_tarballs, []]), "GiB")
  Float input_size = isize_1 + isize_2 + isize_3 + isize_4
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(10.0 + input_size * 5.0),
    cpu_cores: 1,
    preemptible_tries: 1,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }
  
  command <<<
    set -euxo pipefail
    
    # Prep output directory tree
    mkdir ~{prefix}_SV_VCF_QC_output/
    mkdir ~{prefix}_SV_VCF_QC_output/data/
    mkdir ~{prefix}_SV_VCF_QC_output/plots/
    mkdir ~{prefix}_SV_VCF_QC_output/plots/main_plots/
    mkdir ~{prefix}_SV_VCF_QC_output/plots/supplementary_plots/
    mkdir ~{prefix}_SV_VCF_QC_output/plots/supplementary_plots/vcf_summary_plots/
    for tarball_fname in ~{sep=" " plot_qc_site_level_external_benchmarking_tarballs}; do
      dname="$( basename -s '.tar.gz' $tarball_fname )_site_level_benchmarking_plots/"
      mkdir ~{prefix}_SV_VCF_QC_output/plots/supplementary_plots/$dname
    done
    mkdir ~{prefix}_SV_VCF_QC_output/plots/supplementary_plots/per_sample_plots/
    if ~{defined(plot_qc_per_family_tarball)}; then
      mkdir ~{prefix}_SV_VCF_QC_output/plots/supplementary_plots/sv_inheritance_plots/
    fi
    for tarball_fname in ~{sep=" " plot_qc_per_sample_external_benchmarking_tarballs}; do
      dname="$( basename -s '.tar.gz' $tarball_fname )_per_sample_benchmarking_plots/"
      mkdir ~{prefix}_SV_VCF_QC_output/plots/supplementary_plots/$dname
    done

    # Process VCF-wide stats
    cp ~{vcf_stats} \
      ~{prefix}_SV_VCF_QC_output/data/~{prefix}.VCF_sites.stats.bed.gz
    cp ~{vcf_stats_idx} \
      ~{prefix}_SV_VCF_QC_output/data/~{prefix}.VCF_sites.stats.bed.gz.tbi

    # Process VCF-wide plots
    tar -xzvf ~{plot_qc_vcfwide_tarball}
    cp plotQC_vcfwide_output/main_plots/* \
      ~{prefix}_SV_VCF_QC_output/plots/main_plots/
    cp plotQC_vcfwide_output/supporting_plots/vcf_summary_plots/* \
      ~{prefix}_SV_VCF_QC_output/plots/supplementary_plots/vcf_summary_plots/

    # Process site-level external benchmarking plots
    if ~{defined(plot_qc_site_level_external_benchmarking_tarballs)}; then
      for tarball_fname in ~{sep=" " plot_qc_site_level_external_benchmarking_tarballs}; do
        bname="$( basename -s '.tar.gz' $tarball_fname \
                  | sed -e 's/^~{prefix}\.//g' -e 's/\.wPlots$//g' )"
        dname="$( basename -s '.tar.gz' $tarball_fname )_site_level_benchmarking_plots/"
        tar -xzvf $tarball_fname
        cp $bname/data/* \
          ~{prefix}_SV_VCF_QC_output/data/ || true
        cp $bname/plots/*.ALL/main_plots/*.callset_benchmarking.png \
          ~{prefix}_SV_VCF_QC_output/plots/main_plots/ || true
        cp -r $bname/plots/* \
          ~{prefix}_SV_VCF_QC_output/plots/supplementary_plots/$dname || true
      done
    fi

    # Process per-sample plots
    tar -xzvf ~{plot_qc_per_sample_tarball}
    cp ~{prefix}_perSample_plots/main_plots/* \
      ~{prefix}_SV_VCF_QC_output/plots/main_plots/
    cp ~{prefix}_perSample_plots/supporting_plots/per_sample_plots/* \
      ~{prefix}_SV_VCF_QC_output/plots/supplementary_plots/per_sample_plots/

    # Process per-family plots
    if ~{defined(plot_qc_per_family_tarball)}; then
      tar -xzvf ~{plot_qc_per_family_tarball}
      cp ~{prefix}_perFamily_plots/main_plots/* \
        ~{prefix}_SV_VCF_QC_output/plots/main_plots/ || true
      cp ~{prefix}_perFamily_plots/supporting_plots/sv_inheritance_plots/* \
        ~{prefix}_SV_VCF_QC_output/plots/supplementary_plots/sv_inheritance_plots/ || true
    fi

    # Process per-sample external benchmarking plots
    if ~{defined(plot_qc_per_sample_external_benchmarking_tarballs)}; then
      for tarball_fname in ~{sep=" " plot_qc_per_sample_external_benchmarking_tarballs}; do
        bname="$( basename -s '.tar.gz' $tarball_fname )" 
        dname="$bname""_per_sample_benchmarking_plots/"
        tar -xzvf $tarball_fname
        cp $bname/main_plots/* \
          ~{prefix}_SV_VCF_QC_output/plots/main_plots/ || true
        cp $bname/supporting_plots/* \
          ~{prefix}_SV_VCF_QC_output/plots/supplementary_plots/$dname || true
      done
    fi

    # Process misc files
    if ~{defined(cleaned_fam_file)}; then
      cp ~{cleaned_fam_file} \
        ~{prefix}_SV_VCF_QC_output/data/~{prefix}.cleaned_trios.fam || true
    fi
    cp ~{samples_list} \
      ~{prefix}_SV_VCF_QC_output/data/~{prefix}.samples_analyzed.list

    # Compress final output
    tar -czvf ~{prefix}_SV_VCF_QC_output.tar.gz \
      ~{prefix}_SV_VCF_QC_output
  >>>

  output {
    File vcf_qc_tarball = "~{prefix}_SV_VCF_QC_output.tar.gz"
  }
}

