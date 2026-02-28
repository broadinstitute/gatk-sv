version 1.0

import "CollectQcVcfWide.wdl" as vcfwideqc
import "CollectQcPerSample.wdl" as persample
import "CollectSiteLevelBenchmarking.wdl" as sitebench
import "CollectPerSampleBenchmarking.wdl" as samplebench
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "Utils.wdl" as Utils

workflow QcAnnotations {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs

    Boolean vcf_format_has_cn = true
    File? ped_file
    File? sample_renaming_tsv
    Int max_trios = 1000

    String prefix
    Int sv_per_shard
    Int samples_per_shard
    Boolean do_per_sample_qc = true
    Array[Array[String]]? site_level_comparison_datasets
    Array[Array[String]]? sample_level_comparison_datasets

    Int? random_seed
    Int? max_gq
    Int? downsample_qc_per_sample

    File primary_contigs_fai
    File ref_fa
    File ref_fai

    String sv_base_mini_docker
    String gatk_sv_lr_docker

    RuntimeAttr? runtime_override_plot_qc_vcf_wide
    RuntimeAttr? runtime_override_site_level_benchmark_plot
    RuntimeAttr? runtime_override_plot_qc_per_sample
    RuntimeAttr? runtime_override_plot_qc_per_family
    RuntimeAttr? runtime_override_per_sample_benchmark_plot
    RuntimeAttr? runtime_override_identify_duplicates
    RuntimeAttr? runtime_override_merge_duplicates
    RuntimeAttr? runtime_override_sanitize_outputs
    RuntimeAttr? runtime_override_subset_vcf
    RuntimeAttr? runtime_override_rename_vcf_samples
    RuntimeAttr? runtime_override_merge_vcfwide_stat_shards
    RuntimeAttr? runtime_override_merge_vcf_2_bed
    RuntimeAttr? runtime_override_preprocess_vcf
    RuntimeAttr? runtime_override_collect_sharded_vcf_stats
    RuntimeAttr? runtime_override_svtk_vcf_2_bed
    RuntimeAttr? runtime_override_scatter_vcf
    RuntimeAttr? runtime_override_merge_subvcf_stat_shards
    RuntimeAttr? runtime_override_site_level_benchmark
    RuntimeAttr? runtime_override_merge_site_level_benchmark
    RuntimeAttr? runtime_override_collect_vids_per_sample
    RuntimeAttr? runtime_override_split_samples_list
    RuntimeAttr? runtime_override_tar_shard_vid_lists
    RuntimeAttr? runtime_override_merge_sharded_per_sample_vid_lists
    RuntimeAttr? runtime_override_benchmark_samples
    RuntimeAttr? runtime_override_split_shuffled_list
    RuntimeAttr? runtime_override_merge_and_tar_shard_benchmarks
  }

  Array[String] contigs = transpose(read_tsv(primary_contigs_fai))[0]

  scatter (contig in contigs) {
    call vcfwideqc.CollectQcVcfWide {
      input:
        vcfs=vcfs,
        contig=contig,
        sv_per_shard=sv_per_shard,
        ref_fa=ref_fa,
        ref_fai=ref_fai,
        prefix="~{prefix}.~{contig}",
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=gatk_sv_lr_docker,
        runtime_override_preprocess_vcf=runtime_override_preprocess_vcf,
        runtime_override_collect_sharded_vcf_stats=runtime_override_collect_sharded_vcf_stats,
        runtime_override_svtk_vcf_2_bed=runtime_override_svtk_vcf_2_bed,
        runtime_override_scatter_vcf=runtime_override_scatter_vcf,
        runtime_override_merge_subvcf_stat_shards=runtime_override_merge_subvcf_stat_shards,
        runtime_override_merge_svtk_vcf_2_bed=runtime_override_merge_vcf_2_bed
    }
  }

  call MiniTasks.ConcatBeds as MergeVcfWideStats {
    input:
      shard_bed_files=CollectQcVcfWide.vcf_stats,
      prefix="~{prefix}.VCF_sites.stats",
      index_output=true,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_merge_vcfwide_stat_shards
  }

  call MiniTasks.ConcatBeds as MergeVcf2Bed {
    input:
      shard_bed_files=CollectQcVcfWide.vcf2bed_out,
      prefix="~{prefix}.vcf2bed",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_merge_vcf_2_bed
  }

  call PlotQcVcfWide {
    input:
      vcf_stats=MergeVcfWideStats.merged_bed_file,
      samples_list=CollectQcVcfWide.samples_list[0],
      prefix=prefix,
      sv_pipeline_qc_docker=gatk_sv_lr_docker,
      runtime_attr_override=runtime_override_plot_qc_vcf_wide
  }

  if (defined(site_level_comparison_datasets)) {
    scatter (comparison_dataset_info in select_first([site_level_comparison_datasets, [[], []]]) ) {
      call sitebench.CollectSiteLevelBenchmarking {
        input:
          vcf_stats=MergeVcfWideStats.merged_bed_file,
          prefix=prefix,
          contigs=contigs,
          benchmark_url=comparison_dataset_info[1],
          benchmark_name=comparison_dataset_info[0],
          sv_pipeline_qc_docker=gatk_sv_lr_docker,
          sv_base_mini_docker=sv_base_mini_docker,
          runtime_override_site_level_benchmark=runtime_override_site_level_benchmark,
          runtime_override_merge_site_level_benchmark=runtime_override_merge_site_level_benchmark
      }

      call PlotSiteLevelBenchmarking {
        input:
          benchmarking_tarball=CollectSiteLevelBenchmarking.benchmarking_results_tarball,
          tarball_dir_name=CollectSiteLevelBenchmarking.tarball_dir_name,
          benchmark_name=comparison_dataset_info[0],
          sv_pipeline_qc_docker=gatk_sv_lr_docker,
          runtime_attr_override=runtime_override_site_level_benchmark_plot
      }
    }
  }

  if (do_per_sample_qc) {
    call MiniTasks.SplitUncompressed as SplitSamplesList {
      input:
        whole_file=CollectQcVcfWide.samples_list[0],
        lines_per_shard=samples_per_shard,
        shard_prefix="~{prefix}.list_shard.",
        sv_pipeline_docker=gatk_sv_lr_docker,
        runtime_attr_override=runtime_override_split_samples_list
    }

    scatter ( shard in SplitSamplesList.shards ) {
      call persample.CollectQcPerSample {
        input:
          vcfs=vcfs,
          vcf_format_has_cn=vcf_format_has_cn,
          samples_list=shard,
          prefix=prefix,
          sv_base_mini_docker=sv_base_mini_docker,
          sv_pipeline_docker=gatk_sv_lr_docker,
          runtime_override_collect_vids_per_sample=runtime_override_collect_vids_per_sample,
          runtime_override_merge_sharded_per_sample_vid_lists=runtime_override_merge_sharded_per_sample_vid_lists
      }
    }

    call TarShardVidLists {
      input:
        in_tarballs=CollectQcPerSample.vid_lists,
        folder_name="~{prefix}_perSample_VIDs_merged",
        tarball_prefix="~{prefix}_perSample_VIDs",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_tar_shard_vid_lists
    }

    Int max_gq_ = select_first([max_gq, 99])
    call PlotQcPerSample {
      input:
        vcf_stats=MergeVcfWideStats.merged_bed_file,
        samples_list=CollectQcVcfWide.samples_list[0],
        per_sample_tarball=TarShardVidLists.vid_lists,
        prefix=prefix,
        max_gq=max_gq_,
        downsample_qc_per_sample=downsample_qc_per_sample,
        sv_pipeline_qc_docker=gatk_sv_lr_docker,
        runtime_attr_override=runtime_override_plot_qc_per_sample
    }

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
          sv_pipeline_qc_docker=gatk_sv_lr_docker,
          runtime_attr_override=runtime_override_plot_qc_per_family
      }
    }

    if (defined(sample_level_comparison_datasets)) {
      scatter ( comparison_dataset_info in select_first([sample_level_comparison_datasets, [[], []]]) ) {
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
            sv_pipeline_docker=gatk_sv_lr_docker,
            sv_pipeline_qc_docker=gatk_sv_lr_docker,
            runtime_override_benchmark_samples=runtime_override_benchmark_samples,
            runtime_override_split_shuffled_list=runtime_override_split_shuffled_list,
            runtime_override_merge_and_tar_shard_benchmarks=runtime_override_merge_and_tar_shard_benchmarks
        }

        call PlotPerSampleBenchmarking {
          input:
            per_sample_benchmarking_tarball=CollectPerSampleBenchmarking.benchmarking_results_tarball,
            samples_list=CollectQcVcfWide.samples_list[0],
            comparison_set_name=comparison_dataset_info[0],
            prefix=prefix,
            sv_pipeline_qc_docker=gatk_sv_lr_docker,
            runtime_attr_override=runtime_override_per_sample_benchmark_plot
        }
      }
    }
  }

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

  output {
    File sv_vcf_qc_output = SanitizeOutputs.vcf_qc_tarball
    File vcf2bed_output = MergeVcf2Bed.merged_bed_file
  }
}

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
    set -euo pipefail
    
    /opt/sv-pipeline/scripts/vcf_qc/plot_sv_vcf_distribs.R \
      -N $( cat ~{samples_list} | sort | uniq | wc -l ) \
      -S /opt/sv-pipeline/scripts/vcf_qc/SV_colors.txt \
      ~{vcf_stats} \
      plotQC_vcfwide_output/

    tar -czvf ~{prefix}.plotQC_vcfwide_output.tar.gz \
      plotQC_vcfwide_output
  >>>

  output {
    File plots_tarball = "~{prefix}.plotQC_vcfwide_output.tar.gz"
  }
}

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
    set -euo pipefail

    mkdir "~{tar_folder_name}"

    while read tarball_path; do
      tar -xzvf "$tarball_path" --directory ~{tar_folder_name}/
    done < ~{write_lines(in_tarballs)}

    tar -czvf "~{outfile_name}" "~{tar_folder_name}"
  >>>

  output {
    File vid_lists = outfile_name
  }
}

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
    set -euo pipefail
      
    /opt/sv-pipeline/scripts/vcf_qc/plotQC.external_benchmarking.helper.sh \
      ~{benchmarking_tarball} \
      ~{benchmark_name} \
      ~{tarball_dir_name}

    tar -czvf ~{tarball_dir_name}.wPlots.tar.gz \
      ~{tarball_dir_name}
  >>>

  output {
    File tarball_wPlots = "~{tarball_dir_name}.wPlots.tar.gz"
  }
}

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
    set -euo pipefail
    
    mkdir ~{prefix}_perSample/

    mkdir tmp_untar/
    tar -xvzf ~{per_sample_tarball} \
      --directory tmp_untar/
    find tmp_untar/ -name "*.VIDs_genotypes.txt.gz" | while read FILE; do
      mv $FILE ~{prefix}_perSample/
    done

    /opt/sv-pipeline/scripts/vcf_qc/plot_sv_perSample_distribs.R \
      -S /opt/sv-pipeline/scripts/vcf_qc/SV_colors.txt \
      ~{vcf_stats} \
      ~{samples_list} \
      ~{prefix}_perSample/ \
      ~{prefix}_perSample_plots/ \
      --maxgq ~{max_gq} \
      ~{"--downsample " + downsample_qc_per_sample}

    tar -czvf ~{prefix}.plotQC_perSample.tar.gz \
      ~{prefix}_perSample_plots
  >>>

  output {
    File perSample_plots_tarball = "~{prefix}.plotQC_perSample.tar.gz"
  }
}

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
    mem_gb: 15,
    disk_gb: 100,
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
    set -euo pipefail
    
    /opt/sv-pipeline/scripts/vcf_qc/cleanFamFile.sh \
      ~{samples_list} \
      ~{ped_file} \
      ~{prefix}.cleaned.fam
    rm ~{ped_file} ~{samples_list}

    n_fams=$( grep -Ev "^#" ~{prefix}.cleaned.fam | wc -l ) || true
    echo -e "DETECTED $n_fams FAMILIES"
    if [ $n_fams -gt 0 ]; then
      mkdir ~{prefix}_perSample/

      mkdir tmp_untar/
      tar -xvzf ~{per_sample_tarball} \
        --directory tmp_untar/
      for FILE in $( find tmp_untar/ -name "*.VIDs_genotypes.txt.gz" ); do
        mv -v $FILE ~{prefix}_perSample/
      done

      n_trios=$( grep -Ev "^#" ~{prefix}.cleaned.fam \
                 | awk '{ if ($2 != "0" && $2 != "." && \
                              $3 != "0" && $3 != "." && \
                              $4 != "0" && $4 != ".") print $0 }' \
                 | wc -l )
      echo -e "DETECTED $n_trios COMPLETE TRIOS"
      
      cp ~{prefix}.cleaned.fam cleaned.subset.fam

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

    echo -e "COMPRESSING RESULTS AS A TARBALL"
    tar -czvf ~{prefix}.plotQC_perFamily.tar.gz \
      ~{prefix}_perFamily_plots
  >>>

  output {
    File perFamily_plots_tarball = "~{prefix}.plotQC_perFamily.tar.gz"
    File cleaned_fam_file = "~{prefix}.cleaned.fam"
  }
}

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
    set -euo pipefail
    
    mkdir tmp_untar/
    tar -xvzf ~{per_sample_benchmarking_tarball} \
      --directory tmp_untar/
    mkdir results/

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
    
    /opt/sv-pipeline/scripts/vcf_qc/plot_perSample_benchmarking.R \
      -c ~{comparison_set_name} \
      results/ \
      ~{samples_list} \
      /opt/sv-pipeline/scripts/vcf_qc/SV_colors.txt \
      ~{prefix}.~{comparison_set_name}_perSample_benchmarking_plots/

    tar -czvf ~{prefix}.~{comparison_set_name}_perSample_benchmarking_plots.tar.gz \
      ~{prefix}.~{comparison_set_name}_perSample_benchmarking_plots
  >>>

  output {
    File perSample_plots_tarball = "~{prefix}.~{comparison_set_name}_perSample_benchmarking_plots.tar.gz"
    File samples_plotted = "~{prefix}.plotted_samples.list"
  }
}

task SanitizeOutputs {
  input {
    String prefix
    File samples_list
    File vcf_stats
    File vcf_stats_idx
    File plot_qc_vcfwide_tarball
    Array[File]? plot_qc_site_level_external_benchmarking_tarballs
    File? plot_qc_per_sample_tarball
    File? plot_qc_per_family_tarball
    File? cleaned_fam_file
    Array[File]? plot_qc_per_sample_external_benchmarking_tarballs
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float isize_1 = size([samples_list, vcf_stats, vcf_stats_idx, plot_qc_vcfwide_tarball], "GiB")
  Float isize_2 = size(select_first([plot_qc_site_level_external_benchmarking_tarballs, []]), "GiB")
  Float isize_3 = size(select_first([plot_qc_per_family_tarball, []]), "GiB")
  Float isize_4 = size(select_first([plot_qc_per_sample_external_benchmarking_tarballs, []]), "GiB")
  Float isize_5 = size(select_first([plot_qc_per_sample_tarball, []]), "GiB")
  Float input_size = isize_1 + isize_2 + isize_3 + isize_4 + isize_5
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
    set -euo pipefail
    
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
    if ~{defined(plot_qc_per_sample_tarball)}; then
      mkdir ~{prefix}_SV_VCF_QC_output/plots/supplementary_plots/per_sample_plots/
    fi
    if ~{defined(plot_qc_per_family_tarball)}; then
      mkdir ~{prefix}_SV_VCF_QC_output/plots/supplementary_plots/sv_inheritance_plots/
    fi
    for tarball_fname in ~{sep=" " plot_qc_per_sample_external_benchmarking_tarballs}; do
      dname="$( basename -s '.tar.gz' $tarball_fname )_per_sample_benchmarking_plots/"
      mkdir ~{prefix}_SV_VCF_QC_output/plots/supplementary_plots/$dname
    done

    cp ~{vcf_stats} ~{prefix}_SV_VCF_QC_output/data/~{prefix}.VCF_sites.stats.bed.gz
    cp ~{vcf_stats_idx} ~{prefix}_SV_VCF_QC_output/data/~{prefix}.VCF_sites.stats.bed.gz.tbi

    tar -xzvf ~{plot_qc_vcfwide_tarball}
    cp plotQC_vcfwide_output/main_plots/* \
      ~{prefix}_SV_VCF_QC_output/plots/main_plots/
    cp plotQC_vcfwide_output/supporting_plots/vcf_summary_plots/* \
      ~{prefix}_SV_VCF_QC_output/plots/supplementary_plots/vcf_summary_plots/

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

    if ~{defined(plot_qc_per_sample_tarball)}; then
      tar -xzvf ~{plot_qc_per_sample_tarball}
      cp ~{prefix}_perSample_plots/main_plots/* \
        ~{prefix}_SV_VCF_QC_output/plots/main_plots/
      cp ~{prefix}_perSample_plots/supporting_plots/per_sample_plots/* \
        ~{prefix}_SV_VCF_QC_output/plots/supplementary_plots/per_sample_plots/
    fi

    if ~{defined(plot_qc_per_family_tarball)}; then
      tar -xzvf ~{plot_qc_per_family_tarball}
      cp ~{prefix}_perFamily_plots/main_plots/* \
        ~{prefix}_SV_VCF_QC_output/plots/main_plots/ || true
      cp ~{prefix}_perFamily_plots/supporting_plots/sv_inheritance_plots/* \
        ~{prefix}_SV_VCF_QC_output/plots/supplementary_plots/sv_inheritance_plots/ || true
    fi

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

    if ~{defined(cleaned_fam_file)}; then
      cp ~{cleaned_fam_file} \
        ~{prefix}_SV_VCF_QC_output/data/~{prefix}.cleaned_trios.fam || true
    fi
    cp ~{samples_list} \
      ~{prefix}_SV_VCF_QC_output/data/~{prefix}.samples_analyzed.list

    tar -czvf ~{prefix}_SV_VCF_QC_output.tar.gz \
      ~{prefix}_SV_VCF_QC_output
  >>>

  output {
    File vcf_qc_tarball = "~{prefix}_SV_VCF_QC_output.tar.gz"
  }
}
