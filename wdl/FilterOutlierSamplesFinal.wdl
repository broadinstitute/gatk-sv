version 1.0

import "Utils.wdl" as util
import "ShardedCountSVs.wdl" as count
import "PlotSVCountsPerSample.wdl" as plot

workflow FilterOutlierSamplesFinal {
  input {
    String prefix
    Array[File] vcfs  # sharded by contig, in order, allosomes last
    File additional_samples_to_remove  # can be empty file if no additional

    Int N_IQR_cutoff
    Boolean filter_vcf
    Boolean plot_counts
    String? bcftools_preprocessing_options  # for preprocessing prior to SV counting. will not affect output VCFs

    File? sv_counts_in  # SV counts file - if not provided, will create
    File? all_samples_to_remove  # if not provided, will determine outliers to remove. if provided, will just filter vcfs
    
    File autosomes_list
    File contigs_list

    String sv_pipeline_docker
    String sv_base_mini_docker

  }

  Array[Array[String]] autosomes = read_tsv(autosomes_list)
  Array[String] contigs = read_lines(contigs_list)

  if (!defined(sv_counts_in)) {
    scatter ( i in range(length(autosomes)) ) {
      call count.ShardedCountSVs {
        input:
          vcf = vcfs[i],
          bcftools_preprocessing_options = bcftools_preprocessing_options,
          prefix = "~{prefix}.~{autosomes[i][0]}",
          sv_pipeline_docker = sv_pipeline_docker,
          sv_base_mini_docker = sv_base_mini_docker
      }
    }

    call count.SumSVCounts {
      input:
        counts = ShardedCountSVs.sv_counts,
        prefix = "~{prefix}.all_autosomes",
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  File sv_counts_ = select_first([sv_counts_in, SumSVCounts.sv_counts_summed])

  if (!defined(all_samples_to_remove) || plot_counts ) {
    call plot.PlotSVCountsWithCutoff {
      input:
        svcounts = sv_counts_,
        n_iqr_cutoff = N_IQR_cutoff,
        prefix = prefix,
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  if (!defined(all_samples_to_remove)) {
    call util.GetSampleIdsFromVcf {
      input:
        vcf = vcfs[0],
        sv_base_mini_docker = sv_base_mini_docker
    }

    call GetOutliersListAndCount {
      input:
        outlier_samples_with_reason = select_first([PlotSVCountsWithCutoff.outliers_list]),
        cohort_samples = GetSampleIdsFromVcf.out_file,
        additional_samples_to_remove = additional_samples_to_remove,
        prefix = prefix,
        sv_base_mini_docker = sv_base_mini_docker
    }
  }

  File to_remove_ = select_first([all_samples_to_remove, GetOutliersListAndCount.samples_to_remove])

  if (filter_vcf) {
    scatter (i in range(length(contigs))) {
      call util.SubsetVcfBySamplesList {
        input:
          vcf = vcfs[i],
          list_of_samples = to_remove_,
          outfile_name = "~{prefix}.~{contigs[i]}.outliers_removed.vcf.gz",
          remove_samples = true,
          sv_base_mini_docker = sv_base_mini_docker
      }
    }
  }

  
  output {
    File sv_counts = sv_counts_
    File? sv_count_plots = PlotSVCountsWithCutoff.svcount_distrib_plots
    File? outlier_samples_with_reason = PlotSVCountsWithCutoff.outliers_list
    File? outlier_samples = GetOutliersListAndCount.outliers_list
    Int? num_outliers = GetOutliersListAndCount.num_outliers
    File samples_to_remove = to_remove_
    Int? num_to_remove = GetOutliersListAndCount.num_to_remove
    File? updated_sample_list = GetOutliersListAndCount.updated_sample_list
    Int? new_sample_count = GetOutliersListAndCount.new_sample_count
    Array[File]? vcfs_samples_removed = SubsetVcfBySamplesList.vcf_subset
    Array[File]? vcf_indexes_samples_removed = SubsetVcfBySamplesList.vcf_subset_index
  }
}


task GetOutliersListAndCount {
  input {
    File outlier_samples_with_reason
    File cohort_samples
    File additional_samples_to_remove
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: 10 + ceil(size([outlier_samples_with_reason, cohort_samples, additional_samples_to_remove], "GiB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    
    # extract sample id column, remove header, deduplicate outlier samples list
    cut -f1 ~{outlier_samples_with_reason} | sed '1d' | sort | uniq > ~{prefix}.SV_count_outlier_samples.txt
    wc -l < ~{prefix}.SV_count_outlier_samples.txt > outlier_count.txt

    # merge outliers & additional samples to remove and deduplicate
    cat ~{prefix}.SV_count_outlier_samples.txt ~{additional_samples_to_remove} | sort | uniq > ~{prefix}.outliers_plus_additional_to_remove.txt
    wc -l < ~{prefix}.outliers_plus_additional_to_remove.txt > removed_count.txt

    # filter cohort sample list
    fgrep -wvf ~{prefix}.outliers_plus_additional_to_remove.txt ~{cohort_samples} > ~{prefix}.sample_list.filtered.txt || true
    wc -l < ~{prefix}.sample_list.filtered.txt > new_count.txt
  >>>

  output {
    Int num_outliers = read_int("outlier_count.txt")
    File outliers_list = "~{prefix}.SV_count_outlier_samples.txt"
    File samples_to_remove = "~{prefix}.outliers_plus_additional_to_remove.txt"
    Int num_to_remove = read_int("removed_count.txt")
    File updated_sample_list = "~{prefix}.sample_list.filtered.txt"
    Int new_sample_count = read_int("new_count.txt")
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}