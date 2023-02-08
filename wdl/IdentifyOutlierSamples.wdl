version 1.0

import "Structs.wdl"
import "PlotSVCountsPerSample.wdl" as plot_svcounts
import "Utils.wdl" as util
import "CollectQcVcfWide.wdl" as qc_utils
import "FilterOutlierSamplesPostMinGQ.wdl" as legacy

# Identifies a list of SV count outlier samples from one or more input VCFs

workflow IdentifyOutlierSamples {
  input {
    String name  # batch or cohort
    Array[File] vcfs
    File? sv_counts  # SV counts file from PlotSVCountsPerSample - if not provided, will create
    Int N_IQR_cutoff
    File? outlier_cutoff_table
    String? vcf_identifier  # required (enter algorithm here) if providing outlier_cutoff_table, optional otherwise to add to file prefixes (ie. as a VCF identifier)
    String? bcftools_preprocessing_options
    Boolean plot_counts = false
    Array[Pair[String, File]]? sample_subsets # if provided, will identify outliers separately within each subset. Expected format is array of pairs, where pair.left is the subset name and pair.right is a text file with all relevant sample IDs
    String sv_pipeline_docker
    String sv_base_mini_docker
    String linux_docker
    RuntimeAttr? runtime_attr_ids_from_vcf
    RuntimeAttr? runtime_override_preprocess_vcf
    RuntimeAttr? runtime_attr_identify_outliers
    RuntimeAttr? runtime_attr_cat_outliers
    RuntimeAttr? runtime_attr_subset_counts
    RuntimeAttr? runtime_attr_count_svs
    RuntimeAttr? runtime_attr_plot_svcounts
  }

  String prefix = if (defined(vcf_identifier)) then "~{name}_~{vcf_identifier}" else name

  if (!defined(sample_subsets)) {
    call util.GetSampleIdsFromVcf as GetSamplesList {
      input:
        vcf = vcfs[0],
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_ids_from_vcf
    }
  }
  Array[Pair[String, File]] subsets_to_eval = select_first([sample_subsets, [("ALL", select_first([GetSamplesList.out_file]))]])

  # Collect SV counts for each VCF in parallel unless sv_counts is provided
  if (!defined(sv_counts)) {

    # Process each VCF in parallel
    scatter ( vcf in vcfs ) {
      # Preprocess VCF with bcftools, if optioned
      if (defined(bcftools_preprocessing_options)) {
        call qc_utils.PreprocessVcf {
          input:
            vcf = vcf,
            prefix = basename(vcf, ".vcf.gz") + ".preprocessed",
            bcftools_preprocessing_options = select_first([bcftools_preprocessing_options]),
            sv_base_mini_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_override_preprocess_vcf
        }
      }
      File prepped_vcf = select_first([PreprocessVcf.outvcf, vcf])

      call plot_svcounts.CountSVsPerSamplePerType as CountPerVcf {
        input:
          vcf = prepped_vcf,
          prefix = prefix,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_count_svs
      }
    }

    # Combine counts across all VCFs
    call legacy.CombineCounts as Combine {
      input:
        svcounts = CountPerVcf.sv_counts,
        prefix = prefix,
        sv_pipeline_docker = sv_pipeline_docker
    }
  }
  File final_counts = select_first([sv_counts, Combine.summed_svcounts])

  # If a precomputed outlier table is provided, directly apply those cutoffs
  if (defined(outlier_cutoff_table)) {
    call IdentifyOutliersByCutoffTable {
      input:
        svcounts = final_counts,
        outlier_cutoff_table = select_first([outlier_cutoff_table]),
        outfile = "${prefix}_outliers.txt",
        algorithm = select_first([vcf_identifier]),
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_identify_outliers
    }
  }
  # Otherwise, dynamically determine outliers based on N_IQR
  if (!defined(outlier_cutoff_table)) {

    # Determine outliers for each sample subset, if optioned (e.g., PCR+ vs. PCR-)
    scatter ( subset_info in subsets_to_eval ) {
      call SubsetCounts {
        input:
          svcounts = final_counts,
          samples_list = subset_info.right,
          outfile = "${prefix}.${subset_info.left}.counts.tsv",
          linux_docker = linux_docker,
          runtime_attr_override = runtime_attr_subset_counts
      }
      call IdentifyOutliersByIQR {
        input:
          svcounts = SubsetCounts.counts_subset,
          N_IQR_cutoff = N_IQR_cutoff,
          outfile = "${prefix}.${subset_info.left}.IQR_outliers.txt",
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_identify_outliers
      }
      if (plot_counts) {
        call plot_svcounts.PlotSVCountsWithCutoff as PlotCountsPerSubset {
            input:
              svcounts = SubsetCounts.counts_subset,
              n_iqr_cutoff = N_IQR_cutoff,
              prefix = "${prefix}.${subset_info.left}",
              sv_pipeline_docker = sv_pipeline_docker,
              runtime_attr_override = runtime_attr_plot_svcounts
        }
      }
    }

    # Merge outlier lists across all subsets
    call CatOutliers as CatSubsets {
      input:
        outliers = select_all(IdentifyOutliersByIQR.outliers_list),
        batch = prefix,
        linux_docker = linux_docker,
        runtime_attr_override = runtime_attr_cat_outliers
    }
  }

  # Merge list of outliers
  call CatOutliers {
    input:
      outliers = select_all([CatSubsets.outliers_file, IdentifyOutliersByCutoffTable.outliers_list]),
      batch = prefix,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_cat_outliers
  }

  output {
    File outlier_samples_file = CatOutliers.outliers_file
    Array[String] outlier_samples_list = CatOutliers.outliers_list
    File sv_counts_file = final_counts
    Array[File?]? outlier_plots = PlotCountsPerSubset.svcount_distrib_plots
  }
}

task IdentifyOutliersByIQR {
  input {
    File svcounts
    Int N_IQR_cutoff
    String outfile
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File outliers_list = "${outfile}"
  }
  command <<<

    set -euo pipefail

    # Return list of samples exceeding cutoff for at least one sv class
    /opt/sv-pipeline/03_variant_filtering/scripts/get_outliers_from_svcounts.helper.R \
      ~{svcounts} \
      ~{N_IQR_cutoff} \
      ~{outfile}
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

task IdentifyOutliersByCutoffTable {
  input {
    File svcounts
    File outlier_cutoff_table
    String outfile
    String algorithm
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File outliers_list = "${outfile}"
  }
  command <<<

    set -euo pipefail

    # Return list of samples exceeding cutoff for at least one sv class
    /opt/sv-pipeline/03_variant_filtering/scripts/get_outliers_from_svcounts.helper_V2.R \
      ~{svcounts} \
      ~{outlier_cutoff_table} \
      ~{outfile} \
      ~{algorithm}

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

# Restrict a file of SV counts per sample to a subset of samples
task SubsetCounts {
  input {
    File svcounts
    File samples_list
    String outfile
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File counts_subset = "${outfile}"
  }

  command <<<

    set -euo pipefail
    head -n1 ~{svcounts} > ~{outfile}
    fgrep -wf ~{samples_list} ~{svcounts} >> ~{outfile}
  
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}


# Merge outlier sample lists across algorithms
task CatOutliers {
  input {
    Array[File] outliers
    String batch
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    Array[String] outliers_list = read_lines("${batch}.outliers.samples.list")
    File outliers_file = "${batch}.outliers.samples.list"
  }
  command <<<

    set -euo pipefail
    while read file; do
      [ -e "$file" ] || continue
      cat $file
    done < ~{write_lines(outliers)} | sort | uniq > ~{batch}.outliers.samples.list
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}
