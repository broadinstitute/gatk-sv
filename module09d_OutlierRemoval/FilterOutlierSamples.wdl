version 1.0

import "Structs.wdl"
import "Utils.wdl" as util
import "IdentifyOutlierSamplesNewTest.wdl" as identify_outliers

# Filter outlier samples by IQR or cutoff table for one or more VCFs
# Recommended to run PlotSVCountsPerSample first to choose cutoff
workflow FilterOutlierSamples {
  input {
    String name  # batch or cohort
    Array[File] vcfs
    Array[File] vcf_idxes
    Array[String] contigs
    File? ctx_vcf
    File? ctx_vcf_idx

    Boolean split_vcf_sites = true

    File? sv_counts  # SV counts file from PlotSVCountsPerSample - if not provided, will create
    Int N_IQR_cutoff
    File? outlier_cutoff_table
    String? vcf_identifier  # required (enter algorithm here) if providing outlier_cutoff_table, otherwise used in some file prefixes
    String? bcftools_preprocessing_options
    Boolean plot_counts = false

    Array[String]? sample_subset_prefixes # if provided, will identify outliers separately within each subset
    Array[String]? sample_subset_lists # if provided, will identify outliers separately within each subset
    Int samples_per_shard = 5000
    String sv_pipeline_docker
    String sv_base_mini_docker
    String linux_docker
    RuntimeAttr? runtime_override_preprocess_vcf
    RuntimeAttr? runtime_attr_identify_outliers
    RuntimeAttr? runtime_attr_subset_vcf
    RuntimeAttr? runtime_attr_cat_outliers
    RuntimeAttr? runtime_attr_subset_counts
    RuntimeAttr? runtime_attr_filter_samples
    RuntimeAttr? runtime_attr_ids_from_vcf
    RuntimeAttr? runtime_attr_count_svs
    RuntimeAttr? runtime_attr_combine_counts
    RuntimeAttr? runtime_attr_plot_svcounts
    RuntimeAttr? runtime_attr_combine_ctx
    RuntimeAttr? runtime_attr_split_vcf_sites
  }

  if (defined(sample_subset_prefixes) && defined(sample_subset_lists)) {
    Array[Pair[String, File]]? sample_subsets = zip(select_first([sample_subset_prefixes]), 
                                                    select_first([sample_subset_lists]))
  }

  call identify_outliers.IdentifyOutlierSamples {
    input:
      vcfs = vcfs,
      name = name,
      sv_counts = sv_counts, 
      N_IQR_cutoff = N_IQR_cutoff,
      outlier_cutoff_table = outlier_cutoff_table,
      vcf_identifier = vcf_identifier,
      bcftools_preprocessing_options = bcftools_preprocessing_options,
      plot_counts = plot_counts,
      sample_subsets = sample_subsets,
      samples_per_shard = samples_per_shard,
      sv_pipeline_docker = sv_pipeline_docker,
      sv_base_mini_docker = sv_base_mini_docker,
      linux_docker = linux_docker,
      runtime_attr_ids_from_vcf = runtime_attr_ids_from_vcf,
      runtime_override_preprocess_vcf = runtime_override_preprocess_vcf,
      runtime_attr_identify_outliers = runtime_attr_identify_outliers,
      runtime_attr_cat_outliers = runtime_attr_cat_outliers,
      runtime_attr_subset_counts = runtime_attr_subset_counts,
      runtime_attr_count_svs = runtime_attr_count_svs,
      runtime_attr_combine_counts = runtime_attr_combine_counts,
      runtime_attr_plot_svcounts = runtime_attr_plot_svcounts
  }

  scatter (i in range(length(vcfs))){
    call util.SubsetVcfBySamplesList {
      input:
        vcf = vcfs[i],
        list_of_samples = IdentifyOutlierSamples.outlier_samples_file,
        outfile_name = basename(vcfs[i], ".vcf.gz") + ".${name}.outliers_removed.vcf.gz",
        remove_samples = true,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_subset_vcf
    }

    if (defined(ctx_vcf)){
          call CombineWithCtx{
              input:
                  vcf = SubsetVcfBySamplesList.vcf_subset,
                  vcf_idx = SubsetVcfBySamplesList.vcf_subset_index,
                  contig = contigs[i],
                  ctx_vcf = ctx_vcf,
                  ctx_vcf_idx = ctx_vcf_idx,
                  sv_pipeline_docker = sv_pipeline_docker,
                  runtime_attr_override = runtime_attr_combine_ctx
          }
      }

    File output_step1_vcf = select_first([CombineWithCtx.combined_vcf, SubsetVcfBySamplesList.vcf_subset])
    File output_step1_vcf_idx = select_first([CombineWithCtx.combined_vcf_idx, SubsetVcfBySamplesList.vcf_subset_index])

    if (split_vcf_sites){
        call SplitVcfSites{
            input:
              vcf = output_step1_vcf,
              vcf_idx = output_step1_vcf_idx,
              sv_base_mini_docker = sv_base_mini_docker,
              runtime_attr_override = runtime_attr_split_vcf_sites
      }
    }

  }
  
  call util.GetSampleIdsFromVcf {
    input:
      vcf = vcfs[0],
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_ids_from_vcf
  }

  # Write new list of samples without outliers
  call FilterSampleList {
    input:
      original_samples = GetSampleIdsFromVcf.out_array,
      outlier_samples = IdentifyOutlierSamples.outlier_samples_list,
      batch = name,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_filter_samples
  }

  output {
    Array[File] outlier_filtered_vcfs = output_step1_vcf
    Array[File] outlier_filtered_vcf_idxs = output_step1_vcf_idx
    Array[File?] outlier_filtered_sites = SplitVcfSites.vcf_site_vcf
    Array[File?] outlier_filtered_sites_idxs = SplitVcfSites.vcf_site_vcf_idx
    Array[String] filtered_samples_list = FilterSampleList.filtered_samples_list
    File filtered_samples_file = FilterSampleList.filtered_samples_file
    Array[String] outlier_samples_excluded = IdentifyOutlierSamples.outlier_samples_list
    File outlier_samples_excluded_file = IdentifyOutlierSamples.outlier_samples_file
    File sv_counts_file = IdentifyOutlierSamples.sv_counts_file
  }
}


# Write new list of samples per batch after outlier filtering
task FilterSampleList {
  input {
    Array[String] original_samples
    Array[String] outlier_samples
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
    Array[String] filtered_samples_list = read_lines("${batch}.outliers_excluded.samples.list")
    File filtered_samples_file = "${batch}.outliers_excluded.samples.list"
  }
  command <<<

    fgrep -wvf ~{write_lines(outlier_samples)} ~{write_lines(original_samples)} > ~{batch}.outliers_excluded.samples.list
  
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

task CombineWithCtx{
    input{
        File vcf
        File vcf_idx
        File? ctx_vcf
        File? ctx_vcf_idx
        String contig
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 200,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(vcf,".vcf.gz")

    command<<<
        set -euo pipefail

        bcftools view ~{ctx_vcf} ~{contig} -O z -o ~{contig}.ctx.vcf.gz
        tabix -p vcf ~{contig}.ctx.vcf.gz

        bcftools view -h ~{vcf} | awk '{if ($1=="#CHROM") print}' | cut -f10- | sed -e 's/\t/\n/g' > sample_list
        bcftools view -S sample_list ~{contig}.ctx.vcf.gz -O z -o ~{contig}.ctx.outlier_removed.vcf.gz
        tabix -p vcf ~{contig}.ctx.outlier_removed.vcf.gz

        echo "~{vcf}" > vcf_list
        echo "~{contig}.ctx.outlier_removed.vcf.gz" >> vcf_list

        bcftools concat -a --allow-overlaps --output-type z --file-list vcf_list --output  "~{prefix}.with_ctx.vcf.gz"
        tabix -p vcf  "~{prefix}.with_ctx.vcf.gz"

    >>>

    output{
        File out_vcf_list = "vcf_list"
        File combined_vcf = "~{prefix}.with_ctx.vcf.gz"
        File combined_vcf_idx = "~{prefix}.with_ctx.vcf.gz.tbi"
    }

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

task SplitVcfSites{

    input{
        File vcf
        File vcf_idx
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 200,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(vcf,".vcf.gz")

    command<<<
        set -euo pipefail
        zcat ~{vcf} | cut -f1-8 | bgzip > ~{prefix}.sites.gz
        tabix -p vcf ~{prefix}.sites.gz
    >>>

    output{
        File vcf_site_vcf = "~{prefix}.sites.gz"
        File vcf_site_vcf_idx = "~{prefix}.sites.gz.tbi"
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




