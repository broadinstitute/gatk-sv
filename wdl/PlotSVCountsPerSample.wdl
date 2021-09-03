version 1.0

import "Structs.wdl"

# Workflow to count SVs per sample per type and 
workflow PlotSVCountsPerSample {
  input {
    String prefix
    Array[File?] vcfs  # in order of vcf_identifiers array. To skip one, use null keyword
    Int N_IQR_cutoff
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_count_svs
    RuntimeAttr? runtime_attr_plot_svcounts
    RuntimeAttr? runtime_attr_cat_outliers_preview
  }

  scatter (vcf in vcfs) {
    if (defined(vcf)) {
      String vcf_name = basename(select_first([vcf]), ".vcf.gz")
      call CountSVsPerSamplePerType {
        input:
          vcf = select_first([vcf]),
          prefix = vcf_name,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_count_svs
      }

      call PlotSVCountsWithCutoff {
        input:
          svcounts = CountSVsPerSamplePerType.sv_counts,
          n_iqr_cutoff = N_IQR_cutoff,
          prefix = vcf_name,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_plot_svcounts
      }
    }
  }

  # Merge list of outliers for previewing
  call CatOutliersPreview {
    input:
      outliers = select_all(PlotSVCountsWithCutoff.outliers_list),
      prefix = prefix,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_cat_outliers_preview
  }


  output {
    Array[File] sv_counts = select_all(CountSVsPerSamplePerType.sv_counts)
    Array[File] sv_count_plots = select_all(PlotSVCountsWithCutoff.svcount_distrib_plots)
    File outlier_samples_preview = CatOutliersPreview.outliers_preview_file
  }
}

task CountSVsPerSamplePerType {
  input {
    File vcf
    String prefix
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
    File sv_counts = "~{prefix}.svcounts.txt"
  }
  command <<<

    set -euo pipefail
    # Count sv per class per sample
    svtk count-svtypes ~{vcf} ~{prefix}.svcounts.txt
  
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

# Generate SV count distribution plots with IQR cutoff plotted
task PlotSVCountsWithCutoff {
  input{
    File svcounts
    Int n_iqr_cutoff
    String prefix
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

  command <<<
    set -euo pipefail
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/determine_svcount_outliers.R  \
      -p "~{prefix}" \
      -I "~{n_iqr_cutoff}" \
      ~{svcounts} \
      "./"
  >>>

  output {
    File svcount_distrib_plots = "~{prefix}.all_SVTYPEs.counts_per_sample.png"
    File outliers_list = "~{prefix}.SV_count_outlier_samples.txt"
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

# Merge outlier sample lists detected across multiple VCFs
task CatOutliersPreview {
  input {
    Array[File] outliers
    String prefix
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
    File outliers_preview_file = "${prefix}.outliers_preview.samples.txt"
  }
  command <<<

    set -euo pipefail
    python3 <<CODE
    outliers_dict = {}
    for file in ["~{sep='", "' outliers}"]:
      with open(file, 'r') as IN:
        for line in IN:
          if line.strip().startswith("#sample"):
            continue
          sample, reason = line.strip().split("\t")
          if sample in outliers_dict:
            outliers_dict[sample].add(reason)
          else:
            outliers_dict[sample] = {reason}
    with open("~{prefix}.outliers_preview.samples.txt", 'w') as OUT:
      OUT.write("#sample\treason\n")
      for sample in sorted(list(outliers_dict.keys())):
        OUT.write(sample + "\t" + ",".join(outliers_dict[sample]) + "\n")
    CODE
  
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


