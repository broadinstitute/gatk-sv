version 1.0 

import "Structs.wdl"

# This is an analysis WDL to identify & filter outliers from VCFs 
# after minGQ filtering at the end of the Talkowski SV pipeline

# Treats PCR+ and PCR- samples separately

workflow FilterOutlierSamplesPostMinGQ {
  input{
    File vcf
    File vcf_idx
    File? pcrplus_samples_list
    Int? n_iqr_cutoff_pcrplus
    Int n_iqr_cutoff_pcrminus
    String prefix
    File autosomes_list
    String sv_pipeline_docker
  }
  Array[Array[String]] contigs=read_tsv(autosomes_list)
  Boolean PCRPLUS = defined(pcrplus_samples_list)

  # Write original list of unfiltered samples and split by PCR status
  call WriteSamplesList {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      pcrplus_samples_list=pcrplus_samples_list,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker
  }

  # Get count of biallelic autosomal variants per sample
  scatter ( contig in contigs ) {
    call CountSvtypes {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        prefix=prefix,
        contig=contig[0],
        sv_pipeline_docker=sv_pipeline_docker
    }
  }
  call CombineCounts {
    input:
      svcounts=CountSvtypes.sv_counts,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker
  }

  # Get outliers
  if (PCRPLUS) {
    call IdentifyOutliers as identify_PCRPLUS_outliers {
      input:
        svcounts=CombineCounts.summed_svcounts,
        n_iqr_cutoff=select_first([n_iqr_cutoff_pcrplus]),
        samples_list=WriteSamplesList.plus_samples_list,
        prefix="~{prefix}.PCRPLUS",
        sv_pipeline_docker=sv_pipeline_docker
    }
  }
  call IdentifyOutliers as identify_PCRMINUS_outliers {
    input:
      svcounts=CombineCounts.summed_svcounts,
      n_iqr_cutoff=n_iqr_cutoff_pcrminus,
      samples_list=WriteSamplesList.minus_samples_list,
      prefix="~{prefix}.PCRMINUS",
      sv_pipeline_docker=sv_pipeline_docker
  }

  # Exclude outliers from vcf
  call ExcludeOutliers {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      plus_outliers_list=identify_PCRPLUS_outliers.outliers_list,
      minus_outliers_list=identify_PCRMINUS_outliers.outliers_list,
      outfile="~{prefix}.outliers_removed.vcf.gz",
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker
  }

  # Write new list of samples without outliers
  call FilterSampleList {
    input:
      original_samples_list=WriteSamplesList.samples_list,
      outlier_samples=ExcludeOutliers.merged_outliers_list,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker
  }

  # Final outputs
  output {
    File vcf_noOutliers = ExcludeOutliers.vcf_no_outliers
    File vcf_noOutliers_idx = ExcludeOutliers.vcf_no_outliers_idx
    File nooutliers_samples_list = FilterSampleList.filtered_samples_list
    File excluded_samples_list = ExcludeOutliers.merged_outliers_list
    File svcounts_per_sample_data = CombineCounts.summed_svcounts
    File? svcounts_per_sample_plots_PCRPLUS = identify_PCRPLUS_outliers.svcount_distrib_plots
    File svcounts_per_sample_plots_PCRMINUS = identify_PCRMINUS_outliers.svcount_distrib_plots
  }
}


# Write original list of samples
task WriteSamplesList {
  input{
    File vcf
    File vcf_idx
    File? pcrplus_samples_list
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
   RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    tabix -H ~{vcf} | fgrep -v "##" \
    | cut -f10- | sed 's/\t/\n/g' > "~{prefix}.samples.list"
    if [ ! -z "~{pcrplus_samples_list}" ];then
      fgrep -wf ~{pcrplus_samples_list} "~{prefix}.samples.list" \
      > "~{prefix}.PCRPLUS.samples.list" || true
      fgrep -wvf ~{pcrplus_samples_list} "~{prefix}.samples.list" \
      > "~{prefix}.PCRMINUS.samples.list" || true
    else
      cp ~{prefix}.samples.list "~{prefix}.PCRMINUS.samples.list"
      touch "~{prefix}.PCRPLUS.samples.list"
    fi
  >>>

  output {
    File samples_list = "~{prefix}.samples.list"
    File plus_samples_list = "~{prefix}.PCRPLUS.samples.list"
    File minus_samples_list = "~{prefix}.PCRMINUS.samples.list"
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


# Count biallelic SV per sample for a single chromosome
task CountSvtypes {
  input{
    File vcf
    File vcf_idx
    String prefix
    String contig
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    tabix --print-header "~{vcf}" "~{contig}" \
    | fgrep -v "MULTIALLELIC" \
    | fgrep -v "PESR_GT_OVERDISPERSION" \
    | svtk count-svtypes --no-header stdin \
    | awk -v OFS="\t" -v chr="~{contig}" '{ print $0, chr }' \
    > "~{prefix}.~{contig}.svcounts.txt"
  >>>

  output {
    File sv_counts = "~{prefix}.~{contig}.svcounts.txt"
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


# Combine SV count files across chromosomes
task CombineCounts {
  input{
    Array[File] svcounts
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 30,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    while read file; do
      cat "$file"
    done < ~{write_lines(svcounts)} \
    > merged_svcounts.txt
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/sum_svcounts_perSample.R \
      merged_svcounts.txt \
      "~{prefix}.summed_svcounts_per_sample.txt"
  >>>


  output {
    File summed_svcounts = "~{prefix}.summed_svcounts_per_sample.txt"
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


# Identify the list of outlier samples & generate distribution plots
task IdentifyOutliers {
  input{
    File svcounts
    Int n_iqr_cutoff
    File samples_list
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 20,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    # Subset input data to specified samples
    sed -n '1p' ~{svcounts} > filtered_counts.input.txt
    sed '1d' ~{svcounts} | fgrep -wf ~{samples_list} >> filtered_counts.input.txt
    # Return list of samples exceeding cutoff for at least one sv class
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/determine_svcount_outliers.R  \
      -p "~{prefix}" \
      -I "~{n_iqr_cutoff}" \
      filtered_counts.input.txt \
      "~{prefix}_svcount_outlier_plots/"
    cat "~{prefix}_svcount_outlier_plots/~{prefix}.SV_count_outlier_samples.txt" \
      | fgrep -v "#" | cut -f1 | sort | uniq \
      > "~{prefix}.SV_count_outliers.samples.list"
    tar -cvzf "~{prefix}_svcount_outlier_plots.tar.gz" "~{prefix}_svcount_outlier_plots/"
  >>>

  output {
    File outliers_list = "~{prefix}.SV_count_outliers.samples.list"
    File svcount_distrib_plots = "~{prefix}_svcount_outlier_plots.tar.gz"
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


# Exclude outliers from VCF
task ExcludeOutliers {
  input{
    File vcf
    File vcf_idx
    File? plus_outliers_list
    File minus_outliers_list
    String outfile
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  command <<<
    set -euo pipefail
    cat ~{plus_outliers_list} ~{minus_outliers_list} \
      | sort -Vk1,1 | uniq \
      > "~{prefix}.SV_count_outliers.samples.list" || true
    tabix -H ~{vcf} | fgrep -v "##" | \
      sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' | \
      fgrep -wf "~{prefix}.SV_count_outliers.samples.list" | cut -f2 > \
      indexes_to_exclude.txt || true
    if [ $( cat indexes_to_exclude.txt | wc -l ) -gt 0 ]; then
      zcat ~{vcf} | \
      cut --complement -f$( cat indexes_to_exclude.txt | paste -s -d, ) | \
      bgzip -c \
      > "~{prefix}.subsetted_preEmptyRemoval.vcf.gz" || true
      /opt/sv-pipeline/scripts/drop_empty_records.py \
        "~{prefix}.subsetted_preEmptyRemoval.vcf.gz" \
        stdout | \
      bgzip -c > ~{outfile} || true
    else
      cp ~{vcf} ~{outfile}
    fi
    tabix -p vcf -f "~{outfile}"
  >>>

  output {
    File merged_outliers_list = "~{prefix}.SV_count_outliers.samples.list"
    File vcf_no_outliers = "~{outfile}"
    File vcf_no_outliers_idx = "~{outfile}.tbi"
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


# Write new list of samples per prefix after outlier filtering
task FilterSampleList {
  input{
    File original_samples_list
    File outlier_samples
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
    fgrep -wvf ~{outlier_samples} ~{original_samples_list} > \
    ~{prefix}.outliers_excluded.samples.list || true
  >>>

  output {
    File filtered_samples_list = "~{prefix}.outliers_excluded.samples.list"
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
