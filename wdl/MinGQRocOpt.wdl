version 1.0

import "Structs.wdl"

workflow MinGQRocOpt {
  input{
    File trio_tarball
    String prefix
    File trios_list
    File conditions_table
    Int maxSVperTrio
    Float roc_max_fdr
    Int roc_min_gq
    Int roc_max_gq
    Int roc_step_gq
    Int min_sv_per_proband_per_condition
    String sv_pipeline_docker
    String sv_base_mini_docker
  }
  # Scatter over each condition and send the trio data for ROC optimization
  Array[Array[String]] conditions = read_tsv(conditions_table)
  scatter ( condition in conditions ) {
    # Subset variants to condition of interest & merge across trios
    # Also computes median & Q1/Q3 variants per proband
    # If median > min_sv_per_proband_per_condition, also runs ROC
    call FilterMergeVariantsWithROC as roc_single {
      input:
        trio_tarball=trio_tarball,
        prefix="~{prefix}",
        sv_pipeline_docker=sv_pipeline_docker,
        trios_list=trios_list,
        condition_id=condition[0],
        minSVLEN=condition[1],
        maxSVLEN=condition[2],
        minAF=condition[3],
        maxAF=condition[4],
        includeSVTYPE=condition[5],
        excludeSVTYPE=condition[6],
        includeFILTER=condition[7],
        excludeFILTER=condition[8],
        includeEV=condition[9],
        excludeEV=condition[10],
        maxSVperTrio=maxSVperTrio,
        roc_max_fdr=roc_max_fdr,
        roc_min_gq=roc_min_gq,
        roc_max_gq=roc_max_gq,
       roc_step_gq=roc_step_gq,
        min_sv_per_proband_per_condition=min_sv_per_proband_per_condition
    }
  }

  # Merge across conditions
  call CombineRocOptResults as combine {
  input:
    sv_base_mini_docker=sv_base_mini_docker,
    condition_optimizations=roc_single.roc_optimal,
    condition_distrib_stats=roc_single.distrib_stats,
    prefix="~{prefix}"
  }

  # Outputs
  output {
    File roc_optimal_merged = combine.combined_optimizations
    File distrib_stats_merged = combine.combined_distrib_stats
  }
}


# Subset variants to meet a given set of conditions, merge across trios, 
# and run ROC if condition has enough variants per sample
task FilterMergeVariantsWithROC {
  input{
    File trio_tarball
    String prefix
    File trios_list
    String condition_id
    String minSVLEN
    String maxSVLEN
    String minAF
    String maxAF
    String includeSVTYPE
    String excludeSVTYPE
    String includeFILTER
    String excludeFILTER
    String includeEV
    String excludeEV
    Int maxSVperTrio
    Float roc_max_fdr
    Int roc_min_gq
    Int roc_max_gq
    Int roc_step_gq
    Int min_sv_per_proband_per_condition
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 8,
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    tar -xzvf ~{trio_tarball}
    find -name "trio_variant_info.txt.gz" > trio_dat_list.txt
    #Iterate over families and process them one at a time
    while read famdat; do
      #Subset to variants matching condition
      /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/subset_minGQ_trio_data.R \
        --min.size "~{minSVLEN}" \
        --max.size "~{maxSVLEN}" \
        --min.freq "~{minAF}" \
        --max.freq "~{maxAF}" \
        --svtype.include "~{includeSVTYPE}" \
        --svtype.exclude "~{excludeSVTYPE}" \
        --filter.include "~{includeFILTER}" \
        --filter.exclude "~{excludeFILTER}" \
        --ev.include "~{includeEV}" \
        --ev.exclude "~{excludeEV}" \
        --max.variants "~{maxSVperTrio}" \
        "$famdat" /dev/stdout
    done < trio_dat_list.txt \
    | gzip -c \
    > "~{prefix}.~{condition_id}.trio_variants.txt.gz"
    #Compute median # of filtered calls per trio
    if [ $( zcat "~{prefix}.~{condition_id}.trio_variants.txt.gz" | wc -l ) -gt 0 ]; then
      /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/helper_median_counts_per_trio.R \
        --ID "~{condition_id}" \
        "~{prefix}.~{condition_id}.trio_variants.txt.gz" \
        "~{trios_list}" \
        "~{prefix}.~{condition_id}.perTrio_distrib_stats.txt"
      med=$( fgrep -v "#" "~{prefix}.~{condition_id}.perTrio_distrib_stats.txt" | cut -f2 )
    else
      echo -e "#condition\thetsPerProband_median\thetsPerProband_Q1\thetsPerProband_Q2\n~{condition_id}\t0\t0\t0" \
      > "~{prefix}.~{condition_id}.perTrio_distrib_stats.txt"
      med=0
    fi
    #Run ROC if enough variants per proband
    echo -e "FINISHED FILTERING. FOUND $med MEDIAN QUALIFYING VARIANTS PER CHILD."
    if [ "$med" -gt ~{min_sv_per_proband_per_condition} ]; then
      echo -e "STARTING ROC OPTIMIZATION."
      /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/optimize_minGQ_ROC_v2.R \
        --prefix "~{condition_id}" \
        --fdr "~{roc_max_fdr}" \
        --minGQ "~{roc_min_gq}" \
        --maxGQ "~{roc_max_gq}" \
        --step "~{roc_step_gq}" \
        "~{prefix}.~{condition_id}.trio_variants.txt.gz" \
        "~{trios_list}" \
        "./"
      gzip "~{condition_id}.minGQ_ROC.data.txt"
    else
      echo -e "TOO FEW VARIANTS FOR ROC OPTIMIZATION."
      touch "~{condition_id}.minGQ_ROC.data.txt.gz"
      touch "~{condition_id}.minGQ_ROC.optimal.txt"
      touch "~{condition_id}.minGQ_ROC.plot.pdf"
    fi
  >>>

  output {
    File distrib_stats = "~{prefix}.~{condition_id}.perTrio_distrib_stats.txt"
    File roc_data = "~{condition_id}.minGQ_ROC.data.txt.gz"
    File roc_optimal = "~{condition_id}.minGQ_ROC.optimal.txt"
    File roc_plot = "~{condition_id}.minGQ_ROC.plot.pdf"
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


# Merge all ROC optimal cutoffs into single file for tree reconstruction
task CombineRocOptResults {
  input{
    Array[File] condition_optimizations
    Array[File] condition_distrib_stats
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 8,
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    find . -name "*.minGQ_ROC.optimal.txt" \
    | xargs -I {} cat {} | fgrep -v "#" | sort -Vk1,1 \
    > "~{prefix}.minGQ_condition_opts.txt" ||true
    find / -name "*.perTrio_distrib_stats.txt" \
    | xargs -I {} cat {} | fgrep -v "#" | sort -Vk1,1 \
    > "~{prefix}.minGQ_condition_distrib_stats.txt" ||true
  >>>

  output {
    File combined_optimizations = "~{prefix}.minGQ_condition_opts.txt"
    File combined_distrib_stats = "~{prefix}.minGQ_condition_distrib_stats.txt"
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
