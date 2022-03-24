version 1.0

import "Structs.wdl"

task SplitVCF {
  input {
    File vcf
    String batch
    String algorithm
    String chrom
    Int split_size
    Int suffix_len
    String sv_base_mini_docker
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
    Array[File] split_vcfs = glob("headered/*.vcf.gz")
  }
  command <<<

    set -euo pipefail
    tabix -p vcf ~{vcf}
    mkdir splits
    tabix ~{vcf} ~{chrom} | split -a ~{suffix_len} -d -l ~{split_size} - splits/~{batch}.~{algorithm}.split.

    if [ ! -f splits/~{batch}.~{algorithm}.split.`printf '%0~{suffix_len}d' 0` ]; then
      touch splits/~{batch}.~{algorithm}.split.`printf '%0~{suffix_len}d' 0`
    fi

    #Add headers to each split
    mkdir headered
    for split in splits/~{batch}.~{algorithm}.split.*; do
      [ -e "$split" ] || continue
      splitname=$(basename $split)
      cat <(zcat ~{vcf} | sed -n -e '/^#/p') ${split} | bgzip -c > headered/${splitname}.vcf.gz
    done
  
  >>>
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

task GetCommonVCF {
  input{
      File vcf
      Int cnv_size_cutoff
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
    File common_vcf = "common_SVs.vcf.gz"
    File common_vcf_tbi = "common_SVs.vcf.gz.tbi"
  }

  command <<<

    set -euxo pipefail

    { zcat ~{vcf} || true; } | grep -m1 "^#CHROM" | sed -e 's/\t/\n/g' | tail -n+10 > sample_name

    declare int sample_counts=$(wc -l sample_name | awk '{print $1}')
    sample_count_cutoff=`expr $sample_counts / 2 - 1`

    svtk vcf2bed ~{vcf} info.bed
    awk '{print $6}' info.bed | awk -F"," '{print NF-1}' > sample_counts
    paste info.bed sample_counts > info.V2.bed
    awk '{if ($NF > '$sample_count_cutoff' && $5!="DEL" && $5!="DUP") print}' info.V2.bed | cut -f4 > common_SVID
    awk '{if ($NF > '$sample_count_cutoff' && $5=="DEL" && $3-$2<~{cnv_size_cutoff}) print}' info.V2.bed | cut -f4 >> common_SVID
    awk '{if ($NF > '$sample_count_cutoff' && $5=="DUP" && $3-$2<~{cnv_size_cutoff}) print}' info.V2.bed | cut -f4 >> common_SVID

    /usr/local/bin/bcftools filter -i 'ID=@common_SVID' ~{vcf} > common_SVs.vcf

    vcf-sort common_SVs.vcf | bgzip > common_SVs.vcf.gz
    tabix common_SVs.vcf.gz
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

task MergeAllosomes {
  input {
    File male_test
    File female_test
    File male_only_ids_list
    String chrom
    String male_only_expr
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
    File merged_test = "${basename(male_test)}.merged.csv"
  }
  command <<<

    set -euo pipefail
    python3 <<CODE
    import pandas as pd
    males = pd.read_table("~{male_test}")
    females = pd.read_table("~{female_test}")
    male_only_ids = set()
    with open("~{male_only_ids_list}", 'r') as male_only_file:
      for line in male_only_file:
        male_only_ids.add(line.strip())
    if "~{chrom}" == 'Y' or "~{chrom}" == 'chrY':
      males.to_csv("~{basename(male_test)}.merged.csv", sep='\t', index=False, na_rep='NA')
    else:
      male_only = ~{male_only_expr}
      females.loc[male_only] = males
      females.to_csv("~{basename(male_test)}.merged.csv", sep='\t', index=False, na_rep='NA')
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

# Combine per-chromosome stat results into single table
task MergeStats {
  input {
    Array[File] stats
    String prefix
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

  command <<<
    set -euo pipefail
    touch ~{prefix}.stats
    for f in ~{sep=" "  stats}; do sed -n 1p $f >~{prefix}.stats; break; done
    for f in ~{sep=" "  stats}; do sed 1d $f >>~{prefix}.stats; done
  >>>

  output {
    File merged_stats = "~{prefix}.stats"
  }
  
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
