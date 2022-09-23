##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

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
    String sv_mini_docker
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
    tabix -p vcf ~{vcf};
    mkdir splits
    tabix ~{vcf} ~{chrom} | split -a ~{suffix_len} -d -l ~{split_size} - splits/~{batch}.~{algorithm}.split.

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
    docker: sv_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

task MergeAllosomes {
  input {
    File male_test
    File female_test
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
    while read split; do
      sed -e '1d' $split;
    done < ~{write_lines(stats)} | cat <(head -n1 ~{stats[0]}) - > ~{prefix}.stats
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
