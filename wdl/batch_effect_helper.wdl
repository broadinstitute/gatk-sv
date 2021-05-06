version 1.0

import "Structs.wdl"

workflow check_batch_effects {
  input{
    File freq_table
    String batch1
    String batch2
    String prefix
    Int variants_per_shard
    String sv_pipeline_docker
  }
  # Shard frequency table
  call ShardTable {
    input:
      freq_table=freq_table,
      variants_per_shard=variants_per_shard,
      sv_pipeline_docker=sv_pipeline_docker
  }

  # Scatter over shards and compute AF correlations for each variant
  scatter ( shard in ShardTable.shards ) {
    call CompareBatches {
      input:
        freq_table=shard,
        batch1=batch1,
        batch2=batch2,
        prefix=prefix,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  # Combine shards, perform bonferroni correction to determine significant batch effects, and plot AF correlation scatter
  call CombineShards {
    input:
      freq_tables=CompareBatches.results,
      batch1=batch1,
      batch2=batch2,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker
  }

  # Outputs
  output {
    File comparison_table = CombineShards.merged_table
    File batch_effect_variants = CombineShards.batch_effect_variants
    File scatterplots_tarball = CombineShards.correlation_scatterplots_tarball
  }
}


# Shard a frequency table into an even number of evenly sized shards
task ShardTable {
  input{
    File freq_table
    Int variants_per_shard
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    #Split variant lines
    zcat ~{freq_table} | sed '1d' | \
    split -l ~{variants_per_shard} --numeric-suffixes=00001 -a 5 /dev/stdin freq_table_shard_ || true
    #Add header & gzip each shard
    zcat ~{freq_table} | sed -n '1p' > header.txt
    maxshard=$( find / -name "freq_table_shard_*" | awk -v FS="_" '{ print $NF }' \
                | sort -Vrk1,1 | sed -n '1p' || true )
    for i in $( seq -w 00001 "$maxshard" ); do
      cat header.txt "freq_table_shard_$i" \
      | gzip -c \
      > "freq_table_shard_$i.txt.gz" || true
    done
  >>>

  output {
    Array[File] shards = glob("freq_table_shard_*.txt.gz")
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


# Compare AF stats per variant between a pair of batches
task CompareBatches {
  input{
    File freq_table
    String batch1
    String batch2
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/find_batch_effects.shard_helper.R \
      ~{freq_table} \
      "~{batch1}" \
      "~{batch2}" \
      "~{prefix}.~{batch1}_vs_~{batch2}.results.txt"
    gzip "~{prefix}.~{batch1}_vs_~{batch2}.results.txt"
  >>>

  output {
    File results = "~{prefix}.~{batch1}_vs_~{batch2}.results.txt.gz"
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


# Merge sharded comparison results and perform analysis for batch effects
task CombineShards {
  input{
    Array[File] freq_tables
    String batch1
    String batch2
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    #Write header
    zcat ~{freq_tables[0]} | sed -n '1p' > header.txt || true
    #Iterate over files and cat
    while read file; do
      zcat "$file" | sed '1d'
    done < ~{write_lines(freq_tables)} \
    | cat header.txt - \
    | gzip -c \
    > "~{prefix}.~{batch1}_vs_~{batch2}.AF_comparison_table.txt.gz" || true
    #Analyze
    mkdir "~{batch1}_vs_~{batch2}"
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/find_batch_effects.R \
      "~{prefix}.~{batch1}_vs_~{batch2}.AF_comparison_table.txt.gz" \
      "~{batch1}" \
      "~{batch2}" \
      "~{batch1}_vs_~{batch2}/~{prefix}" ||true
    gzip -f "~{batch1}_vs_~{batch2}/~{prefix}.~{batch1}_vs_~{batch2}.freq_table_wBonferroni.txt"
    tar -czvf "~{batch1}_vs_~{batch2}.tar.gz" \
      "~{batch1}_vs_~{batch2}"
  >>>

  output {
    File merged_table = "~{batch1}_vs_~{batch2}/~{prefix}.~{batch1}_vs_~{batch2}.freq_table_wBonferroni.txt.gz"
    File batch_effect_variants = "~{batch1}_vs_~{batch2}/~{prefix}.~{batch1}_vs_~{batch2}.batch_effect_variants.txt"
    File correlation_scatterplots_tarball = "~{batch1}_vs_~{batch2}.tar.gz"
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

