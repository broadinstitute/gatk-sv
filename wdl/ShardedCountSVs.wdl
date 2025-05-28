version 1.0

import "PlotSVCountsPerSample.wdl" as plot
import "CollectQcVcfWide.wdl" as collect
import "TasksMakeCohortVcf.wdl" as tasks

workflow ShardedCountSVs {
  input {
    File vcf
    String prefix

    String? bcftools_preprocessing_options
    Int records_per_shard = 20000

    String sv_pipeline_docker
    String sv_base_mini_docker

  }

  call tasks.ScatterVcf {
    input:
      vcf = vcf,
      records_per_shard = records_per_shard,
      prefix = "~{prefix}.scatter",
      sv_pipeline_docker = sv_pipeline_docker
  }

  scatter (i in range(length(ScatterVcf.shards))) {
    if (defined(bcftools_preprocessing_options)) {
      call collect.PreprocessVcf {
        input:
          vcf = ScatterVcf.shards[i],
          bcftools_preprocessing_options=select_first([bcftools_preprocessing_options]),
          prefix = "~{prefix}.shard_~{i}.preprocessed",
          sv_base_mini_docker = sv_base_mini_docker
      }
    }

    call plot.CountSVsPerSamplePerType {
      input:
        vcf = select_first([PreprocessVcf.outvcf, ScatterVcf.shards[i]]),
        prefix = "~{prefix}.shard_~{i}",
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  call SumSVCounts {
    input:
      counts = CountSVsPerSamplePerType.sv_counts,
      prefix = prefix,
      sv_pipeline_docker = sv_pipeline_docker
  }

  output {
    File sv_counts = SumSVCounts.sv_counts_summed
  }
}

task SumSVCounts {
  input {
    Array[File] counts
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10 + ceil(size(counts, "GiB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    python <<CODE
import pandas as pd
counts = dict()
for file_path in ["~{sep='", "' counts}"]:
  with open(file_path, 'r') as inp:
    inp.readline()  # skip header
    for line in inp:
      sample, svtype, count = line.strip().split('\t')
      if sample in counts:
        if svtype in counts[sample]:
          counts[sample][svtype] += int(count)
        else:
          counts[sample][svtype] = int(count)
      else:
        counts[sample] = dict()
        counts[sample][svtype] = int(count)
with open("~{prefix}.sv_counts_summed.tsv", 'w') as out:
  out.write("sample\tsvtype\tcount\n")
  for sample in counts:
    for svtype in counts[sample]:
      out.write(f"{sample}\t{svtype}\t{counts[sample][svtype]}\n")
CODE

  >>>

  output {
    File sv_counts_summed = "~{prefix}.sv_counts_summed.tsv"
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
