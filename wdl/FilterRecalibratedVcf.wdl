version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow FilterRecalibratedVcf {
  input {
    # VCF produced by RecalibrateGq, note that GT filtering here will be applied on top of
    # any already performed by RecalibrateGq.
    File vcf

    File ploidy_table
    String output_prefix

    # For testing
    File? filter_script

    Int records_per_shard

    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_override_scatter
    RuntimeAttr? runtime_attr_override_filter_vcf
    RuntimeAttr? runtime_override_concat_shards
  }

  call tasks_cohort.ScatterVcf {
    input:
      vcf=vcf,
      records_per_shard=records_per_shard,
      prefix="~{output_prefix}.scatter",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_override_scatter
  }

  scatter ( i in range(length(ScatterVcf.shards)) ) {
    call FilterVcf {
      input:
        vcf=ScatterVcf.shards[i],
        ploidy_table=ploidy_table,
        output_prefix="~{output_prefix}.filter.shard_~{i}",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_override_filter_vcf
    }
  }

  call tasks_cohort.ConcatVcfs {
    input:
      vcfs=FilterVcf.out,
      naive=true,
      outfile_prefix="~{output_prefix}.filtered",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_shards
  }

  output {
    File filtered_vcf = ConcatVcfs.concat_vcf
    File filtered_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}

task FilterVcf {
  input {
    File vcf
    File ploidy_table
    File? script
    String? args
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(100 + size(vcf, "GB") * 2),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{output_prefix}.vcf.gz"
    File out_index = "~{output_prefix}.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

    # Convert format
    python ~{default="/opt/sv-pipeline/scripts/apply_sl_filter.py" script} \
      --vcf ~{vcf} \
      --out ~{output_prefix}.vcf.gz \
      --ploidy-table ~{ploidy_table} \
      ~{args}

    tabix ~{output_prefix}.vcf.gz
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
