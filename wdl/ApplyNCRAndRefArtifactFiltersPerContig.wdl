version 1.0

import "TasksMakeCohortVcf.wdl" as tasks
import "Structs.wdl"

workflow ApplyNCRAndRefArtifactFiltersPerContig {
  input {
    File vcf
    String prefix

    File ploidy_table
    Int records_per_shard = 20000

    Float? no_call_rate_cutoff
    Boolean filter_reference_artifacts = true
    Boolean remove_zero_carrier_sites = true

    File? apply_filters_script

    String sv_pipeline_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_scatter_vcf
    RuntimeAttr? runtime_attr_apply_filters
    RuntimeAttr? runtime_attr_concat_vcfs
  }

  call tasks.ScatterVcf {
    input:
      vcf = vcf,
      records_per_shard = records_per_shard,
      prefix = "~{prefix}.scatter",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_scatter_vcf
  }

  scatter (i in range(length(ScatterVcf.shards))) {
    call ApplyFilters {
      input:
        vcf = ScatterVcf.shards[i],
        prefix = "~{prefix}.shard_~{i}",
        ploidy_table = ploidy_table,
        no_call_rate_cutoff = no_call_rate_cutoff,
        filter_reference_artifacts = filter_reference_artifacts,
        remove_zero_carrier_sites = remove_zero_carrier_sites,
        apply_filters_script = apply_filters_script,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_apply_filters
    }
  }

  call tasks.ConcatVcfs {
    input:
      vcfs = ApplyFilters.filtered_vcf,
      naive = true,
      outfile_prefix = prefix,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_vcfs
  }

  output {
    File filtered_vcf = ConcatVcfs.concat_vcf
    File filtered_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}


task ApplyFilters {
  input {
    File vcf
    String prefix
    File ploidy_table
    Float? no_call_rate_cutoff
    Boolean filter_reference_artifacts
    Boolean remove_zero_carrier_sites
    File? apply_filters_script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 5,
    disk_gb: ceil(10.0 + 2 * size(vcf, "GiB")),
    boot_disk_gb: 30,
    preemptible_tries: 1,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python ~{select_first([apply_filters_script, "/opt/sv-pipeline/scripts/apply_ncr_and_ref_artifact_filters.py"])} \
      --vcf ~{vcf} \
      --out ~{prefix}.vcf.gz \
      --ploidy-table ~{ploidy_table} \
      --ncr-threshold ~{no_call_rate_cutoff} \
      ~{if (filter_reference_artifacts) then "--filter-reference-artifacts" else ""} \
      ~{if (remove_zero_carrier_sites) then "--remove-zero-carrier-sites" else ""}

  >>>

  output {
    File filtered_vcf = "~{prefix}.vcf.gz"
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

