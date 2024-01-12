version 1.0

import "Structs.wdl"

# Post-hoc fix for CleanVcf runs that used the earliest version of ReshardVcf without concatenating
#  the high SR background and bothsides SR support lists prior to CleanVcf.

# IMPORTANT: This workflow is intended as a hotfix for this particular issue, and will only apply flags
#  to INS/CPX records. Be sure to review logs for info/warning messages.

workflow ApplySRFlagsHotfix {
  input {
    Array[File] cleaned_vcfs
    Array[File] complex_resolve_vcfs
    Array[File] complex_resolve_bothside_pass_lists
    Array[File] complex_resolve_background_fail_lists

    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_apply_sr_flags
  }


  scatter (i in range(length(cleaned_vcfs))) {
    String prefix = basename(cleaned_vcfs[i], ".vcf.gz") + ".apply_sr_flags"
    call ApplySRFlagsTask {
      input:
        prefix = prefix,
        cleaned_vcf = cleaned_vcfs[i],
        complex_resolve_vcf = complex_resolve_vcfs[i],
        complex_resolve_bothside_pass_lists = complex_resolve_bothside_pass_lists,
        complex_resolve_background_fail_lists = complex_resolve_background_fail_lists,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_apply_sr_flags
    }
  }

  output {
    Array[File] apply_sr_flags_vcfs = ApplySRFlagsTask.out
    Array[File] apply_sr_flags_vcf_indexes = ApplySRFlagsTask.out_index
  }
}


task ApplySRFlagsTask {
  input {
    String prefix
    File cleaned_vcf
    File complex_resolve_vcf
    Array[File] complex_resolve_bothside_pass_lists
    Array[File] complex_resolve_background_fail_lists
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               mem_gb: 3.75,
                               disk_gb: ceil(10.0 + (size(cleaned_vcf, "GiB") * 2) + size(complex_resolve_vcf, "GiB")),
                               cpu_cores: 1,
                               preemptible_tries: 3,
                               max_retries: 1,
                               boot_disk_gb: 10
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<

    set -euo pipefail

    while read SHARD; do
      if [ -n "$SHARD" ]; then
        cat "$SHARD"
      fi
    done < ~{write_lines(complex_resolve_bothside_pass_lists)} > bothside.txt

    while read SHARD; do
      if [ -n "$SHARD" ]; then
        cat "$SHARD"
      fi
    done < ~{write_lines(complex_resolve_background_fail_lists)} > background.txt

    python /opt/sv-pipeline/scripts/apply_sr_flags.py \
      --cleaned-vcf ~{cleaned_vcf} \
      --cpx-resolve-vcf ~{complex_resolve_vcf} \
      --out ~{prefix}.vcf.gz \
      --bothsides-list bothside.txt \
      --high-background-list background.txt

    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
    File out_index = "~{prefix}.vcf.gz.tbi"
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
