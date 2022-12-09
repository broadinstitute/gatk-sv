version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow FormatVcfForGatk {
  input {
    File vcf
    File ploidy_table
    String prefix

    Int records_per_shard

    Boolean? run_svutils
    Boolean? run_formatter
    String? formatter_args

    # For testing
    File? svtk_to_gatk_script

    String sv_pipeline_docker
    String sv_utils_docker

    RuntimeAttr? runtime_attr_scatter
    RuntimeAttr? runtime_attr_svutils
    RuntimeAttr? runtime_attr_format
    RuntimeAttr? runtime_attr_concat
    RuntimeAttr? runtime_attr_concat_filter
  }

  Boolean run_svutils_ = select_first([run_svutils, true])
  Boolean run_formatter_ = select_first([run_formatter, true])

  call tasks_cohort.ScatterVcf {
    input:
      vcf=vcf,
      records_per_shard=records_per_shard,
      prefix="~{prefix}.scatter",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_scatter
  }

  scatter ( i in range(length(ScatterVcf.shards)) ) {
    if (run_svutils_) {
      call SvutilsFixVcf {
        input:
          vcf=ScatterVcf.shards[i],
          output_prefix="~{prefix}.svutils_~{i}",
          sv_utils_docker=sv_utils_docker,
          runtime_attr_override=runtime_attr_svutils
      }
    }
    if (run_formatter_) {
      call PreprocessVcf {
        input:
          vcf=select_first([SvutilsFixVcf.out, vcf]),
          ploidy_table=ploidy_table,
          args=formatter_args,
          output_prefix="~{prefix}.format_~{i}",
          script=svtk_to_gatk_script,
          sv_pipeline_docker=sv_pipeline_docker,
          runtime_attr_override=runtime_attr_format
      }
    }
  }

  call tasks_cohort.ConcatVcfs as ConcatFormattedVcfs {
    input:
      vcfs=if run_formatter_ then select_all(PreprocessVcf.out) else select_all(SvutilsFixVcf.out),
      naive=true,
      outfile_prefix="~{prefix}.gatk_formatted",
      sv_base_mini_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_concat
  }

  if (run_formatter_) {
    call tasks_cohort.ConcatVcfs as ConcatFilteredRecords {
      input:
        vcfs=select_all(PreprocessVcf.out),
        naive=true,
        outfile_prefix="~{prefix}.filtered_records",
        sv_base_mini_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_concat_filter
    }
  }

  output {
    File formatted_vcf = ConcatFormattedVcfs.concat_vcf
    File formatted_vcf_index = ConcatFormattedVcfs.concat_vcf_idx
    File? filtered_records_vcf = ConcatFilteredRecords.concat_vcf
    File? filtered_records_index = ConcatFilteredRecords.concat_vcf_idx
  }
}

task SvutilsFixVcf {
  input {
    File vcf
    String output_prefix
    String sv_utils_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + size(vcf, "GB") * 2),
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
    sv-utils fix-vcf ~{vcf} ~{output_prefix}.vcf.gz
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_utils_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task PreprocessVcf {
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
                               disk_gb: ceil(10 + size(vcf, "GB") * 2),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{output_prefix}.vcf.gz"
    File out_index = "~{output_prefix}.vcf.gz.tbi"
    File filtered = "~{output_prefix}.filtered_records.vcf.gz"
    File filtered_index = "~{output_prefix}.filtered_records.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

    # Convert format
    python ~{default="/opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py" script} \
      --vcf ~{vcf} \
      --out tmp.vcf.gz \
      --filter-out ~{output_prefix}.filtered_records.vcf.gz \
      --ploidy-table ~{ploidy_table} \
      ~{args}

    # TODO Filter invalid records with SVLEN=0, only needed for legacy runs that used svtk cluster in ClusterBatch
    bcftools view --no-version -i 'INFO/SVLEN="." || INFO/SVLEN>0' tmp.vcf.gz -Oz -o ~{output_prefix}.vcf.gz

    tabix ~{output_prefix}.vcf.gz
    tabix ~{output_prefix}.filtered_records.vcf.gz
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
