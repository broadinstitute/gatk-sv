version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow SVConcordance {
  input {
    File eval_vcf
    File truth_vcf

    File ploidy_table
    String cohort

    Boolean? run_svutils_truth_vcf
    Boolean? run_formatter_truth_vcf
    String? formatter_truth_args

    Boolean? run_svutils_eval_vcf
    Boolean? run_formatter_eval_vcf
    String? formatter_eval_args

    # For testing
    File? svtk_to_gatk_script

    File contig_list
    File reference_dict

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_utils_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_svutils_truth
    RuntimeAttr? runtime_attr_format_truth
    RuntimeAttr? runtime_attr_svutils_eval
    RuntimeAttr? runtime_attr_format_eval
    RuntimeAttr? runtime_attr_sv_concordance
    RuntimeAttr? runtime_attr_postprocess
    RuntimeAttr? runtime_override_concat_shards
  }

  Boolean run_svutils_truth_vcf_ = select_first([run_svutils_truth_vcf, true])
  Boolean run_formatter_truth_vcf_ = select_first([run_formatter_truth_vcf, true])

  Boolean run_svutils_eval_vcf_ = select_first([run_svutils_eval_vcf, true])
  Boolean run_formatter_eval_vcf_ = select_first([run_formatter_eval_vcf, true])

  if (run_svutils_truth_vcf_) {
    call SvutilsFixVcf as SvutilsTruth {
      input:
        vcf=truth_vcf,
        output_prefix="~{cohort}.svutils_truth",
        sv_utils_docker=sv_utils_docker,
        runtime_attr_override=runtime_attr_svutils_truth
    }
  }
  if (run_formatter_truth_vcf_) {
    call PreprocessVcf as FormatTruth {
      input:
        vcf=select_first([SvutilsTruth.out, truth_vcf]),
        ploidy_table=ploidy_table,
        args=formatter_truth_args,
        output_prefix="~{cohort}.format_truth",
        script=svtk_to_gatk_script,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_format_truth
    }
  }

  if (run_svutils_eval_vcf_) {
    call SvutilsFixVcf as SvutilsEval {
      input:
        vcf=eval_vcf,
        output_prefix="~{cohort}.svutils_eval",
        sv_utils_docker=sv_utils_docker,
        runtime_attr_override=runtime_attr_svutils_eval
    }
  }
  if (run_formatter_eval_vcf_) {
    call PreprocessVcf as FormatEval {
      input:
        vcf=select_first([SvutilsEval.out, eval_vcf]),
        ploidy_table=ploidy_table,
        args=formatter_eval_args,
        output_prefix="~{cohort}.format_eval",
        script=svtk_to_gatk_script,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_format_eval
    }
  }

  Array[String] contigs = transpose(read_tsv(contig_list))[0]
  scatter (contig in contigs) {
    call SVConcordanceTask {
      input:
        eval_vcf=select_first([FormatEval.out, SvutilsEval.out, eval_vcf]),
        truth_vcf=select_first([FormatTruth.out, SvutilsTruth.out, truth_vcf]),
        output_prefix="~{cohort}.concordance.~{contig}",
        contig=contig,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_sv_concordance
    }
  }

  call tasks_cohort.ConcatVcfs {
    input:
      vcfs=SVConcordanceTask.out,
      vcfs_idx=SVConcordanceTask.out_index,
      naive=true,
      outfile_prefix="~{cohort}.concordance",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_shards
  }

  output {
    File concordance_vcf = ConcatVcfs.concat_vcf
    File concordance_vcf_index = ConcatVcfs.concat_vcf_idx
    File? filtered_eval_records_vcf = FormatEval.filtered
    File? filtered_eval_records_index =FormatEval.filtered_index
    File? filtered_truth_records_vcf = FormatTruth.filtered
    File? filtered_truth_records_index = FormatTruth.filtered_index
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

task SVConcordanceTask {
  input {
    File truth_vcf
    File eval_vcf
    String output_prefix
    File reference_dict
    String? contig
    String? additional_args
    Float? java_mem_fraction
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    truth_vcf: {
                 localization_optional: true
               }
    eval_vcf:  {
                 localization_optional: true
               }
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + size(eval_vcf, "GB") * 2 + size(truth_vcf, "GB")),
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

    function getJavaMem() {
    # get JVM memory in MiB by getting total memory from /proc/meminfo
    # and multiplying by java_mem_fraction
      cat /proc/meminfo \
        | awk -v MEM_FIELD="$1" '{
          f[substr($1, 1, length($1)-1)] = $2
        } END {
          printf "%dM", f[MEM_FIELD] * ~{default="0.85" java_mem_fraction} / 1024
        }'
    }
    JVM_MAX_MEM=$(getJavaMem MemTotal)
    echo "JVM memory: $JVM_MAX_MEM"

    gatk --java-options "-Xmx${JVM_MAX_MEM}" SVConcordance \
      ~{"-L " + contig} \
      --sequence-dictionary ~{reference_dict} \
      --eval ~{eval_vcf} \
      --truth ~{truth_vcf} \
      -O ~{output_prefix}.vcf.gz \
      --force-biallelic-dups \
      ~{additional_args}
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
