version 1.0

import "FormatVcfForGatk.wdl" as format
import "TasksMakeCohortVcf.wdl" as taskscohort
import "TasksClusterBatch.wdl" as taskscluster
import "TestUtils.wdl" as tu

workflow GenerateBatchMetrics {
  input {
    String batch

    File? depth_vcf
    File? dragen_vcf
    File? melt_vcf
    File? scramble_vcf
    File? wham_vcf
    File? manta_vcf

    File pe_file
    File sr_file
    File baf_file
    File rd_file

    File median_file
    File ped_file

    Int records_per_shard_agg = 10000

    String? additional_gatk_args_agg_pesr
    String? additional_gatk_args_agg_depth

    File? svtk_to_gatk_script

    File primary_contigs_list
    String chr_x
    String chr_y

    Float? java_mem_fraction

    File rmsk
    File segdups
    File reference_dict

    # Module metrics parameters
    # Run module metrics workflow at the end - on by default
    Boolean? run_module_metrics

    File? outlier_sample_ids # sample IDs to exclude from training

    String gatk_docker
    String sv_pipeline_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_create_ploidy
    RuntimeAttr? runtime_attr_aggregate_tests
    RuntimeAttr? runtime_attr_scatter_vcf
    RuntimeAttr? runtime_attr_format
    RuntimeAttr? runtime_attr_concat_vcfs
    RuntimeAttr? runtime_attr_agg_pesr
    RuntimeAttr? runtime_attr_agg_depth
    RuntimeAttr? runtime_attr_metrics_file_metrics
    RuntimeAttr? runtime_attr_annotate_overlap
  }

  Array[String] algorithms = ["depth", "dragen", "melt", "scramble", "wham", "manta"]
  Array[File?] vcfs = [depth_vcf, dragen_vcf, melt_vcf, scramble_vcf, wham_vcf, manta_vcf]
  String prefix = "~{batch}.batch_metrics"

  Array[File] vcfs_ = select_all([depth_vcf, manta_vcf, melt_vcf, scramble_vcf, wham_vcf])
  scatter (i in range(length(vcfs_))) {
    File vcfs_index_ = vcfs_[i] + ".tbi"
  }

  call taskscluster.CreatePloidyTableFromPed {
    input:
      ped_file=ped_file,
      contig_list=primary_contigs_list,
      retain_female_chr_y=false,
      chr_x=chr_x,
      chr_y=chr_y,
      output_prefix="~{batch}.ploidy",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_create_ploidy
  }

  call taskscohort.ConcatVcfs as ConcatInputVcfs {
    input:
      vcfs=vcfs_,
      vcfs_idx=vcfs_index_,
      allow_overlaps=true,
      outfile_prefix="~{prefix}.concat_input_vcfs",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_concat_vcfs
  }

  call taskscohort.ScatterVcf {
    input:
      vcf=ConcatInputVcfs.concat_vcf,
      records_per_shard = records_per_shard_agg,
      prefix = "~{prefix}.scatter_vcf",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_scatter_vcf
  }

  scatter ( i in range(length(ScatterVcf.shards)) ) {
    call format.FormatVcf {
      input:
        vcf=ScatterVcf.shards[i],
        ploidy_table=CreatePloidyTableFromPed.out,
        output_prefix="~{prefix}.format.shard_~{i}",
        script=svtk_to_gatk_script,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_format
    }
    call SVRegionOverlap {
      input:
        vcf = FormatVcf.out,
        vcf_index = FormatVcf.out_index,
        reference_dict = reference_dict,
        output_prefix = "~{prefix}.region_overlap.shard_~{i}",
        region_files = [segdups, rmsk],
        region_file_indexes = [segdups + ".tbi", rmsk + ".tbi"],
        region_names = ["SEGDUP", "RMSK"],
        java_mem_fraction=java_mem_fraction,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_annotate_overlap
    }
    call AggregateSVEvidence {
      input:
        vcf = SVRegionOverlap.out,
        vcf_index = SVRegionOverlap.out_index,
        output_prefix = "~{prefix}.aggregate_pesr.shard_~{i}",
        median_file = median_file,
        ploidy_table=CreatePloidyTableFromPed.out,
        pe_file = pe_file,
        pe_file_index = pe_file + ".tbi",
        sr_file = sr_file,
        sr_file_index = sr_file + ".tbi",
        baf_file = baf_file,
        baf_file_index = baf_file + ".tbi",
        chr_x = chr_x,
        chr_y = chr_y,
        additional_args=additional_gatk_args_agg_pesr,
        java_mem_fraction = java_mem_fraction,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_agg_pesr
    }
    call AggregateDepthEvidence {
      input:
        vcf = AggregateSVEvidence.out,
        vcf_index = AggregateSVEvidence.out_index,
        output_prefix = "~{prefix}.aggregate_depth.shard_~{i}",
        median_file = median_file,
        rd_file = rd_file,
        rd_file_index = rd_file + ".tbi",
        additional_args=additional_gatk_args_agg_depth,
        java_mem_fraction = java_mem_fraction,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_agg_depth
    }
  }

  call taskscohort.ConcatVcfs as ConcatOutputVcfs {
    input:
      vcfs=AggregateDepthEvidence.out,
      vcfs_idx=AggregateDepthEvidence.out_index,
      naive=true,
      outfile_prefix="~{prefix}.concat",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_concat_vcfs
  }

  call AggregateTests {
    input:
      vcf = ConcatOutputVcfs.concat_vcf,
      vcf_index = ConcatOutputVcfs.concat_vcf_idx,
      prefix = "~{prefix}.aggregate_tests",
      outlier_sample_ids = outlier_sample_ids,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_aggregate_tests
  }

  Boolean run_module_metrics_ = if defined(run_module_metrics) then select_first([run_module_metrics]) else true
  if (run_module_metrics_) {
    call tu.MetricsFileMetrics {
      input:
        metrics_file = AggregateTests.out,
        contig_list = select_first([primary_contigs_list]),
        common = false,
        prefix = "GenerateBatchMetrics.~{batch}",
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_metrics_file_metrics
    }
  }

  output {
    File metrics = AggregateTests.out
    File? metrics_file_batchmetrics = MetricsFileMetrics.out
    File ploidy_table = CreatePloidyTableFromPed.out
  }
}

task AggregateSVEvidence {
  input {
    File vcf
    File vcf_index
    String output_prefix

    File median_file
    File ploidy_table
    File? pe_file
    File? pe_file_index
    File? sr_file
    File? sr_file_index
    File? baf_file
    File? baf_file_index

    String chr_x
    String chr_y

    String? additional_args

    Float? java_mem_fraction
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }


  parameter_meta {
    pe_file: {
               localization_optional: true
             }
    sr_file: {
               localization_optional: true
             }
    baf_file: {
                localization_optional: true
              }
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 7.5,
                               disk_gb: ceil(10 + size(vcf, "GB") * 2.5 + size(sr_file, "GB")),
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

    gatk --java-options "-Xmx${JVM_MAX_MEM}" AggregateSVEvidence \
      -V ~{vcf} \
      -O ~{output_prefix}.vcf.gz \
      --median-coverage ~{median_file} \
      --ploidy-table ~{ploidy_table} \
      --x-chromosome-name ~{chr_x} \
      --y-chromosome-name ~{chr_y} \
      ~{"--discordant-pairs-file " + pe_file} \
      ~{"--split-reads-file " + sr_file} \
      ~{"--baf-file " + baf_file} \
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

task AggregateDepthEvidence {
  input {
    File vcf
    File vcf_index
    String output_prefix
    File median_file
    File rd_file
    File rd_file_index

    String? additional_args

    Float? java_mem_fraction
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 7.5,
                               disk_gb: ceil(10 + size(vcf, "GB") * 2.5 + size(rd_file, "GB")),
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

    gatk --java-options "-Xmx${JVM_MAX_MEM}" AggregateDepthEvidence \
      -V ~{vcf} \
      -O ~{output_prefix}.vcf.gz \
      --median-coverage ~{median_file} \
      ~{"--rd-file " + rd_file} \
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

task SVRegionOverlap {
  input {
    File vcf
    File vcf_index
    File reference_dict
    String output_prefix
    Array[File] region_files
    Array[File] region_file_indexes
    Array[String] region_names

    String? region_set_rule
    String? region_merging_rule
    Int? region_padding

    Boolean? suppress_overlap_fraction
    Boolean? suppress_endpoint_counts

    Float? java_mem_fraction

    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + size(vcf, "GB") * 2.0),
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

    gatk --java-options "-Xmx${JVM_MAX_MEM}" SVRegionOverlap \
      -V ~{vcf} \
      -O ~{output_prefix}.vcf.gz \
      --sequence-dictionary ~{reference_dict} \
      --track-intervals ~{sep=" --track-intervals " region_files} \
      --track-name ~{sep=" --track-name " region_names} \
      ~{"--region-set-rule " + region_set_rule} \
      ~{"--region-merging-rule " + region_merging_rule} \
      ~{"--region-padding " + region_padding} \
      --suppress-overlap-fraction ~{default="false" suppress_overlap_fraction} \
      --suppress-endpoint-counts ~{default="false" suppress_endpoint_counts}

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

task AggregateTests {
  input {
    File vcf
    File vcf_index
    String prefix
    File? outlier_sample_ids
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 7.5,
                               disk_gb: ceil(50 + size(vcf, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{prefix}.metrics.tsv"
  }
  command <<<
    /opt/sv-pipeline/02_evidence_assessment/02e_metric_aggregation/scripts/aggregate.py \
      -v ~{vcf} \
      ~{if defined(outlier_sample_ids) then "-o ~{outlier_sample_ids}" else ""} \
      ~{prefix}.metrics.tsv
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