version 1.0

import "Structs.wdl"

# Task to collect SV evidence on a single sample
task CollectSVEvidence {
  input {
    File bam_or_cram_file
    String sample_id

    # next three args are required for cram, but not for bam
    File? reference_fasta
    File? reference_index
    File? reference_dict

    File? primary_contigs_list # restrict output to this list of contigs

    File preprocessed_intervals # intervals for depth evidence collection
    File sd_locs_vcf # loci for site depth evidence collection

    Int site_depth_min_mapq = 6
    Int site_depth_min_baseq = 10
    File? gatk_jar_override
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
      bam_or_cram_file: {
        localization_optional: true
      }
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 30,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int command_mem_mb = ceil(mem_gb * 1000 - 500)

  output {
    File sr_file = "${sample_id}.sr.txt.gz"
    File sr_file_idx = "${sample_id}.sr.txt.gz.tbi"
    File pe_file = "${sample_id}.pe.txt.gz"
    File pe_file_idx = "${sample_id}.pe.txt.gz.tbi"
    File rd_file = "${sample_id}.rd.txt.gz"
    File rd_file_idx = "${sample_id}.rd.txt.gz.tbi"
    File sd_file = "${sample_id}.sd.txt.gz"
    File sd_file_idx = "${sample_id}.sd.txt.gz.tbi"
    File rd_metrics = "${sample_id}.rd_metrics.tsv"
  }
  command <<<

    set -euo pipefail

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_jar_override}

    /gatk/gatk --java-options "-Xmx~{command_mem_mb}m" CollectSVEvidence \
        -I ~{bam_or_cram_file} \
        --sample-name ~{sample_id} \
        -F ~{sd_locs_vcf} \
        -DI ~{preprocessed_intervals} \
        -SR "~{sample_id}.sr.txt.gz" \
        -PE "~{sample_id}.pe.txt.gz" \
        -SD "~{sample_id}.sd.txt.gz" \
        -RD "~{sample_id}.rd.txt.gz" \
        -DS "~{sample_id}.rd_metrics.tsv" \
        --site-depth-min-mapq "~{site_depth_min_mapq}" \
        --site-depth-min-baseq "~{site_depth_min_baseq}" \
        ~{"-R " + reference_fasta} \
        ~{"-L " + primary_contigs_list}

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

