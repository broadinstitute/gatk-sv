version 1.0

import "Structs.wdl"

task SVCluster {
  input {
    Array[File] vcfs
    Array[File] vcf_indexes
    String output_name
    File ref_dict

    String? contig
    String? vid_prefix

    String algorithm = "SINGLE_LINKAGE"
    Boolean disable_dict_validation = true  # Disabling is very slow
    Boolean omit_members = false
    Boolean convert_inv = false
    Boolean enable_cnv = false

    Float? defrag_sample_overlap
    Float? defrag_padding_fraction
    String? breakpoint_summary_strategy

    Float? depth_overlap_fraction
    Float? mixed_overlap_fraction
    Float? pesr_overlap_fraction
    Int? depth_breakend_window
    Int? mixed_breakend_window
    Int? pesr_breakend_window
    Float? depth_sample_overlap
    Float? mixed_sample_overlap
    Float? pesr_sample_overlap

    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    vcfs: {
           localization_optional: true
         }
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

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{output_name}.vcf.gz"
    File out_index = "~{output_name}.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

    # Create arguments file
    touch args.txt
    while read line; do
      echo "-V $line" >> args.txt
    done < ~{write_lines(vcfs)}

    /gatk/gatk --java-options "-Xmx~{java_mem_mb}m" SVCluster \
      --arguments_file args.txt \
      -O ~{output_name}.vcf.gz \
      --algorithm ~{algorithm} \
      --sequence-dictionary ~{ref_dict} \
      ~{if disable_dict_validation then "--disable-sequence-dictionary-validation" else ""} \
      ~{if omit_members then "--omit-members" else ""} \
      ~{"--variant-prefix " + vid_prefix} \
      ~{"-L " + contig} \
      ~{"--breakpoint-summary-strategy " + breakpoint_summary_strategy} \
      ~{"--min-sample-set-fraction-overlap " + defrag_sample_overlap} \
      ~{"--defrag-padding-fraction " + defrag_padding_fraction} \
      ~{"--depth-overlap-fraction " + depth_overlap_fraction} \
      ~{"--mixed-overlap-fraction " + mixed_overlap_fraction} \
      ~{"--pesr-overlap-fraction " + pesr_overlap_fraction} \
      ~{"--depth-breakend-window " + depth_breakend_window} \
      ~{"--mixed-breakend-window " + mixed_breakend_window} \
      ~{"--pesr-breakend-window " + pesr_breakend_window} \
      ~{"--depth-sample-overlap " + depth_sample_overlap} \
      ~{"--mixed-sample-overlap " + mixed_sample_overlap} \
      ~{"--pesr-sample-overlap " + pesr_sample_overlap} \
      ~{if convert_inv then "--convert-inv-to-bnd" else ""} \
      ~{if enable_cnv then "--enable-cnv" else ""}
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}