version 1.0

import "Structs.wdl"
import "BatchEvidenceMerging.wdl" as bem
import "MakeBincovMatrix.wdl" as mbm

workflow GATKSVCombineSampleEvidence {
  input {
    # Batch info
    String batch
    Array[String] samples

    # Preprocessed sample SV VCFs
    Array[File] vcfs

    # PE/SR files
    Array[File] counts
    Array[File] pe_files
    Array[File] sr_files
    File inclusion_bed

    # Global files
    File ped_file
    File genome_file
    File ref_fasta_dict

    # Runtime parameters
    String sv_base_mini_docker
    String sv_base_docker
    String sv_pipeline_docker
    String gatk_docker

    RuntimeAttr? runtime_attr_merge_sv
    RuntimeAttr? runtime_attr_shard_pe
    RuntimeAttr? runtime_attr_merge_pe
    RuntimeAttr? runtime_attr_shard_sr
    RuntimeAttr? runtime_attr_merge_sr
    RuntimeAttr? runtime_attr_set_sample
    RuntimeAttr? runtime_attr_override_make_bincov
  }

  scatter (i in range(length(vcfs))) {
    File vcf_indexes_ = vcfs[i] + ".tbi"
  }

  call SVPreprocessRecords {
    input:
      vcfs = vcfs,
      vcf_indexes = vcf_indexes_,
      batch = batch,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_merge_sv
  }

  call mbm.MakeBincovMatrix {
    input:
      samples = samples,
      count_files = counts,
      batch = batch,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_base_docker = sv_base_docker,
      runtime_attr_override = runtime_attr_override_make_bincov
  }

  call bem.EvidenceMerging {
    input:
      samples = samples,
      PE_files = pe_files,
      SR_files = sr_files,
      inclusion_bed = inclusion_bed,
      genome_file = genome_file,
      batch = batch,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_set_sample = runtime_attr_set_sample,
      runtime_attr_shard_pe = runtime_attr_shard_pe,
      runtime_attr_merge_pe = runtime_attr_merge_pe,
      runtime_attr_shard_sr = runtime_attr_shard_sr,
      runtime_attr_merge_sr = runtime_attr_merge_sr
  }

  output {
    File merged_sv_vcf = SVPreprocessRecords.out
    File merged_sv_vcf_index = SVPreprocessRecords.out

    File merged_rd = MakeBincovMatrix.merged_bincov
    File merged_rd_index = MakeBincovMatrix.merged_bincov_idx
    File merged_sr = EvidenceMerging.merged_SR
    File merged_sr_index = EvidenceMerging.merged_SR_idx
    File merged_pe = EvidenceMerging.merged_PE
    File merged_pe_index = EvidenceMerging.merged_PE_idx
  }
}

task SVPreprocessRecords {
  input {
    Array[File] vcfs
    Array[File] vcf_indexes
    String batch
    String gatk_path = "/gatk/gatk"
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 7.5,
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{batch}.merged.vcf.gz"
    File out_index = "~{batch}.merged.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

    # Create arguments file
    touch args.txt
    while read line; do
      echo "-V $line" >> args.txt
    done < ~{write_lines(vcfs)}

    ~{gatk_path} --java-options "-Xmx~{java_mem_mb}m" SVPreprocessRecords \
      --arguments_file args.txt \
      --output ~{batch}.merged.vcf.gz \
      --disable-sequence-dictionary-validation
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