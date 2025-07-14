version 1.0

import "Structs.wdl"
import "FormatVcfForGatk.wdl" as FormatVcf
import "TasksClusterBatch.wdl" as ClusterTasks
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow ScatterFormatCluster {
  input {
    File vcf1
    File vcf1_index
    File vcf2
    File vcf2_index
    Array[String] contigs
    String prefix

    # FormatVcfForGatk inputs
    File ped_file
    File contig_list
    Int records_per_shard = 40000
    File? contigs_header
    String? formatter_args
    File? bothside_pass_list
    File? background_fail_list
    String? chr_x
    String? chr_y
    File? svtk_to_gatk_script

    # SVCluster inputs
    File ploidy_table
    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    # Optional SVCluster overlap parameters
    Boolean? fast_mode
    Boolean? omit_members
    Boolean? enable_cnv
    Boolean? default_no_call
    String? algorithm
    String? insertion_length_summary_strategy
    String? breakpoint_summary_strategy
    String? alt_allele_summary_strategy
    Float? defrag_padding_fraction
    Float? defrag_sample_overlap
    Float? depth_sample_overlap
    Float? depth_interval_overlap
    Float? depth_size_similarity
    Int? depth_breakend_window
    Float? mixed_sample_overlap
    Float? mixed_interval_overlap
    Float? mixed_size_similarity
    Int? mixed_breakend_window
    Float? pesr_sample_overlap
    Float? pesr_interval_overlap
    Float? pesr_size_similarity
    Int? pesr_breakend_window
    Float? java_mem_fraction
    String? additional_args
    String? variant_prefix

    # Docker images
    String sv_base_mini_docker
    String sv_pipeline_docker
    String gatk_docker

    # Runtime overrides
    RuntimeAttr? runtime_attr_extract_contig
    RuntimeAttr? runtime_attr_create_ploidy
    RuntimeAttr? runtime_attr_scatter
    RuntimeAttr? runtime_attr_format
    RuntimeAttr? runtime_attr_concat_format
    RuntimeAttr? runtime_attr_cluster
    RuntimeAttr? runtime_attr_concat_final
  }

  scatter (contig in contigs) {
    # Extract contig from each VCF
    call ExtractContig as ExtractContig1 {
      input:
        vcf = vcf1,
        vcf_index = vcf1_index,
        contig = contig,
        prefix = "~{prefix}.vcf1.~{contig}",
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_extract_contig
    }

    call ExtractContig as ExtractContig2 {
      input:
        vcf = vcf2,
        vcf_index = vcf2_index,
        contig = contig,
        prefix = "~{prefix}.vcf2.~{contig}",
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_extract_contig
    }

    # Format each contig VCF for GATK
    call FormatVcf.FormatVcfForGatk as FormatVcf1 {
      input:
        vcf = ExtractContig1.out,
        prefix = "~{prefix}.vcf1.~{contig}.formatted",
        ped_file = ped_file,
        records_per_shard = records_per_shard,
        contig_list = contig_list,
        contigs_header = contigs_header,
        formatter_args = formatter_args,
        bothside_pass_list = bothside_pass_list,
        background_fail_list = background_fail_list,
        chr_x = chr_x,
        chr_y = chr_y,
        svtk_to_gatk_script = svtk_to_gatk_script,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_create_ploidy = runtime_attr_create_ploidy,
        runtime_attr_scatter = runtime_attr_scatter,
        runtime_attr_format = runtime_attr_format,
        runtime_override_concat = runtime_attr_concat_format
    }

    call FormatVcf.FormatVcfForGatk as FormatVcf2 {
      input:
        vcf = ExtractContig2.out,
        prefix = "~{prefix}.vcf2.~{contig}.formatted",
        ped_file = ped_file,
        records_per_shard = records_per_shard,
        contig_list = contig_list,
        contigs_header = contigs_header,
        formatter_args = formatter_args,
        bothside_pass_list = bothside_pass_list,
        background_fail_list = background_fail_list,
        chr_x = chr_x,
        chr_y = chr_y,
        svtk_to_gatk_script = svtk_to_gatk_script,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_create_ploidy = runtime_attr_create_ploidy,
        runtime_attr_scatter = runtime_attr_scatter,
        runtime_attr_format = runtime_attr_format,
        runtime_override_concat = runtime_attr_concat_format
    }

    # Cluster the pair of formatted VCFs
    call ClusterTasks.SVCluster {
      input:
        vcfs = [FormatVcf1.gatk_formatted_vcf, FormatVcf2.gatk_formatted_vcf],
        ploidy_table = ploidy_table,
        output_prefix = "~{prefix}.~{contig}.clustered",
        contig = contig,
        fast_mode = fast_mode,
        omit_members = omit_members,
        enable_cnv = enable_cnv,
        default_no_call = default_no_call,
        algorithm = algorithm,
        insertion_length_summary_strategy = insertion_length_summary_strategy,
        breakpoint_summary_strategy = breakpoint_summary_strategy,
        alt_allele_summary_strategy = alt_allele_summary_strategy,
        defrag_padding_fraction = defrag_padding_fraction,
        defrag_sample_overlap = defrag_sample_overlap,
        depth_sample_overlap = depth_sample_overlap,
        depth_interval_overlap = depth_interval_overlap,
        depth_size_similarity = depth_size_similarity,
        depth_breakend_window = depth_breakend_window,
        mixed_sample_overlap = mixed_sample_overlap,
        mixed_interval_overlap = mixed_interval_overlap,
        mixed_size_similarity = mixed_size_similarity,
        mixed_breakend_window = mixed_breakend_window,
        pesr_sample_overlap = pesr_sample_overlap,
        pesr_interval_overlap = pesr_interval_overlap,
        pesr_size_similarity = pesr_size_similarity,
        pesr_breakend_window = pesr_breakend_window,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        reference_dict = reference_dict,
        java_mem_fraction = java_mem_fraction,
        additional_args = additional_args,
        variant_prefix = variant_prefix,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_cluster
    }
  }

  # Concatenate all clustered results
  call MiniTasks.ConcatVcfs {
    input:
      vcfs = SVCluster.out,
      vcfs_idx = SVCluster.out_index,
      allow_overlaps = false,
      naive = false,
      outfile_prefix = "~{prefix}.final_clustered",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_final
  }

  output {
    File final_clustered_vcf = ConcatVcfs.concat_vcf
    File final_clustered_vcf_index = ConcatVcfs.concat_vcf_idx
    Array[File] per_contig_clustered_vcfs = SVCluster.out
    Array[File] per_contig_clustered_vcf_indices = SVCluster.out_index
  }
}

task ExtractContig {
  input {
    File vcf
    File vcf_index
    String contig
    String prefix
    String sv_base_mini_docker
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
    File out = "~{prefix}.vcf.gz"
    File out_index = "~{prefix}.vcf.gz.tbi"
  }

  command <<<
    set -euo pipefail
    bcftools view --regions ~{contig} ~{vcf} -Oz -o ~{prefix}.vcf.gz
    tabix ~{prefix}.vcf.gz
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
} 