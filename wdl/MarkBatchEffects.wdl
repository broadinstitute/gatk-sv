# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Identify and mark variants that exhibit substantial differences in non-ref GT
# rates between different sets of user-specified samples (e.g., batches, cohorts)


version 1.0

import "Utilities.wdl" as Utilities

workflow MarkBatchEffects {
  input {
    File vcf
    File vcf_idx
    File group_membership_tsv             # Two-column .tsv mapping group names (e.g., baches, cohorts) to sample IDs. One line per sample.

    # Parameters for mark_batch_effects.py
    Int min_samples_per_group = 10
    Float lower_nonref_gt_freq = 0.1
    Float upper_nonref_gt_freq = 0.5
    Boolean flag_failing_records = false  # Should records with observed batch effects have a non-PASS FILTER assigned?
    String? custom_filter_id
    String? custom_filter_description
    String? custom_info_flag_id

    # Other options
    File? ref_fai                         # Used to determine contigs for parallelization. By default will take all contigs in vcf header
    String? out_vcf_prefix

    # Resource options
    Float? mark_batch_effects_mem_gb
    Int? mark_batch_effects_cpu_cores
    Int? mark_batch_effects_disk_gb
    Float? merge_mem_gb
    Int? merge_cpu_cores
    Int? merge_disk_gb

    String g2c_pipeline_docker
  }

  # Get contigs for scatter depending on user input
  if (defined(ref_fai)) {
    call Utilities.GetContigsFromFai {
      input:
        ref_fai = select_first([ref_fai]),
        docker = g2c_pipeline_docker
    }
  }
  if (!defined(ref_fai)) {
    call Utilities.GetContigsFromVcfHeader {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        docker = g2c_pipeline_docker
    }
  }
  Array[String] contigs = select_first([GetContigsFromFai.contigs, 
                                        GetContigsFromVcfHeader.contigs])

  # Apply script to each chromosome in parallel
  scatter (contig in contigs) {
    call MarkBatchEffectsSingleContig {
      input:
        vcf = vcf,
        vcf_idx = select_first([vcf_idx, GetContigsFromVcfHeader.vcf_idx_out]),
        group_membership_tsv = group_membership_tsv,
        contig = contig,
        min_n = min_samples_per_group,
        lower_freq = lower_nonref_gt_freq,
        upper_freq = upper_nonref_gt_freq,
        strict = flag_failing_records,
        custom_filter_id = custom_filter_id,
        custom_filter_description = custom_filter_description,
        custom_info_flag_id = custom_info_flag_id,
        mem_gb = mark_batch_effects_mem_gb,
        cpu_cores = mark_batch_effects_cpu_cores,
        disk_gb = mark_batch_effects_disk_gb,
        docker = g2c_pipeline_docker
    }
  }

  # Merge outputs from per-chromosome scatter
  call Utilities.ConcatVcfs {
    input:
      vcfs = MarkBatchEffectsSingleContig.vcf_out,
      vcf_idxs = MarkBatchEffectsSingleContig.vcf_out_idx,
      out_prefix = select_first([out_vcf_prefix, basename(vcf, ".vcf.gz") + ".bfx_marked"]),
      bcftools_concat_options = "--naive",
      mem_gb = merge_mem_gb,
      cpu_cores = merge_cpu_cores,
      disk_gb = merge_disk_gb,
      bcftools_docker = g2c_pipeline_docker
  }

  output {
    File output_vcf = ConcatVcfs.merged_vcf
    File output_vcf_idx = ConcatVcfs.merged_vcf_idx
  }
}


task MarkBatchEffectsSingleContig {
  input {
    File vcf
    File? vcf_idx
    File group_membership_tsv
    String contig

    Int min_n
    Float lower_freq
    Float upper_freq
    Boolean strict
    String? custom_filter_id
    String? custom_filter_description
    String? custom_info_flag_id

    Float mem_gb = 3.5
    Int cpu_cores = 2
    Int? disk_gb

    String docker
  }

  String vcf_basename = basename(vcf, ".vcf.gz")
  String in_vcf_name = vcf_basename + "." + contig + ".in.vcf.gz"
  String out_vcf_name = vcf_basename + "." + contig + ".bfx_marked.vcf.gz"

  Int default_disk_gb = ceil(2.5 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    if [ ~{defined(vcf_idx)} == "false" ]; then
      tabix -p vcf -f ~{vcf}
    fi

    # Subset VCF to contig of interest
    tabix -h ~{vcf} "~{contig}" | bgzip -c > "~{in_vcf_name}"

    # Mark batch effects and update AC/AN/AF
    /opt/pancan_germline_wgs/scripts/variant_filtering/mark_batch_effects.py \
      --input-vcf "~{in_vcf_name}" \
      --group-membership ~{group_membership_tsv} \
      ~{if (strict) then "--strict" else ""} \
      ~{if defined(custom_filter_id) then "--filter-id \"~{custom_filter_id}\"" else ""} \
      ~{if defined(custom_filter_description) then "--filter-description \"~{custom_filter_description}\"" else ""} \
      ~{if defined(custom_info_flag_id) then "--info-flag-id \"~{custom_info_flag_id}\"" else ""} \
      --min-samples ~{min_n} \
      --lower-freq ~{lower_freq} \
      --upper-freq ~{upper_freq} \
    | bcftools +fill-tags - -Oz -o "~{out_vcf_name}" -- -t AC,AN,AF
    tabix -p vcf -f "~{out_vcf_name}"
  >>>

  output {
    File vcf_out = "~{out_vcf_name}"
    File vcf_out_idx = "~{out_vcf_name}.tbi"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: cpu_cores
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    preemptible: 3
  }
}
