version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

import "GenotypeCpxCnvsPerBatch.wdl" as RunDepthGenotypePerBatch
import "TasksMakeCohortVcf.wdl" as MiniTasks

# Workflow to perform depth-based genotyping for a single vcf shard scattered
# across batches on predicted CPX CNVs
workflow GenotypeCpxCnvs {
  input {
    File bin_exclude
    File vcf
    Array[String] batches
    Array[File] coverage_files
    Array[File] rd_depth_sep_cutoff_files
    Array[File] ped_files
    Array[File] median_coverage_files
    Int n_per_split_small
    Int n_per_split_large
    Int n_rd_test_bins
    String prefix
    File ped_file
    String contig
    File ref_dict

    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_get_cpx_cnv_intervals
    RuntimeAttr? runtime_override_parse_genotypes

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_merge_melted_gts

    # overrides for RunDepthGenotypePerBatch
    RuntimeAttr? runtime_override_ids_from_median
    RuntimeAttr? runtime_override_split_bed_by_size
    RuntimeAttr? runtime_override_rd_genotype
    RuntimeAttr? runtime_override_concat_melted_genotypes
  }

  String contig_prefix = prefix + "." + contig

  # Convert VCF to bed of CPX CNV intervals
  call GetCpxCnvIntervals {
    input:
      vcf=vcf,
      prefix=contig_prefix,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_get_cpx_cnv_intervals
  }

  # Scatter over each batch (row) in gt_input_files and run depth genotyping
  scatter (i in range(length(batches))) {
    call RunDepthGenotypePerBatch.GenotypeCpxCnvsPerBatch as GenotypeBatch {
      input:
        bin_exclude=bin_exclude,
        cpx_bed=GetCpxCnvIntervals.cpx_cnv_bed,
        batch=batches[i],
        coverage_file=coverage_files[i],
        rd_depth_sep_cutoff=rd_depth_sep_cutoff_files[i],
        ped_file=ped_files[i],
        median_file=median_coverage_files[i],
        n_per_split_small=n_per_split_small,
        n_per_split_large=n_per_split_large,
        n_rd_test_bins=n_rd_test_bins,
        ref_dict=ref_dict,
        linux_docker=linux_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
        runtime_override_ids_from_median=runtime_override_ids_from_median,
        runtime_override_split_bed_by_size=runtime_override_split_bed_by_size,
        runtime_override_rd_genotype=runtime_override_rd_genotype,
        runtime_override_concat_melted_genotypes=runtime_override_concat_melted_genotypes
    }
  }

  # Merge melted genotypes across all batches
  call MiniTasks.ZcatCompressedFiles as MergeMeltedGts {
    input:
      shards=GenotypeBatch.genotypes,
      outfile_name=contig_prefix + ".CPX_intervals.merged_rd_genos.bed.gz",
      filter_command="sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_merge_melted_gts
  }

  # Parse genotyping results
  call ParseGenotypes {
    input:
      vcf=vcf,
      intervals=GetCpxCnvIntervals.cpx_cnv_bed,
      genotypes=MergeMeltedGts.outfile,
      prefix=contig_prefix,
      ped_file=ped_file,
      contig=contig,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_parse_genotypes
  }

  # Final output
  output {
    File cpx_depth_gt_resolved_vcf = ParseGenotypes.cpx_depth_gt_resolved_vcf
    File cpx_depth_gt_resolved_vcf_idx = ParseGenotypes.cpx_depth_gt_resolved_vcf_idx
    File reclassification_table = ParseGenotypes.reclassification_table
    File interval_genotype_counts_table = ParseGenotypes.gt_counts_table
  }
}


# Get CNV intervals from complex SV for depth genotyping
task GetCpxCnvIntervals {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_file = prefix + ".complex_CNV_intervals_to_test.bed.gz"

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(vcf, "GiB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + 2.0 * compression_factor)),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail
    
    /opt/sv-pipeline/04_variant_resolution/scripts/gather_cpx_intervals_for_rd_gt.sh \
      ~{vcf} \
      ~{output_file}
  >>>

  output {
    File cpx_cnv_bed = output_file
  }
}


# Parse genotyping results
task ParseGenotypes {
  input {
    File vcf
    File intervals
    File genotypes
    File ped_file
    String prefix
    String contig
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size([vcf, intervals, genotypes, ped_file], "GiB")
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0 + 5.0 * input_size,
    disk_gb: ceil(5 + input_size * 60),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail
    
    /opt/sv-pipeline/04_variant_resolution/scripts/process_posthoc_cpx_depth_regenotyping.sh \
      -R ~{prefix}.CPXregenotyping_reclassification_table.~{contig}.txt \
      -G ~{prefix}.CPXregenotyping_raw_genotype_counts_table.~{contig}.txt \
      ~{vcf} \
      ~{intervals} \
      ~{genotypes} \
      ~{ped_file} \
      ~{prefix}.postCPXregenotyping.~{contig}.vcf.gz
    tabix ~{prefix}.postCPXregenotyping.~{contig}.vcf.gz
  >>>

  output {
    File cpx_depth_gt_resolved_vcf = "~{prefix}.postCPXregenotyping.~{contig}.vcf.gz"
    File cpx_depth_gt_resolved_vcf_idx = "~{prefix}.postCPXregenotyping.~{contig}.vcf.gz.tbi"
    File reclassification_table = "~{prefix}.CPXregenotyping_reclassification_table.~{contig}.txt"
    File gt_counts_table = "~{prefix}.CPXregenotyping_raw_genotype_counts_table.~{contig}.txt"
  }
}
