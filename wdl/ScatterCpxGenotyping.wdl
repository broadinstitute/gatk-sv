version 1.0

# based on snapshot 12
# https://portal.firecloud.org/#methods/Talkowski-SV/04b_scatter_CPX_genotyping/12/wdl

# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License

import "GenotypeCpxCnvs.wdl" as GenotypeCpx
import "Tasks0506.wdl" as MiniTasks

# Workflow to perform depth-based genotyping for a single vcf shard scattered 
# across batches on predicted CPX CNVs from 04b
workflow ScatterCpxGenotyping {
  input {
    File bin_exclude
    File vcf
    Int n_master_vcf_shards
    Int n_master_min_vars_per_vcf_shard
    Array[String] batches
    Array[File] coverage_files
    Array[File] rd_depth_sep_cutoff_files
    Array[File] ped_files
    Array[File] median_coverage_files
    Int n_per_split_small
    Int n_per_split_large
    Int n_rd_test_bins
    String prefix
    File merged_ped_file
    String contig
    File ref_dict

    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_ids_from_vcf
    RuntimeAttr? runtime_override_split_vcf_to_genotype
    RuntimeAttr? runtime_override_concat_cpx_cnv_vcfs

    # overrides for GenotypeCpx
    RuntimeAttr? runtime_override_get_cpx_cnv_intervals
    RuntimeAttr? runtime_override_parse_genotypes
    RuntimeAttr? runtime_override_merge_melted_gts
    RuntimeAttr? runtime_override_split_bed_by_size
    RuntimeAttr? runtime_override_rd_genotype
    RuntimeAttr? runtime_override_concat_melted_genotypes
  }

  String contig_prefix = prefix + "." + contig

  # Shard VCF into even slices
  call MiniTasks.SplitVcf as SplitVcfToGenotype {
    input:
      vcf=vcf,
      prefix=contig_prefix + ".shard_",
      n_shards=n_master_vcf_shards,
      min_vars_per_shard=n_master_min_vars_per_vcf_shard,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_split_vcf_to_genotype
  }

  # Scatter genotyping over shards
  scatter ( shard in SplitVcfToGenotype.vcf_shards ) {
    # Run genotyping
    call GenotypeCpx.GenotypeCpxCnvs as GenotypeShard {
      input:
        bin_exclude=bin_exclude,
        vcf=shard,
        batches=batches,
        coverage_files=coverage_files,
        rd_depth_sep_cutoff_files=rd_depth_sep_cutoff_files,
        ped_files=ped_files,
        median_coverage_files=median_coverage_files,
        n_per_split_large=n_per_split_large,
        n_per_split_small=n_per_split_small,
        n_rd_test_bins=n_rd_test_bins,
        prefix=prefix,
        merged_ped_file=merged_ped_file,
        contig=contig,
        ref_dict=ref_dict,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
        runtime_override_get_cpx_cnv_intervals=runtime_override_get_cpx_cnv_intervals,
        runtime_override_ids_from_vcf=runtime_override_ids_from_vcf,
        runtime_override_parse_genotypes=runtime_override_parse_genotypes,
        runtime_override_merge_melted_gts=runtime_override_merge_melted_gts,
        runtime_override_split_bed_by_size=runtime_override_split_bed_by_size,
        runtime_override_rd_genotype=runtime_override_rd_genotype,
        runtime_override_concat_melted_genotypes=runtime_override_concat_melted_genotypes
    }
  }

  # Merge VCF shards
  call MiniTasks.ConcatVcfs as ConcatCpxCnvVcfs {
    input:
      vcfs=GenotypeShard.cpx_depth_gt_resolved_vcf,
      vcfs_idx=GenotypeShard.cpx_depth_gt_resolved_vcf_idx,
      outfile_prefix=contig_prefix + ".resolved",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_cpx_cnv_vcfs
  }

  # Output merged VCF
  output {
    File cpx_depth_gt_resolved_vcf = ConcatCpxCnvVcfs.concat_vcf
    File cpx_depth_gt_resolved_vcf_idx = ConcatCpxCnvVcfs.concat_vcf_idx
  }
 }
