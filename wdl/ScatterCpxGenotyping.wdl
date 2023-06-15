version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

import "GenotypeCpxCnvs.wdl" as GenotypeCpx
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "HailMerge.wdl" as HailMerge

# Workflow to perform depth-based genotyping for a single vcf shard scattered 
# across batches on predicted CPX CNVs
workflow ScatterCpxGenotyping {
  input {
    File bin_exclude
    File vcf
    Int records_per_shard
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

    Boolean use_hail
    String? gcs_project

    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_updates_docker
    String sv_pipeline_docker
    String sv_pipeline_hail_docker
    String sv_pipeline_rdtest_docker

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_split_vcf_to_genotype
    RuntimeAttr? runtime_override_concat_cpx_cnv_vcfs

    RuntimeAttr? runtime_override_preconcat
    RuntimeAttr? runtime_override_hail_merge
    RuntimeAttr? runtime_override_fix_header

    # overrides for GenotypeCpx
    RuntimeAttr? runtime_override_ids_from_median
    RuntimeAttr? runtime_override_get_cpx_cnv_intervals
    RuntimeAttr? runtime_override_parse_genotypes
    RuntimeAttr? runtime_override_merge_melted_gts
    RuntimeAttr? runtime_override_split_bed_by_size
    RuntimeAttr? runtime_override_rd_genotype
    RuntimeAttr? runtime_override_concat_melted_genotypes
  }

  String contig_prefix = prefix + "." + contig

  # Shard VCF into even slices
  call MiniTasks.ScatterVcf as SplitVcfToGenotype {
    input:
      vcf=vcf,
      prefix=contig_prefix,
      records_per_shard=records_per_shard,
      sv_pipeline_docker=sv_pipeline_updates_docker,
      runtime_attr_override=runtime_override_split_vcf_to_genotype
  }

  # Scatter genotyping over shards
  scatter ( shard in SplitVcfToGenotype.shards ) {
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
        ped_file=ped_file,
        contig=contig,
        ref_dict=ref_dict,
        linux_docker=linux_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
        runtime_override_ids_from_median=runtime_override_ids_from_median,
        runtime_override_get_cpx_cnv_intervals=runtime_override_get_cpx_cnv_intervals,
        runtime_override_parse_genotypes=runtime_override_parse_genotypes,
        runtime_override_merge_melted_gts=runtime_override_merge_melted_gts,
        runtime_override_split_bed_by_size=runtime_override_split_bed_by_size,
        runtime_override_rd_genotype=runtime_override_rd_genotype,
        runtime_override_concat_melted_genotypes=runtime_override_concat_melted_genotypes
    }
  }

  if (use_hail) {
    call HailMerge.HailMerge as ConcatCpxCnvVcfsHail {
      input:
        vcfs=GenotypeShard.cpx_depth_gt_resolved_vcf,
        prefix="~{prefix}.regenotyped",
        gcs_project=gcs_project,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_hail_docker=sv_pipeline_hail_docker,
        runtime_override_preconcat=runtime_override_preconcat,
        runtime_override_hail_merge=runtime_override_hail_merge,
        runtime_override_fix_header=runtime_override_fix_header
    }
  }
  if (!use_hail) {
    call MiniTasks.ConcatVcfs as ConcatCpxCnvVcfs {
      input:
        vcfs=GenotypeShard.cpx_depth_gt_resolved_vcf,
        vcfs_idx=GenotypeShard.cpx_depth_gt_resolved_vcf_idx,
        allow_overlaps=true,
        outfile_prefix="~{prefix}.regenotyped",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_concat_cpx_cnv_vcfs
    }
  }

  # Output merged VCF
  output {
    File cpx_depth_gt_resolved_vcf = select_first([ConcatCpxCnvVcfs.concat_vcf, ConcatCpxCnvVcfsHail.merged_vcf])
    File cpx_depth_gt_resolved_vcf_idx = select_first([ConcatCpxCnvVcfs.concat_vcf_idx, ConcatCpxCnvVcfsHail.merged_vcf_index])
  }
 }
