version 1.0

import "ScatterCpxGenotyping.wdl" as GenotypeComplexContig
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "Utils.wdl" as util

workflow GenotypeComplexVariants {
  input {
    String cohort_name
    Array[String] batches
    File ped_file
    Array[File] depth_vcfs

    Boolean merge_vcfs = false
    Int? records_per_shard

    Array[File] complex_resolve_vcfs
    Array[File] complex_resolve_vcf_indexes

    Array[File] bincov_files

    Array[File] depth_gt_rd_sep_files
    Array[File] median_coverage_files

    File bin_exclude
    File contig_list
    File ref_dict

    Boolean use_hail = false
    String? gcs_project

    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_updates_docker
    String sv_pipeline_docker
    String sv_pipeline_hail_docker
    String sv_pipeline_rdtest_docker

    # overrides for mini tasks
    RuntimeAttr? runtime_override_concat

    # overrides for GenotypeComplexContig
    RuntimeAttr? runtime_override_ids_from_median
    RuntimeAttr? runtime_override_split_vcf_to_genotype
    RuntimeAttr? runtime_override_concat_cpx_cnv_vcfs
    RuntimeAttr? runtime_override_get_cpx_cnv_intervals
    RuntimeAttr? runtime_override_parse_genotypes
    RuntimeAttr? runtime_override_merge_melted_gts
    RuntimeAttr? runtime_override_split_bed_by_size
    RuntimeAttr? runtime_override_rd_genotype
    RuntimeAttr? runtime_override_concat_melted_genotypes
    RuntimeAttr? runtime_attr_ids_from_vcf
    RuntimeAttr? runtime_attr_subset_ped
    RuntimeAttr? runtime_override_preconcat
    RuntimeAttr? runtime_override_hail_merge
    RuntimeAttr? runtime_override_fix_header
  }

  scatter (i in range(length(batches))) {
    call util.GetSampleIdsFromVcf {
      input:
        vcf = depth_vcfs[i],
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_ids_from_vcf
    }
    call util.SubsetPedFile {
      input:
        ped_file = ped_file,
        sample_list = GetSampleIdsFromVcf.out_file,
        subset_name = batches[i],
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_subset_ped
    }
  }

  #Scatter per chromosome
  Array[String] contigs = transpose(read_tsv(contig_list))[0]
  scatter ( i in range(length(contigs)) ) {
    String contig = contigs[i]

    #Depth-based genotyping of complex intervals
    call GenotypeComplexContig.ScatterCpxGenotyping {
      input:
        bin_exclude=bin_exclude,
        vcf=complex_resolve_vcfs[i],
        records_per_shard=select_first([records_per_shard, 50000]),
        batches=batches,
        coverage_files=bincov_files,
        rd_depth_sep_cutoff_files=depth_gt_rd_sep_files,
        ped_file=ped_file,
        median_coverage_files=median_coverage_files,
        n_per_split_small=2500,
        n_per_split_large=250,
        n_rd_test_bins=100000,
        prefix="~{cohort_name}.~{contig}",
        contig=contig,
        ped_files=SubsetPedFile.ped_subset_file,
        ref_dict=ref_dict,
        use_hail=use_hail,
        gcs_project=gcs_project,
        linux_docker=linux_docker,
        sv_pipeline_updates_docker=sv_pipeline_updates_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_hail_docker=sv_pipeline_hail_docker,
        sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
        runtime_override_ids_from_median=runtime_override_ids_from_median,
        runtime_override_split_vcf_to_genotype=runtime_override_split_vcf_to_genotype,
        runtime_override_concat_cpx_cnv_vcfs=runtime_override_concat_cpx_cnv_vcfs,
        runtime_override_get_cpx_cnv_intervals=runtime_override_get_cpx_cnv_intervals,
        runtime_override_parse_genotypes=runtime_override_parse_genotypes,
        runtime_override_merge_melted_gts=runtime_override_merge_melted_gts,
        runtime_override_split_bed_by_size=runtime_override_split_bed_by_size,
        runtime_override_rd_genotype=runtime_override_rd_genotype,
        runtime_override_concat_melted_genotypes=runtime_override_concat_melted_genotypes,
        runtime_override_preconcat=runtime_override_preconcat,
        runtime_override_hail_merge=runtime_override_hail_merge,
        runtime_override_fix_header=runtime_override_fix_header
    }
  }

  #Merge resolved vcfs for QC
  if (merge_vcfs) {
    call MiniTasks.ConcatVcfs {
      input:
        vcfs=ScatterCpxGenotyping.cpx_depth_gt_resolved_vcf,
        vcfs_idx=ScatterCpxGenotyping.cpx_depth_gt_resolved_vcf_idx,
        allow_overlaps=true,
        outfile_prefix="~{cohort_name}.complex_genotype",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_concat
    }
  }

  #Final outputs
  output {
    Array[File] complex_genotype_vcfs = ScatterCpxGenotyping.cpx_depth_gt_resolved_vcf
    Array[File] complex_genotype_vcf_indexes = ScatterCpxGenotyping.cpx_depth_gt_resolved_vcf_idx
    File? complex_genotype_merged_vcf = ConcatVcfs.concat_vcf
    File? complex_genotype_merged_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}
