version 1.0

import "ScatterCpxGenotyping.wdl" as GenotypeComplexContig
import "Tasks0506.wdl" as MiniTasks

workflow Module0506ComplexGenotype {
  input {
    String cohort_name
    Array[String] batches
    Array[File] ped_files
    File merged_ped_file

    Boolean merge_vcfs = false

    Array[File] complex_resolve_vcfs
    Array[File] complex_resolve_vcf_indexes

    Array[File] bincov_files

    Array[File] depth_gt_rd_sep_files
    Array[File] median_coverage_files

    File bin_exclude
    File contig_list
    Int max_shards_per_chrom
    File ref_dict

    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker

    # overrides for mini tasks
    RuntimeAttr? runtime_override_ids_from_vcf
    RuntimeAttr? runtime_override_merge_fam_file_list
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
        n_master_vcf_shards=200,
        n_master_min_vars_per_vcf_shard=5000,
        batches=batches,
        coverage_files=bincov_files,
        rd_depth_sep_cutoff_files=depth_gt_rd_sep_files,
        merged_ped_file=merged_ped_file,
        median_coverage_files=median_coverage_files,
        n_per_split_small=2500,
        n_per_split_large=250,
        n_rd_test_bins=100000,
        prefix=cohort_name,
        contig=contig,
        ped_files=ped_files,
        ref_dict=ref_dict,
        linux_docker=linux_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
        runtime_override_ids_from_median=runtime_override_ids_from_median,
        runtime_override_split_vcf_to_genotype=runtime_override_split_vcf_to_genotype,
        runtime_override_ids_from_vcf=runtime_override_ids_from_vcf,
        runtime_override_concat_cpx_cnv_vcfs=runtime_override_concat_cpx_cnv_vcfs,
        runtime_override_get_cpx_cnv_intervals=runtime_override_get_cpx_cnv_intervals,
        runtime_override_parse_genotypes=runtime_override_parse_genotypes,
        runtime_override_merge_melted_gts=runtime_override_merge_melted_gts,
        runtime_override_split_bed_by_size=runtime_override_split_bed_by_size,
        runtime_override_rd_genotype=runtime_override_rd_genotype,
        runtime_override_concat_melted_genotypes=runtime_override_concat_melted_genotypes
    }
  }

  #Merge resolved vcfs for QC
  if (merge_vcfs) {
    call MiniTasks.ConcatVcfs {
      input:
        vcfs=ScatterCpxGenotyping.cpx_depth_gt_resolved_vcf,
        vcfs_idx=ScatterCpxGenotyping.cpx_depth_gt_resolved_vcf_idx,
        merge_sort=true,
        outfile_prefix="~{cohort_name}.0506_complex",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_concat
    }
  }

  #Final outputs
  output {
    Array[File] complex_genotype_vcfs = ScatterCpxGenotyping.cpx_depth_gt_resolved_vcf
    Array[File] complex_genotype_vcf_indexes = ScatterCpxGenotyping.cpx_depth_gt_resolved_vcf_idx
    File? merged_vcf = ConcatVcfs.concat_vcf
    File? merged_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}
