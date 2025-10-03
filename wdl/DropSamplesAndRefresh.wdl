version 1.0

import "SubsetVcfBySamples.wdl" as subset
import "SanitizeHeader.wdl" as sanitize
import "AnnotateVcf.wdl" as anno
import "GetVcfStats.wdl" as stats
import "MainVcfQc.wdl" as qc
import "ApplyNCRAndRefArtifactFilters.wdl" as ncr
import "TasksMakeCohortVcf.wdl" as tasks

workflow DropSamplesAndRefresh {
  input {
    Array[File] vcfs
    File keep_samples
    File related_samples
    String prefix

    File primary_contigs_list
    File ploidy_table

    # NCR
    File apply_filters_script

    # SanitizeHeader
    String drop_fields
    File sample_id_rename_map

    # AnnotateVcf
    File protein_coding_gtf
    Int anno_sv_per_shard
    File ped_file
    File par_bed
    File sample_pop_assignments
    File ref_bed
    File ref_prefix
    Array[String] ref_populations

    # MainVcfQc
    File primary_contigs_fai
    Int qc_sv_per_shard


    String sv_pipeline_docker
    String sv_base_mini_docker
    String gatk_docker
  }

  call subset.SubsetVcfBySamples as DropSamples {
    input:
      vcfs=vcfs,
      list_of_samples=keep_samples,
      remove_samples=false,
      remove_private_sites=true,
      sv_base_mini_docker=sv_base_mini_docker
  }

  call ncr.ApplyNCRAndRefArtifactFilters {
    input:
      vcfs=DropSamples.vcfs_subset,
      apply_filters_script=apply_filters_script,
      primary_contigs_list=primary_contigs_list,
      cohort_id=prefix,
      ploidy_table=ploidy_table,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_base_mini_docker=sv_base_mini_docker
  }

  call subset.SubsetVcfBySamples as SubsetToUnrelated {
    input:
      vcfs=ApplyNCRAndRefArtifactFilters.filtered_vcfs,
      list_of_samples=related_samples,
      remove_samples=true,
      remove_private_sites=true,
      sv_base_mini_docker=sv_base_mini_docker
  }

  # Supply unique inputs (genetic ancestries, external AF files, GTF, GATK docker, etc) through subworkflow inputs
  call anno.AnnotateVcf {
    input:
      vcfs=ApplyNCRAndRefArtifactFilters.filtered_vcfs,
      contig_list=primary_contigs_list,
      prefix=prefix,
      protein_coding_gtf=protein_coding_gtf,
      ped_file=ped_file,
      par_bed=par_bed,
      sample_pop_assignments=sample_pop_assignments,
      ref_bed=ref_bed,
      ref_prefix=ref_prefix,
      population=ref_populations,
      sv_per_shard=anno_sv_per_shard,
      gatk_docker=gatk_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call anno.AnnotateVcf as AnnotateUnrelated{
    input:
      vcfs=SubsetToUnrelated.vcfs_subset,
      contig_list=primary_contigs_list,
      prefix="~{prefix}.unrelated",
      protein_coding_gtf=protein_coding_gtf,
      ped_file=ped_file,
      par_bed=par_bed,
      sample_pop_assignments=sample_pop_assignments,
      ref_bed=ref_bed,
      ref_prefix=ref_prefix,
      population=ref_populations,
      sv_per_shard=anno_sv_per_shard,
      gatk_docker=gatk_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker
  }

  # Supply unique inputs (sample ID rename map, drop fields) through subworkflow inputs
  call sanitize.SanitizeHeader {
    input:
      vcfs=AnnotateVcf.annotated_vcfs,
      prefix=prefix,
      drop_fields=drop_fields,
      sample_id_rename_map=sample_id_rename_map,
      primary_contigs_list=primary_contigs_list,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call sanitize.SanitizeHeader as SanitizeUnrelated {
    input:
      vcfs=AnnotateUnrelated.annotated_vcfs,
      prefix="~{prefix}.unrelated",
      drop_fields=drop_fields,
      sample_id_rename_map=sample_id_rename_map,
      primary_contigs_list=primary_contigs_list,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call tasks.ConcatVcfs {
    input:
      vcfs=SanitizeHeader.vcf_header_sanitized,
      outfile_prefix="~{prefix}.sites_only",
      sites_only=true,
      sv_base_mini_docker=sv_base_mini_docker
  }

  call tasks.ConcatVcfs as ConcatUnrelated {
    input:
      vcfs=SanitizeUnrelated.vcf_header_sanitized,
      outfile_prefix="~{prefix}.unrelated.sites_only",
      sites_only=true,
      sv_base_mini_docker=sv_base_mini_docker
  }

  call stats.GetVcfStats {
    input:
      vcfs=SanitizeHeader.vcf_header_sanitized,
      prefix=prefix,
      contigs_list=primary_contigs_list,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker
  }

  # Supply unique inputs (ie primary contigs fai) directly to subworkflow
  call qc.MainVcfQc as SiteQc {
    input:
      vcfs=SanitizeHeader.vcf_header_sanitized,
      bcftools_preprocessing_options="-i 'FILTER=\"PASS\" || FILTER=\"MULTIALLELIC\"'",
      prefix=prefix,
      do_per_sample_qc=false,
      samples_per_shard=600,
      random_seed=9,
      sv_per_shard=qc_sv_per_shard,
      primary_contigs_fai=primary_contigs_fai,
      sv_per_shard=qc_sv_per_shard,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_qc_docker=sv_pipeline_docker
  }

  call qc.MainVcfQc as UnrelatedQc {
    input:
      vcfs=SanitizeUnrelated.vcf_header_sanitized,
      bcftools_preprocessing_options="-i 'FILTER=\"PASS\" || FILTER=\"MULTIALLELIC\"'",
      prefix="~{prefix}.unrelated",
      do_per_sample_qc=false,
      samples_per_shard=600,
      random_seed=9,
      sv_per_shard=qc_sv_per_shard,
      primary_contigs_fai=primary_contigs_fai,
      sv_per_shard=qc_sv_per_shard,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_qc_docker=sv_pipeline_docker
  }


  output {
    Array[File] refreshed_vcfs = SanitizeHeader.vcf_header_sanitized
    Array[File] refreshed_vcf_indexes = SanitizeHeader.vcf_header_sanitized_index

    Array[File] unrelated_vcfs = SanitizeUnrelated.vcf_header_sanitized
    Array[File] unrelated_vcf_indexes = SanitizeUnrelated.vcf_header_sanitized_index

    File sites_only_vcf = ConcatVcfs.concat_vcf
    File sites_only_vcf_index = ConcatVcfs.concat_vcf_idx

    File sites_only_unrelated_vcf = ConcatUnrelated.concat_vcf
    File sites_only_unrelated_vcf_index = ConcatUnrelated.concat_vcf_idx

    File? per_sample_sv_counts = GetVcfStats.sv_counts
    File sites_info = GetVcfStats.sites_info

    File sites_qc_tarball = SiteQc.sv_vcf_qc_output
    File bed_file = SiteQc.vcf2bed_output
    File duplicate_records = SiteQc.duplicate_records_output
    File duplicate_counts_output = SiteQc.duplicate_counts_output

    File unrelated_sites_qc_tarball = UnrelatedQc.sv_vcf_qc_output
    File unrelated_bed_file = UnrelatedQc.vcf2bed_output

    Array[File] vid_rename_maps = ApplyNCRAndRefArtifactFilters.id_rename_maps
  }
}
