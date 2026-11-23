version 1.0

import "Structs.wdl"
import "FilterGenotypes.wdl" as filter_genotypes
import "MainVcfQc.wdl" as qc

# Multi-contig version of FilterGenotypes. Runs genotype filtering independently on each per-contig vcf,
# then (optionally) runs MainVcfQc once across all filtered per-contig vcfs.
workflow FilterGenotypesByContig {

  input {
    # Per-contig vcfs, ordered to match contig_list
    Array[File] vcfs
    File contig_list
    String output_prefix

    File ploidy_table
    Float no_call_rate_cutoff = 0.05 # Set to 1 to disable NCR filtering
    File sl_cutoff_table # Sl filtering cutoff table
    File? optimized_sl_cutoff_table # Overrides sl_cutoff_table if provided
    String? sl_filter_args # Explicitly set SL arguments - see apply_sl_filter.py

    Int filter_vcf_records_per_shard = 20000
    String? header_drop_fields # Uses FilterGenotypes default if not provided

    # For MainVcfQc
    File primary_contigs_fai
    File? ped_file
    Array[Array[String]]? site_level_comparison_datasets    # Array of two-element arrays, one per dataset, each of format [prefix, gs:// path to directory with one BED per population]
    Array[Array[String]]? sample_level_comparison_datasets  # Array of two-element arrays, one per dataset, each of format [prefix, gs:// path to per-sample tarballs]
    File? sample_renaming_tsv # File with mapping to rename per-sample benchmark sample IDs for compatibility with cohort
    Boolean run_qc = true
    String qc_bcftools_preprocessing_options = "-i 'FILTER=\"PASS\" || FILTER=\"MULTIALLELIC\"'"
    Int qc_sv_per_shard = 2500
    Int qc_samples_per_shard = 600
    RuntimeAttr? runtime_override_plot_qc_per_family
    RuntimeAttr? runtime_override_per_sample_benchmark_plot

    String sv_base_mini_docker
    String sv_pipeline_docker
  }

  Array[String] contigs = read_lines(contig_list)

  scatter (i in range(length(contigs))) {
    call filter_genotypes.FilterGenotypes {
      input:
        vcf=vcfs[i],
        output_prefix="~{output_prefix}.~{contigs[i]}",
        ploidy_table=ploidy_table,
        no_call_rate_cutoff=no_call_rate_cutoff,
        sl_cutoff_table=sl_cutoff_table,
        optimized_sl_cutoff_table=optimized_sl_cutoff_table,
        sl_filter_args=sl_filter_args,
        filter_vcf_records_per_shard=filter_vcf_records_per_shard,
        header_drop_fields=header_drop_fields,
        primary_contigs_fai=primary_contigs_fai,
        run_qc=false,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  if (run_qc) {
    call qc.MainVcfQc {
      input:
        vcfs=FilterGenotypes.filtered_vcf,
        prefix="~{output_prefix}.filter_genotypes",
        ped_file=ped_file,
        bcftools_preprocessing_options=qc_bcftools_preprocessing_options,
        primary_contigs_fai=primary_contigs_fai,
        site_level_comparison_datasets=site_level_comparison_datasets,
        sample_level_comparison_datasets=sample_level_comparison_datasets,
        sample_renaming_tsv=sample_renaming_tsv,
        sv_per_shard=qc_sv_per_shard,
        runtime_override_per_sample_benchmark_plot=runtime_override_per_sample_benchmark_plot,
        runtime_override_plot_qc_per_family=runtime_override_plot_qc_per_family,
        samples_per_shard=qc_samples_per_shard,
        sv_pipeline_qc_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  output {
    Array[File] filtered_vcfs = FilterGenotypes.filtered_vcf
    Array[File] filtered_vcf_indexes = FilterGenotypes.filtered_vcf_index
    File? main_vcf_qc_tarball = MainVcfQc.sv_vcf_qc_output
  }
}
