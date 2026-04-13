version 1.0

import "Structs.wdl"

workflow SVShell {
  input {
    File gcnv_model_tars_list
    File ref_pesr_split_files_list
    File ref_pesr_disc_files_list
    File ref_pesr_sd_files_list
    Array[File] genome_tracks
    File pesr_exclude_intervals
    File ref_panel_vcf
    File ref_panel_bincov_matrix
    File bin_exclude
    File cytobands
    File depth_exclude_list
    File mei_bed
    File manta_region_bed
  }

  Array[File] gcnv_model_tars = read_lines(gcnv_model_tars_list)

  Array[File] ref_pesr_split_files = read_lines(ref_pesr_split_files_list)
  scatter (ref_pesr_split_file in ref_pesr_split_files) {
    File ref_pesr_split_file_index = ref_pesr_split_file + ".tbi"
  }

  Array[File] ref_pesr_disc_files = read_lines(ref_pesr_disc_files_list)
  scatter (ref_pesr_disc_file in ref_pesr_disc_files) {
    File ref_pesr_disc_file_index = ref_pesr_disc_file + ".tbi"
  }

  Array[File] ref_pesr_sd_files = read_lines(ref_pesr_sd_files_list)
  scatter (ref_pesr_sd_file in ref_pesr_sd_files) {
    File ref_pesr_sd_file_index = ref_pesr_sd_file + ".tbi"
  }

  scatter (genome_track in genome_tracks) {
    File genome_track_index = genome_track + ".tbi"
  }

  File pesr_exclude_intervals_index = pesr_exclude_intervals + ".tbi"
  File ref_panel_vcf_index = ref_panel_vcf + ".tbi"
  File ref_panel_bincov_matrix_index = ref_panel_bincov_matrix + ".tbi"
  File bin_exclude_index = bin_exclude + ".tbi"
  File cytobands_index = cytobands + ".tbi"
  File depth_exclude_list_index = depth_exclude_list + ".tbi"
  File mei_bed_index = mei_bed + ".tbi"
  File manta_region_bed_index = manta_region_bed + ".tbi"


  call RunSVShell {
    input:
      gcnv_model_tars = gcnv_model_tars,
      ref_pesr_split_files = ref_pesr_split_files,
      ref_pesr_split_files_indices = ref_pesr_split_file_index,
      ref_pesr_disc_files = ref_pesr_disc_files,
      ref_pesr_disc_files_indices = ref_pesr_disc_file_index,
      ref_pesr_sd_files = ref_pesr_sd_files,
      ref_pesr_sd_files_indices = ref_pesr_sd_file_index,
      genome_tracks = genome_tracks,
      genome_tracks_indices = genome_track_index,
      pesr_exclude_intervals = pesr_exclude_intervals,
      pesr_exclude_intervals_index = pesr_exclude_intervals_index,
      ref_panel_vcf = ref_panel_vcf,
      ref_panel_vcf_index = ref_panel_vcf_index,
      ref_panel_bincov_matrix = ref_panel_bincov_matrix,
      ref_panel_bincov_matrix_index = ref_panel_bincov_matrix_index,
      bin_exclude = bin_exclude,
      bin_exclude_index = bin_exclude_index,
      cytobands = cytobands,
      cytobands_index = cytobands_index,
      depth_exclude_list = depth_exclude_list,
      depth_exclude_list_index = depth_exclude_list_index,
      mei_bed = mei_bed,
      mei_bed_index = mei_bed_index,
      manta_region_bed = manta_region_bed,
      manta_region_bed_index = manta_region_bed_index
  }


  output {
    File inputs_json = RunSVShell.inputs_json
    File outputs_json = RunSVShell.outputs_json
    File final_vcf = RunSVShell.final_vcf
    File final_vcf_idx = RunSVShell.final_vcf_idx
    File pre_cleanup_vcf = RunSVShell.pre_cleanup_vcf
    File pre_cleanup_vcf_idx = RunSVShell.pre_cleanup_vcf_idx
    File stripy_json_output = RunSVShell.stripy_json_output
    File stripy_tsv_output = RunSVShell.stripy_tsv_output
    File stripy_html_output = RunSVShell.stripy_html_output
    File stripy_vcf_output = RunSVShell.stripy_vcf_output
    File metrics_file = RunSVShell.metrics_file
    File qc_file = RunSVShell.qc_file
    File ploidy_matrix = RunSVShell.ploidy_matrix
    File ploidy_plots = RunSVShell.ploidy_plots
    File non_genotyped_unique_depth_calls = RunSVShell.non_genotyped_unique_depth_calls
  }
}

task RunSVShell {
  input {
    String batch
    String sample_id
    File ref_samples_list
    File ref_ped_file
    File genome_file
    File primary_contigs_list
    File primary_contigs_fai
    File reference_fasta
    File reference_index
    File reference_dict
    File autosome_file
    File allosome_file
    File bam_or_cram_file
    File bam_or_cram_index
    File preprocessed_intervals
    File sd_locs_vcf
    File wham_include_list_bed_file
    File reference_bwa_alt
    File reference_bwa_amb
    File reference_bwa_ann
    File reference_bwa_bwt
    File reference_bwa_pac
    File reference_bwa_sa
    Boolean run_vcf_qc
    File wgd_scoring_mask
    Int min_svsize
    File contig_ploidy_model_tar
    Array[File] gcnv_model_tars
    Array[File] ref_pesr_split_files
    Array[File] ref_pesr_split_files_indices
    Array[File] ref_pesr_disc_files
    Array[File] ref_pesr_disc_files_indices
    Array[File] ref_pesr_sd_files
    Array[File] ref_pesr_sd_files_indices
    Array[File] genome_tracks
    Array[File] genome_tracks_indices
    File pesr_exclude_intervals
    File pesr_exclude_intervals_index
    File ref_panel_vcf
    File ref_panel_vcf_index
    File ref_panel_bincov_matrix
    File ref_panel_bincov_matrix_index
    File bin_exclude
    File bin_exclude_index
    File cytobands
    File cytobands_index
    File depth_exclude_list
    File depth_exclude_list_index
    File mei_bed
    File mei_bed_index
    File manta_region_bed
    File manta_region_bed_index
    File HERVK_reference
    File LINE1_reference
    File intron_reference
    File? par_bed
    File rmsk
    File segdups
    Int gcnv_qs_cutoff
    Int ref_copy_number_autosomal_contigs
    File cnmops_exclude_list
    Int matrix_qc_distance
    File? ref_std_manta_vcf_tar
    File? ref_std_scramble_vcf_tar
    File? ref_std_wham_vcf_tar
    File ref_panel_del_bed
    File ref_panel_dup_bed
    Float depth_exclude_overlap_fraction
    Float depth_interval_overlap
    String? depth_clustering_algorithm
    Float pesr_interval_overlap
    String? pesr_clustering_algorithm
    File cutoffs
    File genotyping_rd_table
    File genotyping_pe_table
    File genotyping_sr_table
    Float clean_vcf_min_sr_background_fail_batches
    File clustering_config_part1
    File stratification_config_part1
    File clustering_config_part2
    File stratification_config_part2
    Array[String] clustering_track_names
    Array[File] clustering_track_bed_files
    Int max_shard_size_resolve
    String chr_x
    String chr_y
    File protein_coding_gtf
    File noncoding_bed
    Int annotation_sv_per_shard
    File? external_af_ref_bed
    String? external_af_ref_bed_prefix
    Array[String]? external_af_population
    Int min_pe_cpx
    Int min_pe_ctx
    File gq_recalibrator_model_file
    Array[String] recalibrate_gq_args
    Float no_call_rate_cutoff
    File sl_cutoff_table
    String? sl_filter_args
    File qc_definitions
    File ref_panel_median_cov
    File? outlier_samples_list
    Boolean run_sampleevidence_metrics
    Int genotyping_n_per_split
    Int n_RD_genotype_bins
    Int clean_vcf1b_records_per_shard
    Int clean_vcf5_records_per_shard
    Int clean_vcf_max_shards_per_chrom_clean_vcf_step1
    Int clean_vcf_min_records_per_shard_clean_vcf_step1
    Int clean_vcf_random_seed
    Int clean_vcf_samples_per_clean_vcf_step2_shard
    Int refine_complex_variants_n_per_split
    String sv_shell_docker
    RuntimeAttr? runtime_attr_override

    String sv_shell_docker
    RuntimeAttr? runtime_attr_override
  }

  String final_vcf_filename = sample_id + ".vcf.gz"
  String final_vcf_idx_filename = final_vcf_filename + ".tbi"
  String pre_cleanup_vcf_filename = batch + ".annotated.vcf.gz"
  String pre_cleanup_vcf_idx_filename = pre_cleanup_vcf_filename + ".tbi"
  String stripy_json_filename = sample_id + ".stripy.json"
  String stripy_tsv_filename = sample_id + ".stripy.tsv"
  String stripy_html_filename = sample_id + ".stripy.html"
  String stripy_vcf_filename = sample_id + ".stripy.vcf"
  String metrics_filename = "single_sample." + batch + ".metrics.tsv"
  String qc_filename = "sv_qc." + batch + ".tsv"
  String ploidy_matrix_filename = batch + "_ploidy_matrix.bed.gz"
  String ploidy_plots_filename = batch + "_ploidy_plots.tar.gz"
  String non_genotyped_unique_depth_calls_filename = batch + ".non_genotyped_unique_depth_calls.vcf.gz"

  command <<<
    set -Exeuo pipefail

    export BASE_DIR="${PWD}"
    export SV_SHELL_BASE_DIR="${PWD}/wd"
    export TMPDIR="${PWD}/wd/tmp"
    mkdir -p "${PWD}/wd/tmp"

    touch single_sample_pipeline_inputs.json
    touch single_sample_pipeline_outputs.json
    touch "~{final_vcf_filename}"
    touch "~{final_vcf_idx_filename}"
    touch "~{pre_cleanup_vcf_filename}"
    touch "~{pre_cleanup_vcf_idx_filename}"
    touch "~{stripy_json_filename}"
    touch "~{stripy_tsv_filename}"
    touch "~{stripy_html_filename}"
    touch "~{stripy_vcf_filename}"
    touch "~{metrics_filename}"
    touch "~{qc_filename}"
    touch "~{ploidy_matrix_filename}"
    touch "~{ploidy_plots_filename}"
    touch "~{non_genotyped_unique_depth_calls_filename}"
  >>>

  output {
    File inputs_json = "single_sample_pipeline_inputs.json"
    File outputs_json = "single_sample_pipeline_outputs.json"
    File final_vcf = final_vcf_filename
    File final_vcf_idx = final_vcf_idx_filename
    File pre_cleanup_vcf = pre_cleanup_vcf_filename
    File pre_cleanup_vcf_idx = pre_cleanup_vcf_idx_filename
    File stripy_json_output = stripy_json_filename
    File stripy_tsv_output = stripy_tsv_filename
    File stripy_html_output = stripy_html_filename
    File stripy_vcf_output = stripy_vcf_filename
    File metrics_file = metrics_filename
    File qc_file = qc_filename
    File ploidy_matrix = ploidy_matrix_filename
    File ploidy_plots = ploidy_plots_filename
    File non_genotyped_unique_depth_calls = non_genotyped_unique_depth_calls_filename
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 4,
    mem_gb: 16,
    disk_gb: 400,
    boot_disk_gb: 30,
    preemptible_tries: 0,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_shell_docker
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}