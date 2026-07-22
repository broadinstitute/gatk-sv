version 1.0

import "Structs.wdl"

workflow SVShell {
  input {
    String sample_id
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
    File primary_contigs_fai
    Int min_svsize
    String sv_pipeline_docker

    File? dragen_sv_vcf
    File? dragen_sv_vcf_index
    File? dragen_cnv_vcf
    File? dragen_cnv_vcf_index
    File? ref_std_dragen_vcf_tar
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
      sample_id = sample_id,
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
      manta_region_bed_index = manta_region_bed_index,
      primary_contigs_fai = primary_contigs_fai,
      min_svsize = min_svsize,
      dragen_sv_vcf = dragen_sv_vcf,
      dragen_sv_vcf_index = dragen_sv_vcf_index,
      dragen_cnv_vcf = dragen_cnv_vcf,
      dragen_cnv_vcf_index = dragen_cnv_vcf_index,
      ref_std_dragen_vcf_tar = ref_std_dragen_vcf_tar
  }


  if (defined(RunSVShell.scramble_vcf)) {
    call StandardizeVcf as StandardizeScramble {
      input:
        sample_id = sample_id,
        vcf_path = select_first([RunSVShell.scramble_vcf]),
        caller = "scramble",
        contigs_fai = primary_contigs_fai,
        min_size = min_svsize,
        sv_pipeline_docker = sv_pipeline_docker
    }
    call FormatVcfForGatk as FormatScramble {
      input:
        sample_id = sample_id,
        vcf_path = StandardizeScramble.standardized_vcf,
        ploidy_table = RunSVShell.ploidy_table,
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  if (defined(RunSVShell.wham_vcf)) {
    call StandardizeVcf as StandardizeWham {
      input:
        sample_id = sample_id,
        vcf_path = select_first([RunSVShell.wham_vcf]),
        caller = "wham",
        contigs_fai = primary_contigs_fai,
        min_size = min_svsize,
        sv_pipeline_docker = sv_pipeline_docker
    }
    call FormatVcfForGatk as FormatWham {
      input:
        sample_id = sample_id,
        vcf_path = StandardizeWham.standardized_vcf,
        ploidy_table = RunSVShell.ploidy_table,
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  if (defined(dragen_sv_vcf)) {
    call StandardizeVcf as StandardizeDragenSv {
      input:
        sample_id = sample_id,
        vcf_path = select_first([dragen_sv_vcf]),
        caller = "dragen",
        contigs_fai = primary_contigs_fai,
        min_size = min_svsize,
        sv_pipeline_docker = sv_pipeline_docker
    }
    call FormatVcfForGatk as FormatDragenSv {
      input:
        sample_id = sample_id,
        vcf_path = StandardizeDragenSv.standardized_vcf,
        ploidy_table = RunSVShell.ploidy_table,
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  if (defined(dragen_cnv_vcf)) {
    call StandardizeVcf as StandardizeDragenCnv {
      input:
        sample_id = sample_id,
        vcf_path = select_first([dragen_cnv_vcf]),
        caller = "dragen",
        contigs_fai = primary_contigs_fai,
        min_size = min_svsize,
        sv_pipeline_docker = sv_pipeline_docker
    }
    call FormatVcfForGatk as FormatDragenCnv {
      input:
        sample_id = sample_id,
        vcf_path = StandardizeDragenCnv.standardized_vcf,
        ploidy_table = RunSVShell.ploidy_table,
        sv_pipeline_docker = sv_pipeline_docker
    }
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
    File? manta_vcf = RunSVShell.manta_vcf
    File? manta_vcf_idx = RunSVShell.manta_vcf_idx
    File? scramble_vcf = RunSVShell.scramble_vcf
    File? scramble_vcf_idx = RunSVShell.scramble_vcf_idx
    File? wham_vcf = RunSVShell.wham_vcf
    File? wham_vcf_idx = RunSVShell.wham_vcf_idx
    File? scramble_clusters = RunSVShell.scramble_clusters
    File? scramble_table = RunSVShell.scramble_table
    File coverage_counts = RunSVShell.coverage_counts
    File coverage_counts_idx = RunSVShell.coverage_counts_idx
    File merged_dels = RunSVShell.merged_dels
    File merged_dups = RunSVShell.merged_dups
    File ploidy_table = RunSVShell.ploidy_table
    File? scramble_vcf_formatted = FormatScramble.formatted_vcf
    File? scramble_vcf_formatted_index = FormatScramble.formatted_vcf_index
    File? wham_vcf_formatted = FormatWham.formatted_vcf
    File? wham_vcf_formatted_index = FormatWham.formatted_vcf_index
    File? dragen_sv_vcf_formatted = FormatDragenSv.formatted_vcf
    File? dragen_sv_vcf_formatted_index = FormatDragenSv.formatted_vcf_index
    File? dragen_cnv_vcf_formatted = FormatDragenCnv.formatted_vcf
    File? dragen_cnv_vcf_formatted_index = FormatDragenCnv.formatted_vcf_index
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
    File? dragen_sv_vcf
    File? dragen_sv_vcf_index
    File? dragen_cnv_vcf
    File? dragen_cnv_vcf_index
    File? ref_std_dragen_vcf_tar

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
  String manta_vcf_filename = sample_id + ".manta.vcf.gz"
  String manta_vcf_idx_filename = manta_vcf_filename + ".tbi"
  String scramble_vcf_filename = sample_id + ".scramble.vcf.gz"
  String scramble_vcf_idx_filename = scramble_vcf_filename + ".tbi"
  String wham_vcf_filename = sample_id + ".wham.vcf.gz"
  String wham_vcf_idx_filename = wham_vcf_filename + ".tbi"
  String scramble_clusters_filename = sample_id + ".scramble_clusters.tsv.gz"
  String scramble_table_filename = sample_id + ".scramble.tsv.gz"
  String coverage_counts_filename = sample_id + ".counts.tsv.gz"
  String coverage_counts_idx_filename = coverage_counts_filename + ".tbi"
  String merged_dels_filename = batch + ".DEL.bed.gz"
  String merged_dups_filename = batch + ".DUP.bed.gz"
  String ploidy_table_filename = sample_id + "." + batch + ".ploidy.tsv"

  command <<<
    set -Exeuo pipefail

    export BASE_DIR="${PWD}"
    export SV_SHELL_BASE_DIR="${PWD}/wd"
    export TMPDIR="${PWD}/wd/tmp"
    mkdir -p "${PWD}/wd/tmp"

    jq -n \
      --arg batch "~{batch}" \
      --arg sample_id "~{sample_id}" \
      --arg ref_samples_list "~{ref_samples_list}" \
      --arg ref_ped_file "~{ref_ped_file}" \
      --arg genome_file "~{genome_file}" \
      --arg primary_contigs_list "~{primary_contigs_list}" \
      --arg primary_contigs_fai "~{primary_contigs_fai}" \
      --arg reference_fasta "~{reference_fasta}" \
      --arg reference_index "~{reference_index}" \
      --arg reference_dict "~{reference_dict}" \
      --arg ref_panel_vcf "~{ref_panel_vcf}" \
      --arg autosome_file "~{autosome_file}" \
      --arg allosome_file "~{allosome_file}" \
      --arg bam_or_cram_file "~{bam_or_cram_file}" \
      --arg bam_or_cram_index "~{bam_or_cram_index}" \
      --arg preprocessed_intervals "~{preprocessed_intervals}" \
      --arg manta_region_bed "~{manta_region_bed}" \
      --arg sd_locs_vcf "~{sd_locs_vcf}" \
      --arg wham_include_list_bed_file "~{wham_include_list_bed_file}" \
      --arg reference_bwa_alt "~{reference_bwa_alt}" \
      --arg reference_bwa_amb "~{reference_bwa_amb}" \
      --arg reference_bwa_ann "~{reference_bwa_ann}" \
      --arg reference_bwa_bwt "~{reference_bwa_bwt}" \
      --arg reference_bwa_pac "~{reference_bwa_pac}" \
      --arg reference_bwa_sa "~{reference_bwa_sa}" \
      --argjson run_vcf_qc ~{run_vcf_qc} \
      --arg wgd_scoring_mask "~{wgd_scoring_mask}" \
      --argjson min_svsize ~{min_svsize} \
      --arg contig_ploidy_model_tar "~{contig_ploidy_model_tar}" \
      --arg gcnv_model_tars_list "~{write_lines(gcnv_model_tars)}" \
      --arg ref_panel_bincov_matrix "~{ref_panel_bincov_matrix}" \
      --arg ref_pesr_disc_files_list "~{write_lines(ref_pesr_disc_files)}" \
      --arg ref_pesr_split_files_list "~{write_lines(ref_pesr_split_files)}" \
      --arg ref_pesr_sd_files_list "~{write_lines(ref_pesr_sd_files)}" \
      --argjson ref_copy_number_autosomal_contigs ~{ref_copy_number_autosomal_contigs} \
      --argjson gcnv_qs_cutoff ~{gcnv_qs_cutoff} \
      --arg cnmops_exclude_list "~{cnmops_exclude_list}" \
      --argjson matrix_qc_distance ~{matrix_qc_distance} \
      --arg ref_std_manta_vcf_tar "~{ref_std_manta_vcf_tar}" \
      --arg ref_std_scramble_vcf_tar "~{ref_std_scramble_vcf_tar}" \
      --arg ref_std_wham_vcf_tar "~{ref_std_wham_vcf_tar}" \
      --arg ref_panel_del_bed "~{ref_panel_del_bed}" \
      --arg ref_panel_dup_bed "~{ref_panel_dup_bed}" \
      --arg depth_exclude_list "~{depth_exclude_list}" \
      --argjson depth_exclude_overlap_fraction ~{depth_exclude_overlap_fraction} \
      --argjson depth_interval_overlap ~{depth_interval_overlap} \
      --arg depth_clustering_algorithm "~{select_first([depth_clustering_algorithm, ""])}" \
      --arg pesr_exclude_intervals "~{pesr_exclude_intervals}" \
      --argjson pesr_interval_overlap ~{pesr_interval_overlap} \
      --arg pesr_clustering_algorithm "~{select_first([pesr_clustering_algorithm, ""])}" \
      --arg cutoffs "~{cutoffs}" \
      --arg genotyping_rd_table "~{genotyping_rd_table}" \
      --arg genotyping_pe_table "~{genotyping_pe_table}" \
      --arg genotyping_sr_table "~{genotyping_sr_table}" \
      --arg bin_exclude "~{bin_exclude}" \
      --argjson clean_vcf_min_sr_background_fail_batches ~{clean_vcf_min_sr_background_fail_batches} \
      --arg clustering_config_part1 "~{clustering_config_part1}" \
      --arg stratification_config_part1 "~{stratification_config_part1}" \
      --arg clustering_config_part2 "~{clustering_config_part2}" \
      --arg stratification_config_part2 "~{stratification_config_part2}" \
      --argjson clustering_track_names "$(jq -R . < ~{write_lines(clustering_track_names)} | jq -s .)" \
      --argjson clustering_track_bed_files "$(jq -R . < ~{write_lines(clustering_track_bed_files)} | jq -s .)" \
      --arg cytobands "~{cytobands}" \
      --arg mei_bed "~{mei_bed}" \
      --argjson max_shard_size_resolve ~{max_shard_size_resolve} \
      --arg chr_x "~{chr_x}" \
      --arg chr_y "~{chr_y}" \
      --arg protein_coding_gtf "~{protein_coding_gtf}" \
      --arg noncoding_bed "~{noncoding_bed}" \
      --argjson annotation_sv_per_shard ~{annotation_sv_per_shard} \
      --arg external_af_ref_bed "~{select_first([external_af_ref_bed, ""])}" \
      --arg external_af_ref_bed_prefix "~{select_first([external_af_ref_bed_prefix, ""])}" \
      --argjson external_af_population "$(jq -R . < ~{write_lines(select_first([external_af_population, []]))} | jq -s .)" \
      --argjson min_pe_cpx ~{min_pe_cpx} \
      --argjson min_pe_ctx ~{min_pe_ctx} \
      --arg gq_recalibrator_model_file "~{gq_recalibrator_model_file}" \
      --argjson recalibrate_gq_args "$(jq -R . < ~{write_lines(recalibrate_gq_args)} | jq -s .)" \
      --argjson genome_tracks "$(jq -R . < ~{write_lines(genome_tracks)} | jq -s .)" \
      --argjson no_call_rate_cutoff ~{no_call_rate_cutoff} \
      --arg sl_cutoff_table "~{sl_cutoff_table}" \
      --arg sl_filter_args "~{select_first([sl_filter_args, ""])}" \
      --arg qc_definitions "~{qc_definitions}" \
      --arg ref_panel_median_cov "~{ref_panel_median_cov}" \
      --arg manta_region_bed_index "~{manta_region_bed_index}" \
      --arg HERVK_reference "~{HERVK_reference}" \
      --arg LINE1_reference "~{LINE1_reference}" \
      --arg intron_reference "~{intron_reference}" \
      --arg outlier_samples_list "~{outlier_samples_list}" \
      --arg par_bed "~{par_bed}" \
      --argjson run_sampleevidence_metrics ~{run_sampleevidence_metrics} \
      --arg rmsk "~{rmsk}" \
      --arg segdups "~{segdups}" \
      --argjson genotyping_n_per_split ~{genotyping_n_per_split} \
      --argjson n_RD_genotype_bins ~{n_RD_genotype_bins} \
      --argjson clean_vcf1b_records_per_shard ~{clean_vcf1b_records_per_shard} \
      --argjson clean_vcf5_records_per_shard ~{clean_vcf5_records_per_shard} \
      --argjson clean_vcf_max_shards_per_chrom_clean_vcf_step1 ~{clean_vcf_max_shards_per_chrom_clean_vcf_step1} \
      --argjson clean_vcf_min_records_per_shard_clean_vcf_step1 ~{clean_vcf_min_records_per_shard_clean_vcf_step1} \
      --argjson clean_vcf_random_seed ~{clean_vcf_random_seed} \
      --argjson clean_vcf_samples_per_clean_vcf_step2_shard ~{clean_vcf_samples_per_clean_vcf_step2_shard} \
      --argjson "RefineComplexVariants.n_per_split" ~{refine_complex_variants_n_per_split} \
      --arg dragen_sv_vcf "~{select_first([dragen_sv_vcf, ""])}" \
      --arg dragen_sv_vcf_index "~{select_first([dragen_sv_vcf_index, ""])}" \
      --arg dragen_cnv_vcf "~{select_first([dragen_cnv_vcf, ""])}" \
      --arg dragen_cnv_vcf_index "~{select_first([dragen_cnv_vcf_index, ""])}" \
      --arg ref_std_dragen_vcf_tar "~{select_first([ref_std_dragen_vcf_tar, ""])}" \
      '$ARGS.named | with_entries(select(.value != "" and .value != null))' > "${SV_SHELL_BASE_DIR}/single_sample_pipeline_inputs.json"

    bash /opt/sv_shell/single_sample_pipeline.sh \
      "${SV_SHELL_BASE_DIR}/single_sample_pipeline_inputs.json" \
      "${SV_SHELL_BASE_DIR}/single_sample_pipeline_outputs.json"

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
    touch "~{manta_vcf_filename}"
    touch "~{manta_vcf_idx_filename}"
    touch "~{scramble_vcf_filename}"
    touch "~{scramble_vcf_idx_filename}"
    touch "~{wham_vcf_filename}"
    touch "~{wham_vcf_idx_filename}"
    touch "~{scramble_clusters_filename}"
    touch "~{scramble_table_filename}"
    touch "~{coverage_counts_filename}"
    touch "~{coverage_counts_idx_filename}"
    touch "~{merged_dels_filename}"
    touch "~{merged_dups_filename}"
    touch "~{ploidy_table_filename}"

    cp "${SV_SHELL_BASE_DIR}/single_sample_pipeline_inputs.json" "${BASE_DIR}/"
    cp "${SV_SHELL_BASE_DIR}/single_sample_pipeline_outputs.json" "${BASE_DIR}/"

    final_vcf_path=$(jq -r '.final_vcf' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    mv "${final_vcf_path}" "${BASE_DIR}/~{final_vcf_filename}"
    mv "${final_vcf_path}.tbi" "${BASE_DIR}/~{final_vcf_idx_filename}"

    pre_cleanup_vcf_path=$(jq -r '.pre_cleanup_vcf' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    mv "${pre_cleanup_vcf_path}" "${BASE_DIR}/~{pre_cleanup_vcf_filename}"
    mv "${pre_cleanup_vcf_path}.tbi" "${BASE_DIR}/~{pre_cleanup_vcf_idx_filename}"

    stripy_json_path=$(jq -r '.stripy_json_output' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    mv "${stripy_json_path}" "${BASE_DIR}/~{stripy_json_filename}"

    stripy_tsv_path=$(jq -r '.stripy_tsv_output' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    mv "${stripy_tsv_path}" "${BASE_DIR}/~{stripy_tsv_filename}"

    stripy_html_path=$(jq -r '.stripy_html_output' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    mv "${stripy_html_path}" "${BASE_DIR}/~{stripy_html_filename}"

    stripy_vcf_path=$(jq -r '.stripy_vcf_output' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    mv "${stripy_vcf_path}" "${BASE_DIR}/~{stripy_vcf_filename}"

    metrics_path=$(jq -r '.metrics_file' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    mv "${metrics_path}" "${BASE_DIR}/~{metrics_filename}"

    qc_path=$(jq -r '.qc_file' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    mv "${qc_path}" "${BASE_DIR}/~{qc_filename}"

    ploidy_matrix_path=$(jq -r '.ploidy_matrix' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    mv "${ploidy_matrix_path}" "${BASE_DIR}/~{ploidy_matrix_filename}"

    ploidy_plots_path=$(jq -r '.ploidy_plots' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    mv "${ploidy_plots_path}" "${BASE_DIR}/~{ploidy_plots_filename}"

    non_genotyped_path=$(jq -r '.non_genotyped_unique_depth_calls' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    mv "${non_genotyped_path}" "${BASE_DIR}/~{non_genotyped_unique_depth_calls_filename}"

    manta_vcf_path=$(jq -r '.manta_vcf' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    if [[ -n "${manta_vcf_path}" && "${manta_vcf_path}" != "null" ]]; then
      mv "${manta_vcf_path}" "${BASE_DIR}/~{manta_vcf_filename}"
      mv "${manta_vcf_path}.tbi" "${BASE_DIR}/~{manta_vcf_idx_filename}"
    fi

    scramble_vcf_path=$(jq -r '.scramble_vcf' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    if [[ -n "${scramble_vcf_path}" && "${scramble_vcf_path}" != "null" ]]; then
      mv "${scramble_vcf_path}" "${BASE_DIR}/~{scramble_vcf_filename}"
      mv "${scramble_vcf_path}.tbi" "${BASE_DIR}/~{scramble_vcf_idx_filename}"
    fi

    wham_vcf_path=$(jq -r '.wham_vcf' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    if [[ -n "${wham_vcf_path}" && "${wham_vcf_path}" != "null" ]]; then
      mv "${wham_vcf_path}" "${BASE_DIR}/~{wham_vcf_filename}"
      mv "${wham_vcf_path}.tbi" "${BASE_DIR}/~{wham_vcf_idx_filename}"
    fi

    scramble_clusters_path=$(jq -r '.scramble_clusters' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    if [[ -n "${scramble_clusters_path}" && "${scramble_clusters_path}" != "null" ]]; then
      mv "${scramble_clusters_path}" "${BASE_DIR}/~{scramble_clusters_filename}"
    fi

    scramble_table_path=$(jq -r '.scramble_table' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    if [[ -n "${scramble_table_path}" && "${scramble_table_path}" != "null" ]]; then
      mv "${scramble_table_path}" "${BASE_DIR}/~{scramble_table_filename}"
    fi

    coverage_counts_path=$(jq -r '.coverage_counts' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    mv "${coverage_counts_path}" "${BASE_DIR}/~{coverage_counts_filename}"
    coverage_counts_idx_path=$(jq -r '.coverage_counts_idx' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    mv "${coverage_counts_idx_path}" "${BASE_DIR}/~{coverage_counts_idx_filename}"

    merged_dels_path=$(jq -r '.merged_dels' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    mv "${merged_dels_path}" "${BASE_DIR}/~{merged_dels_filename}"

    merged_dups_path=$(jq -r '.merged_dups' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    mv "${merged_dups_path}" "${BASE_DIR}/~{merged_dups_filename}"

    ploidy_table_path=$(jq -r '.ploidy_table' "${BASE_DIR}/single_sample_pipeline_outputs.json")
    mv "${ploidy_table_path}" "${BASE_DIR}/~{ploidy_table_filename}"

    df -h
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
    File? manta_vcf = manta_vcf_filename
    File? manta_vcf_idx = manta_vcf_idx_filename
    File? scramble_vcf = scramble_vcf_filename
    File? scramble_vcf_idx = scramble_vcf_idx_filename
    File? wham_vcf = wham_vcf_filename
    File? wham_vcf_idx = wham_vcf_idx_filename
    File? scramble_clusters = scramble_clusters_filename
    File? scramble_table = scramble_table_filename
    File coverage_counts = coverage_counts_filename
    File coverage_counts_idx = coverage_counts_idx_filename
    File merged_dels = merged_dels_filename
    File merged_dups = merged_dups_filename
    File ploidy_table = ploidy_table_filename
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

task StandardizeVcf {
  input {
    String sample_id
    File vcf_path
    String caller
    File contigs_fai
    Int min_size
    String sv_pipeline_docker
  }

  command <<<
    set -eu -o pipefail

    svtk standardize \
      --sample-names ~{sample_id} \
      --contigs ~{contigs_fai} \
      --min-size ~{min_size} \
      ~{vcf_path} \
      ~{sample_id}.std.vcf.gz \
      ~{caller}

    tabix -p vcf ~{sample_id}.std.vcf.gz
  >>>

  output {
    File standardized_vcf = "~{sample_id}.std.vcf.gz"
    File standardized_vcf_index = "~{sample_id}.std.vcf.gz.tbi"
  }

  runtime {
    cpu: 1
    memory: "2 GiB"
    disks: "local-disk 20 HDD"
    docker: sv_pipeline_docker
  }
}

task FormatVcfForGatk {
  input {
    String sample_id
    File vcf_path
    File ploidy_table
    String sv_pipeline_docker
  }

  command <<<
    set -eu -o pipefail

    python /opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py \
      --vcf ~{vcf_path} \
      --out ~{sample_id}.fmt.vcf.gz \
      --ploidy-table ~{ploidy_table} \
      --fix-end

    tabix -p vcf ~{sample_id}.fmt.vcf.gz
  >>>

  output {
    File formatted_vcf = "~{sample_id}.fmt.vcf.gz"
    File formatted_vcf_index = "~{sample_id}.fmt.vcf.gz.tbi"
  }

  runtime {
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk 20 SSD"
    docker: sv_pipeline_docker
  }
}