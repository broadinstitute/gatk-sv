version 1.0

import "Structs.wdl"

workflow SVShell {

  input {
    File ref_samples_list
  }

  Array[File] ref_samples = read_lines(ref_samples_list)

  call RunSVShell {
    input:
      ref_samples_list = ref_samples
  }

  output {
    File test = RunSVShell.test
#    File final_vcf = RunSVShell.final_vcf
#    File pre_cleanup_vcf = RunSVShell.pre_cleanup_vcf
#    File stripy_json_output = RunSVShell.stripy_json_output
#    File metrics_file = RunSVShell.metrics_file
#    File ploidy_matrix = RunSVShell.ploidy_matrix
  }
}

task RunSVShell {
  input {
    String batch
    String sample_id
    Array[File] ref_samples_list
    File ref_ped_file
    File genome_file
    File primary_contigs_list
    File primary_contigs_fai
    File reference_fasta
    File reference_index
    File reference_dict
    File ref_panel_vcf
    File autosome_file
    File allosome_file
    File bam_or_cram_file
    File bam_or_cram_index
    File preprocessed_intervals
    File manta_region_bed
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
    File gcnv_model_tars_list
    File ref_panel_bincov_matrix
    File ref_pesr_disc_files_list
    File ref_pesr_split_files_list
    File ref_pesr_sd_files_list
    Int ref_copy_number_autosomal_contigs
    Int gcnv_qs_cutoff
    File cnmops_exclude_list
    Int matrix_qc_distance
    File? ref_std_manta_vcf_tar
    File? ref_std_scramble_vcf_tar
    File? ref_std_wham_vcf_tar
    File ref_panel_del_bed
    File ref_panel_dup_bed
    File depth_exclude_list
    Float depth_exclude_overlap_fraction
    Float depth_interval_overlap
    String? depth_clustering_algorithm
    File pesr_exclude_intervals
    Float pesr_interval_overlap
    String? pesr_clustering_algorithm
    File cutoffs
    File genotyping_rd_table
    File genotyping_pe_table
    File genotyping_sr_table
    File bin_exclude
    Float clean_vcf_min_sr_background_fail_batches
    File clustering_config_part1
    File stratification_config_part1
    File clustering_config_part2
    File stratification_config_part2
    Array[String] clustering_track_names
    Array[File] clustering_track_bed_files
    File cytobands
    File mei_bed
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
    Array[File] genome_tracks
    Float no_call_rate_cutoff
    File sl_cutoff_table
    String? sl_filter_args
    File qc_definitions
    String sv_shell_docker
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -Exeuo pipefail

    export SV_SHELL_BASE_DIR="${PWD}/wd"
    export TMPDIR="${PWD}/wd/tmp"
    mkdir -p "${PWD}/wd/tmp"

    df -h

    jq -n \
      --arg batch "~{batch}" \
      --arg sample_id "~{sample_id}" \
      --arg ref_samples_list "~{write_lines(ref_samples_list)}" \
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
      --arg gcnv_model_tars_list "~{gcnv_model_tars_list}" \
      --arg ref_panel_bincov_matrix "~{ref_panel_bincov_matrix}" \
      --arg ref_pesr_disc_files_list "~{ref_pesr_disc_files_list}" \
      --arg ref_pesr_split_files_list "~{ref_pesr_split_files_list}" \
      --arg ref_pesr_sd_files_list "~{ref_pesr_sd_files_list}" \
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
      '$ARGS.named | with_entries(select(.value != "" and .value != null))' > "${SV_SHELL_BASE_DIR}/single_sample_pipeline_inputs.json"

#    bash single_sample_pipeline.sh \
#      "${SV_SHELL_BASE_DIR}/single_sample_pipeline_inputs.json" \
#      "${SV_SHELL_BASE_DIR}/single_sample_pipeline_outputs.json"

    echo "----------------------"
    echo "${PWD}"
    cp "${SV_SHELL_BASE_DIR}/single_sample_pipeline_inputs.json" .
    ls

  >>>

  output {
    File test = "single_sample_pipeline_inputs.json"
#    Map[String, String] manifest = read_json("pipeline_outputs_manifest.json")
#
#    File final_vcf = manifest["final_vcf"]
#    File pre_cleanup_vcf = manifest["pre_cleanup_vcf"]
#    File stripy_json_output = manifest["stripy_json_output"]
#    File stripy_tsv_output = manifest["stripy_tsv_output"]
#    File stripy_html_output = manifest["stripy_html_output"]
#    File stripy_vcf_output = manifest["stripy_vcf_output"]
#    File metrics_file = manifest["metrics_file"]
#    File qc_file = manifest["qc_file"]
#    File ploidy_matrix = manifest["ploidy_matrix"]
#    File ploidy_plots = manifest["ploidy_plots"]
#    File non_genotyped_unique_depth_calls = manifest["non_genotyped_unique_depth_calls"]
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 4,
    mem_gb: 32,
    disk_gb: 500,
    boot_disk_gb: 30,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_shell_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }
}