version 1.0

import "Structs.wdl"
import "BatchEvidenceMerging.wdl" as bem
import "CNMOPS.wdl" as cnmops
import "CollectCoverage.wdl" as cov
import "DepthPreprocessing.wdl" as dpn
import "MakeBincovMatrix.wdl" as mbm
import "MatrixQC.wdl" as mqc
import "MedianCov.wdl" as mc
import "GatherBatchEvidenceMetrics.wdl" as metrics
import "PESRPreprocessing.wdl" as pp
import "GermlineCNVCase.wdl" as gcnv
import "PloidyEstimation.wdl" as pe
import "TinyResolve.wdl" as tiny
import "Utils.wdl" as util

# Batch-level workflow:
#   - Merge sample evidence data into a single batch
#   - Run cnMOPS
#   - Run gCNV
#   - Run MedianCoverage

workflow GatherBatchEvidence {
  input {
    # Batch info
    String batch
    Array[String] samples
    Array[String]? ref_panel_samples

    # Optional QC tasks
    Boolean run_matrix_qc

    # Global files
    File ped_file
    File genome_file
    File primary_contigs_fai            # .fai file of included contigs
    File ref_dict

    # PE/SR/BAF/bincov files
    # If neither SD_files nor ref_panel_SD_files is present, BAF_files must be supplied
    # If BAF_files is absent, SD_files and/or ref_panel_SD_files and sd_locs_vcf must be supplied
    Array[File] counts
    File? ref_panel_bincov_matrix
    File? bincov_matrix
    File? bincov_matrix_index
    Boolean subset_primary_contigs = false  # PE/SR/BAF files will be subsetted to primary contigs only (for legacy files with bad sorting)
    Boolean rename_samples = false  # Rename samples in PE/SR/BAF to IDs in the "samples" array (always done for RD)
    Array[File?]? BAF_files         # Required for MatrixQC
    Array[File] PE_files
    Array[File]? ref_panel_PE_files
    Array[File] SR_files
    Array[File]? ref_panel_SR_files
    Array[File]? SD_files	# required unless BAF_files or ref_panel_SD_files is supplied
    Array[File]? ref_panel_SD_files	# required unless BAF_files or SD_files is supplied
    File? sd_locs_vcf	# must be same sd_locs_vcf that was presented to GatherSampleEvidence

    # Condense read counts
    Int? min_interval_size
    Int? max_interval_size

    # gCNV inputs
    File contig_ploidy_model_tar
    Array[File] gcnv_model_tars

    File? gatk4_jar_override
    Float? gcnv_p_alt
    Float? gcnv_cnv_coherence_length
    Int? gcnv_max_copy_number

    Float? gcnv_mapping_error_rate
    Float? gcnv_sample_psi_scale
    Float? gcnv_depth_correction_tau
    String? gcnv_copy_number_posterior_expectation_mode
    Int? gcnv_active_class_padding_hybrid_mode

    Float? gcnv_learning_rate
    Float? gcnv_adamax_beta_1
    Float? gcnv_adamax_beta_2
    Int? gcnv_log_emission_samples_per_round
    Float? gcnv_log_emission_sampling_median_rel_error
    Int? gcnv_log_emission_sampling_rounds
    Int? gcnv_max_advi_iter_first_epoch
    Int? gcnv_max_advi_iter_subsequent_epochs
    Int? gcnv_min_training_epochs
    Int? gcnv_max_training_epochs
    Float? gcnv_initial_temperature
    Int? gcnv_num_thermal_advi_iters
    Int? gcnv_convergence_snr_averaging_window
    Float? gcnv_convergence_snr_trigger_threshold
    Int? gcnv_convergence_snr_countdown_window
    Int? gcnv_max_calling_iters
    Float? gcnv_caller_update_convergence_threshold
    Float? gcnv_caller_internal_admixing_rate
    Float? gcnv_caller_external_admixing_rate
    Boolean? gcnv_disable_annealing

    Float? ploidy_sample_psi_scale
    Int ref_copy_number_autosomal_contigs
    Array[String]? allosomal_contigs

    Boolean run_ploidy = false

    # Option to add first sample to the ped file (for single sample mode); run_ploidy must be true
    Boolean append_first_sample_to_ped = false

    Int gcnv_qs_cutoff              # QS filtering cutoff
    Float? defragment_max_dist

    # SV tool calls
    Array[File]? manta_vcfs        # Manta VCF
    Array[File]? melt_vcfs         # Melt VCF
    Array[File]? scramble_vcfs     # Scramble VCF
    Array[File]? wham_vcfs         # Wham VCF
    Int min_svsize                  # Minimum SV length to include

    # CNMops files
    File cnmops_chrom_file
    File cnmops_exclude_list
    File cnmops_allo_file
    Int? cnmops_large_min_size      # minimum size call to be detected by CNMOPS running in large mode

    # Resolve files 
    File cytoband
    File mei_bed
    # QC files
    Int matrix_qc_distance

    # Module metrics parameters
    # Run module metrics workflow at the end - off by default for GatherBatchEvidence because of runtime/expense
    Boolean? run_module_metrics
    String? sv_pipeline_base_docker  # required if run_module_metrics = true
    File? primary_contigs_list  # required if run_module_metrics = true

    # baseline files are optional for metrics workflow
    # run ClusterBatch for vcf metrics
    File? baseline_merged_dels
    File? baseline_merged_dups
    File? baseline_median_cov
    
    # Runtime parameters
    String sv_base_mini_docker
    String sv_base_docker
    String sv_pipeline_docker
    String sv_pipeline_qc_docker
    String linux_docker
    String condense_counts_docker
    String gatk_docker
    String? gcnv_gatk_docker
    String cnmops_docker

    RuntimeAttr? median_cov_runtime_attr        # Memory ignored, use median_cov_mem_gb_per_sample
    Float? median_cov_mem_gb_per_sample

    RuntimeAttr? evidence_merging_bincov_runtime_attr
    RuntimeAttr? runtime_attr_bem

    RuntimeAttr? cnmops_sample10_runtime_attr
    RuntimeAttr? cnmops_sample3_runtime_attr

    RuntimeAttr? ploidy_score_runtime_attr
    RuntimeAttr? ploidy_build_runtime_attr
    RuntimeAttr? runtime_attr_subset_ped
    RuntimeAttr? runtime_attr_validate_ped
    RuntimeAttr? add_sample_to_ped_runtime_attr
    RuntimeAttr? condense_counts_runtime_attr
    RuntimeAttr? preprocess_calls_runtime_attr
    RuntimeAttr? depth_merge_set_runtime_attr
    RuntimeAttr? depth_merge_sample_runtime_attr
    RuntimeAttr? cnmops_ped_runtime_attr
    RuntimeAttr? cnmops_clean_runtime_attr
    RuntimeAttr? matrix_qc_pesrbaf_runtime_attr
    RuntimeAttr? matrix_qc_rd_runtime_attr
    RuntimeAttr? runtime_attr_tiny_untar
    RuntimeAttr? runtime_attr_tiny_resolve

    RuntimeAttr? runtime_attr_ploidy
    RuntimeAttr? runtime_attr_case
    RuntimeAttr? runtime_attr_postprocess
    RuntimeAttr? runtime_attr_explode
  }

  Array[String] all_samples = flatten(select_all([samples, ref_panel_samples]))
  Array[File] all_PE_files = flatten(select_all([PE_files, ref_panel_PE_files]))
  Array[File] all_SR_files = flatten(select_all([SR_files, ref_panel_SR_files]))
  Array[File] all_SD_files = flatten(select_all([SD_files, ref_panel_SD_files]))

  Array[String] samples_batch = select_first([ref_panel_samples, samples])
  call util.ValidatePedFile {
    input:
      ped_file = ped_file,
      sample_list = write_lines(samples_batch),
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_validate_ped
  }

  call util.SubsetPedFile {
    input:
      ped_file = ValidatePedFile.output_ped,
      sample_list = write_lines(samples_batch),
      subset_name = batch,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_subset_ped
  }

  if (append_first_sample_to_ped) {
    call AddCaseSampleToPed {
      input:
        ref_ped_file = SubsetPedFile.ped_subset_file,
        ploidy_plots = select_first([Ploidy.ploidy_plots]),
        sample_id = samples[0],
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = add_sample_to_ped_runtime_attr
    }
  }

  Float median_cov_mem_gb = select_first([median_cov_mem_gb_per_sample, 0.5]) * length(all_samples) + 7.5
  call mc.MedianCov as MedianCov {
    input:
      bincov_matrix = merged_bincov_,
      cohort_id = batch,
      sv_pipeline_qc_docker = sv_pipeline_qc_docker,
      runtime_attr = median_cov_runtime_attr,
      mem_gb_override = median_cov_mem_gb
  }

  Boolean run_module_metrics_ = if defined(run_module_metrics) then select_first([run_module_metrics]) else false
  if (run_module_metrics_) {
    call metrics.GatherBatchEvidenceMetrics {
      input:
        name = batch,
        samples = samples,
        merged_BAF = BatchEvidenceMerging.merged_BAF,
        merged_SR = BatchEvidenceMerging.merged_SR,
        merged_PE = BatchEvidenceMerging.merged_PE,
        merged_bincov = merged_bincov_,
        merged_dels = MergeDepth.del,
        merged_dups = MergeDepth.dup,
        median_cov = MedianCov.medianCov,
        baseline_merged_dels = baseline_merged_dels,
        baseline_merged_dups = baseline_merged_dups,
        baseline_median_cov = baseline_median_cov,
        contig_list = select_first([primary_contigs_list]),
        sv_pipeline_base_docker = select_first([sv_pipeline_base_docker]),
        linux_docker = linux_docker
    }
  }

  output {
    File median_cov = MedianCov.medianCov
  }
}

task AddCaseSampleToPed {
  input {
    File ref_ped_file
    File ploidy_plots
    String sample_id
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File combined_ped_file = "combined_ped_file.ped"
  }

  command <<<

    set -euo pipefail

    tar xzf ~{ploidy_plots} -C .
    RECORD=$(gunzip -c ploidy_est/sample_sex_assignments.txt.gz | { grep -w "^~{sample_id}" || true; })
    if [ -z "$RECORD" ]; then
      >&2 echo "Error: Sample ~{sample_id} not found in ploidy calls"
      exit 1
    fi
    SEX=$(echo "$RECORD" | cut -f2)

    awk -v sample=~{sample_id} '$2 == sample { print "ERROR: A sample with the name "sample" is already present in the ped file." > "/dev/stderr"; exit 1; }' < ~{ref_ped_file}
    awk -v sample=~{sample_id} -v sex=$SEX '{print} END {OFS="\t"; print "case_sample",sample,"0","0",sex,"1" }' < ~{ref_ped_file} > combined_ped_file.ped
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
