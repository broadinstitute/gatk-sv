version 1.0

import "GatherSampleEvidenceBatch.wdl" as sampleevidence
import "EvidenceQC.wdl" as evidenceqc
import "GATKSVPipelinePhase1.wdl" as phase1
import "GenotypeBatch.wdl" as genotypebatch
import "RegenotypeCNVs.wdl" as regenocnvs
import "MakeCohortVcf.wdl" as makecohortvcf
import "Utils.wdl" as utils
import "Structs.wdl"
import "TestUtils.wdl" as tu

# GATK SV Pipeline batch mode
# Runs GatherSampleEvidence, EvidenceQC, GatherBatchEvidence, ClusterBatch, GenerateBatchMetrics, FilterBatch, GenotypeBatch, RegenotypeCNVs,
# and MakeCohortVcf (CombineBatches, ResolveComplexVariants, GenotypeComplexVariants, and GenotypeComplexVariants)

workflow GATKSVPipelineBatch {
  input {
    # Batch data
    String name
    Array[String] samples

    # Required unless caller and evidence outputs are provided (below)
    Array[File]? bam_or_cram_files
    Array[File]? bam_or_cram_indexes

    # Optionally provide calls and evidence (override caller flags below)
    Array[File]? counts_files_input
    Array[File]? pe_files_input
    Array[File]? sr_files_input
    Array[File]? sd_files_input
    Array[File?]? baf_files_input
    Array[File]? manta_vcfs_input
    Array[File]? melt_vcfs_input
    Array[File]? scramble_vcfs_input
    Array[File]? wham_vcfs_input

    # Enable different callers
    Boolean use_manta = true
    Boolean use_melt = true
    Boolean use_scramble = true
    Boolean use_wham = true

    # Merge contig vcfs at each stage of MakeCohortVcf for QC
    Boolean makecohortvcf_merge_cluster_vcfs = false
    Boolean makecohortvcf_merge_complex_resolve_vcfs = false
    Boolean makecohortvcf_merge_complex_genotype_vcfs = false

    # Global files
    File ped_file
    File genome_file
    File primary_contigs_list
    File primary_contigs_fai
    File reference_fasta
    File reference_index    # Index (.fai), must be in same dir as fasta
    File reference_dict     # Dictionary (.dict), must be in same dir as fasta
    File autosome_file      # fai of autosomal contigs
    File allosome_file      # fai of allosomal contigs
    String chr_x
    String chr_y

    # gCNV
    File contig_ploidy_model_tar
    Array[File] gcnv_model_tars

    # PlotSVCountsPerSample metrics from ClusterBatch in GATKSVPipelinePhase1
    Int? N_IQR_cutoff_plotting
    
    File? outlier_cutoff_table
    File qc_definitions

    # Run module metrics - all modules on by default for batch WDL
    Boolean? run_sampleevidence_metrics
    Boolean? run_batchevidence_metrics = true  # GatherBatchEvidenceMetrics is off by default standalone but on for batch WDL
    Boolean? run_clusterbatch_metrics
    Boolean? run_batchmetrics_metrics
    Boolean? run_filterbatch_metrics
    Boolean? run_genotypebatch_metrics
    Boolean? run_makecohortvcf_metrics

    File? baseline_sampleevidence_metrics
    File? baseline_batchevidence_metrics
    File? baseline_clusterbatch_metrics
    File? baseline_batchmetrics_metrics
    File? baseline_filterbatch_metrics
    File? baseline_genotypebatch_metrics
    File? baseline_makecohortvcf_metrics

    String sv_base_mini_docker
    String sv_base_docker
    String sv_pipeline_docker
    String sv_pipeline_hail_docker
    String sv_pipeline_updates_docker
    String sv_pipeline_rdtest_docker
    String sv_pipeline_base_docker
    String sv_pipeline_qc_docker
    String linux_docker
    String cnmops_docker
    String gatk_docker
    String? gatk_docker_pesr_override
    String? gcnv_gatk_docker
    String condense_counts_docker
    String genomes_in_the_cloud_docker
    String samtools_cloud_docker
    String? manta_docker
    String? melt_docker
    String? scramble_docker
    String? wham_docker
    String cloud_sdk_docker

    # Batch metrics
    RuntimeAttr? runtime_attr_cat_metrics
    RuntimeAttr? runtime_attr_plot_metrics

    # Do not use
    Array[File]? NONE_ARRAY_
    String? NONE_STRING_
  }

  Boolean collect_coverage_ = !defined(counts_files_input)
  Boolean collect_pesr_ = !(defined(pe_files_input) && defined(sr_files_input) && defined(sd_files_input))

  String? manta_docker_ = if (!defined(manta_vcfs_input) && use_manta) then manta_docker else NONE_STRING_
  String? melt_docker_ = if (!defined(melt_vcfs_input) && use_melt) then melt_docker else NONE_STRING_
  String? scramble_docker_ = if (!defined(scramble_vcfs_input) && use_scramble) then scramble_docker else NONE_STRING_
  String? wham_docker_ = if (!defined(wham_vcfs_input) && use_wham) then wham_docker else NONE_STRING_

  Boolean run_sampleevidence = collect_coverage_ || collect_pesr_ || defined(manta_docker_) || defined(melt_docker_) || defined(scramble_docker_) || defined(wham_docker_)

  if (run_sampleevidence) {
    call sampleevidence.GatherSampleEvidenceBatch {
      input:
        bam_or_cram_files=select_first([bam_or_cram_files]),
        bam_or_cram_indexes=bam_or_cram_indexes,
        collect_coverage=collect_coverage_,
        collect_pesr=collect_pesr_,
        sample_ids=samples,
        primary_contigs_list=primary_contigs_list,
        reference_fasta=reference_fasta,
        reference_index=reference_index,
        reference_dict=reference_dict,
        run_module_metrics = run_sampleevidence_metrics,
        primary_contigs_fai = primary_contigs_fai,
        batch = name,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        linux_docker = linux_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        manta_docker=manta_docker_,
        melt_docker=melt_docker_,
        scramble_docker=scramble_docker_,
        wham_docker=wham_docker_,
        gatk_docker=gatk_docker,
        gatk_docker_pesr_override = gatk_docker_pesr_override,
        genomes_in_the_cloud_docker=genomes_in_the_cloud_docker,
        samtools_cloud_docker=samtools_cloud_docker,
        cloud_sdk_docker = cloud_sdk_docker
    }
  }

  Array[File] counts_files_ = if collect_coverage_ then select_all(select_first([GatherSampleEvidenceBatch.coverage_counts])) else select_first([counts_files_input])
  Array[File] pe_files_ = if collect_pesr_ then select_all(select_first([GatherSampleEvidenceBatch.pesr_disc])) else select_first([pe_files_input])
  Array[File] sr_files_ = if collect_pesr_ then select_all(select_first([GatherSampleEvidenceBatch.pesr_split])) else select_first([sr_files_input])
  Array[File] sd_files_ = if collect_pesr_ then select_all(select_first([GatherSampleEvidenceBatch.pesr_sd])) else select_first([sd_files_input])

  if (use_manta) {
    Array[File] manta_vcfs_ = if defined(manta_vcfs_input) then select_first([manta_vcfs_input]) else select_all(select_first([GatherSampleEvidenceBatch.manta_vcf]))
  }
  if (use_melt) {
    Array[File] melt_vcfs_ = if defined(melt_vcfs_input) then select_first([melt_vcfs_input]) else select_all(select_first([GatherSampleEvidenceBatch.melt_vcf]))
  }
  if (use_scramble) {
    Array[File] scramble_vcfs_ = if defined(scramble_vcfs_input) then select_first([scramble_vcfs_input]) else select_all(select_first([GatherSampleEvidenceBatch.scramble_vcf]))
  }
  if (use_wham) {
    Array[File] wham_vcfs_ = if defined(wham_vcfs_input) then select_first([wham_vcfs_input]) else select_all(select_first([GatherSampleEvidenceBatch.wham_vcf]))
  }

  call evidenceqc.EvidenceQC as EvidenceQC {
    input:
      batch=name,
      samples=samples,
      genome_file=genome_file,
      counts=counts_files_,
      run_ploidy = false,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_qc_docker=sv_pipeline_qc_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_base_docker=sv_base_docker
  }

  call phase1.GATKSVPipelinePhase1 {
    input:
      batch=name,
      samples=samples,
      ped_file=ped_file,
      genome_file=genome_file,
      contigs=primary_contigs_fai,
      reference_fasta=reference_fasta,
      reference_index=reference_index,
      reference_dict=reference_dict,
      chr_x=chr_x,
      chr_y=chr_y,
      contig_ploidy_model_tar=contig_ploidy_model_tar,
      gcnv_model_tars=gcnv_model_tars,
      BAF_files=baf_files_input,
      counts=counts_files_,
      bincov_matrix=EvidenceQC.bincov_matrix,
      bincov_matrix_index=EvidenceQC.bincov_matrix_index,
      N_IQR_cutoff_plotting = N_IQR_cutoff_plotting,
      PE_files=pe_files_,
      SR_files=sr_files_,
      SD_files=sd_files_,
      manta_vcfs=manta_vcfs_,
      melt_vcfs=melt_vcfs_,
      scramble_vcfs=scramble_vcfs_,
      wham_vcfs=wham_vcfs_,

      cnmops_chrom_file=autosome_file,
      cnmops_allo_file=allosome_file,
      allosome_contigs=allosome_file,
      autosome_contigs=autosome_file,
      run_batchevidence_metrics = run_batchevidence_metrics,
      run_clusterbatch_metrics = run_clusterbatch_metrics,
      run_batchmetrics_metrics = run_batchmetrics_metrics,
      run_filterbatch_metrics = run_filterbatch_metrics,
      primary_contigs_list = primary_contigs_list,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_base_docker=sv_base_docker,
      sv_pipeline_base_docker=sv_pipeline_base_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
      sv_pipeline_qc_docker=sv_pipeline_qc_docker,
      linux_docker=linux_docker,
      cnmops_docker=cnmops_docker,
      gatk_docker=gatk_docker,
      gcnv_gatk_docker=gcnv_gatk_docker,
      condense_counts_docker=condense_counts_docker
  }

  call genotypebatch.GenotypeBatch as GenotypeBatch {
    input:
      batch_pesr_vcf=select_first([GATKSVPipelinePhase1.filtered_pesr_vcf]),
      batch_depth_vcf=select_first([GATKSVPipelinePhase1.filtered_depth_vcf]),
      cohort_pesr_vcf=select_first([GATKSVPipelinePhase1.filtered_pesr_vcf]),
      cohort_depth_vcf=select_first([GATKSVPipelinePhase1.filtered_depth_vcf]),
      batch=name,
      rf_cutoffs=GATKSVPipelinePhase1.cutoffs,
      medianfile=GATKSVPipelinePhase1.median_cov,
      coveragefile=GATKSVPipelinePhase1.merged_bincov,
      coveragefile_index=GATKSVPipelinePhase1.merged_bincov_index,
      discfile=GATKSVPipelinePhase1.merged_PE,
      discfile_index=GATKSVPipelinePhase1.merged_PE_index,
      splitfile=GATKSVPipelinePhase1.merged_SR,
      splitfile_index=GATKSVPipelinePhase1.merged_SR_index,
      ref_dict=reference_dict,
      run_module_metrics = run_genotypebatch_metrics,
      primary_contigs_list = primary_contigs_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
      linux_docker=linux_docker
  }

  call regenocnvs.RegenotypeCNVs as RegenotypeCNVs {
    input:
      depth_vcfs=[GenotypeBatch.genotyped_depth_vcf],
      batch_depth_vcfs=[select_first([GATKSVPipelinePhase1.filtered_depth_vcf])],
      cohort_depth_vcf=select_first([GATKSVPipelinePhase1.filtered_depth_vcf]),
      batches=[name],
      cohort=name,
      medianfiles=[GATKSVPipelinePhase1.median_cov],
      coveragefiles=[GATKSVPipelinePhase1.merged_bincov],
      coveragefile_idxs=[GATKSVPipelinePhase1.merged_bincov_index],
      RD_depth_sepcutoffs=[select_first([GenotypeBatch.trained_genotype_depth_depth_sepcutoff])],
      contig_list=primary_contigs_list,
      regeno_coverage_medians=[GenotypeBatch.regeno_coverage_medians],
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
      sv_pipeline_base_docker=sv_pipeline_base_docker
  }
  

  call makecohortvcf.MakeCohortVcf as MakeCohortVcf {
    input:
      merge_cluster_vcfs = makecohortvcf_merge_cluster_vcfs,
      merge_complex_resolve_vcfs = makecohortvcf_merge_complex_resolve_vcfs,
      merge_complex_genotype_vcfs = makecohortvcf_merge_complex_genotype_vcfs,
      raw_sr_bothside_pass_files=[GenotypeBatch.sr_bothside_pass],
      raw_sr_background_fail_files=[GenotypeBatch.sr_background_fail],
      ped_file=ped_file,
      pesr_vcfs=[GenotypeBatch.genotyped_pesr_vcf],
      depth_vcfs=RegenotypeCNVs.regenotyped_depth_vcfs,
      contig_list=primary_contigs_fai,
      allosome_fai=allosome_file,
      ref_dict=reference_dict,
      chr_x=chr_x,
      chr_y=chr_y,
      disc_files=[GATKSVPipelinePhase1.merged_PE],
      bincov_files=[GATKSVPipelinePhase1.merged_bincov],
      cohort_name=name,
      rf_cutoff_files=[GATKSVPipelinePhase1.cutoffs],
      batches=[name],
      depth_gt_rd_sep_files=[select_first([GenotypeBatch.trained_genotype_depth_depth_sepcutoff])],
      median_coverage_files=[GATKSVPipelinePhase1.median_cov],
      run_module_metrics = run_makecohortvcf_metrics,
      primary_contigs_list = primary_contigs_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      linux_docker=linux_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_hail_docker=sv_pipeline_hail_docker,
      sv_pipeline_updates_docker=sv_pipeline_updates_docker,
      sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
      sv_pipeline_qc_docker=sv_pipeline_qc_docker,
      sv_base_mini_docker=sv_base_mini_docker
  }

  call tu.CatMetrics as CatBatchMetrics {
      input:
        prefix = "batch_sv." + name,
        metric_files = select_all([GatherSampleEvidenceBatch.metrics_file_sampleevidence, GATKSVPipelinePhase1.metrics_file_batchevidence, GATKSVPipelinePhase1.metrics_file_clusterbatch, GATKSVPipelinePhase1.metrics_file_batchmetrics, GATKSVPipelinePhase1.metrics_file_filterbatch, GenotypeBatch.metrics_file_genotypebatch, MakeCohortVcf.metrics_file_makecohortvcf]),
        linux_docker = linux_docker,
        runtime_attr_override = runtime_attr_cat_metrics
    }

  Array[File] defined_baseline_metrics = select_all([baseline_sampleevidence_metrics, baseline_batchevidence_metrics, baseline_clusterbatch_metrics, baseline_batchmetrics_metrics, baseline_filterbatch_metrics, baseline_genotypebatch_metrics, baseline_makecohortvcf_metrics])
  if (length(defined_baseline_metrics) > 0) {
    call tu.CatMetrics as CatBaselineMetrics {
      input:
        prefix = "baseline." + name,
        metric_files = defined_baseline_metrics,
        linux_docker = linux_docker,
        runtime_attr_override = runtime_attr_cat_metrics
    }
    call tu.PlotMetrics {
      input:
        name = name,
        samples = samples,
        test_metrics = CatBatchMetrics.out,
        base_metrics = CatBaselineMetrics.out,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_plot_metrics
    }
  }

  call utils.RunQC as BatchQC {
    input:
      name = name,
      metrics = CatBatchMetrics.out,
      qc_definitions = qc_definitions,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  scatter (i in range(length(samples))) {
    File pe_files_index_ = pe_files_[i] + ".tbi"
    File sr_files_index_ = sr_files_[i] + ".tbi"
  }

  if (defined(manta_vcfs_)) {
    scatter (i in range(length(samples))) {
      File manta_vcfs_index_ = select_first([manta_vcfs_])[i] + ".tbi"
    }
  }
  if (defined(melt_vcfs_)) {
    scatter (i in range(length(samples))) {
      File melt_vcfs_index_ = select_first([melt_vcfs_])[i] + ".tbi"
    }
  }
  if (defined(wham_vcfs_)) {
    scatter (i in range(length(samples))) {
      File wham_vcfs_index_ = select_first([wham_vcfs_])[i] + ".tbi"
    }
  }

  output {
    File clean_vcf = MakeCohortVcf.vcf
    File clean_vcf_index = MakeCohortVcf.vcf_index
    File metrics_file_batch = CatBatchMetrics.out
    File qc_file = BatchQC.out
    File master_vcf_qc = MakeCohortVcf.vcf_qc
    File? metrics_file_makecohortvcf = MakeCohortVcf.metrics_file_makecohortvcf
    File final_sample_list = GATKSVPipelinePhase1.batch_samples_postOutlierExclusion_file
    File final_sample_outlier_list = GATKSVPipelinePhase1.outlier_samples_excluded_file

    # Additional outputs for creating a reference panel
    Array[File] counts = counts_files_
    Array[File] PE_files = pe_files_
    Array[File] PE_files_index = pe_files_index_
    Array[File] SR_files = sr_files_
    Array[File] SR_files_index = sr_files_index_
    Array[File]? manta_vcfs = manta_vcfs_
    Array[File]? manta_vcfs_index = manta_vcfs_index_
    Array[File]? melt_vcfs = melt_vcfs_
    Array[File]? melt_vcfs_index = melt_vcfs_index_
    Array[File]? wham_vcfs = wham_vcfs_
    Array[File]? wham_vcfs_index = wham_vcfs_index_

    File medianfile = GATKSVPipelinePhase1.median_cov
    File merged_coverage_file = GATKSVPipelinePhase1.merged_bincov
    File merged_coverage_file_index = GATKSVPipelinePhase1.merged_bincov_index
    File merged_baf_file = GATKSVPipelinePhase1.merged_BAF
    File merged_baf_file_index = GATKSVPipelinePhase1.merged_BAF_index
    File merged_disc_file = GATKSVPipelinePhase1.merged_PE
    File merged_disc_file_index = GATKSVPipelinePhase1.merged_PE_index
    File merged_split_file = GATKSVPipelinePhase1.merged_SR
    File merged_split_file_index = GATKSVPipelinePhase1.merged_SR_index

    File del_bed = GATKSVPipelinePhase1.merged_dels
    File del_bed_index = GATKSVPipelinePhase1.merged_dels + ".tbi"
    File dup_bed = GATKSVPipelinePhase1.merged_dups
    File dup_bed_index = GATKSVPipelinePhase1.merged_dups + ".tbi"

    File? std_manta_vcf_tar = GATKSVPipelinePhase1.std_manta_vcf_tar
    File? std_melt_vcf_tar = GATKSVPipelinePhase1.std_melt_vcf_tar
    File? std_scramble_vcf_tar = GATKSVPipelinePhase1.std_scramble_vcf_tar
    File? std_wham_vcf_tar = GATKSVPipelinePhase1.std_wham_vcf_tar

    File merged_depth_vcf = GATKSVPipelinePhase1.depth_vcf
    File merged_depth_vcf_index = GATKSVPipelinePhase1.depth_vcf_index
    File? merged_manta_vcf = GATKSVPipelinePhase1.manta_vcf
    File? merged_manta_vcf_index = GATKSVPipelinePhase1.manta_vcf_index
    File? merged_melt_vcf = GATKSVPipelinePhase1.melt_vcf
    File? merged_melt_vcf_index = GATKSVPipelinePhase1.melt_vcf_index
    File? merged_wham_vcf = GATKSVPipelinePhase1.wham_vcf
    File? merged_wham_vcf_index = GATKSVPipelinePhase1.wham_vcf_index
    Array[File] ?clustered_sv_counts = GATKSVPipelinePhase1.clustered_sv_counts
    Array[File]? clustered_sv_count_plots = GATKSVPipelinePhase1.clustered_sv_count_plots
    File? clustered_outlier_samples_preview = GATKSVPipelinePhase1.clustered_outlier_samples_preview
    File? clustered_outlier_samples_with_reason = GATKSVPipelinePhase1.clustered_outlier_samples_with_reason
    Int? clustered_num_outlier_samples = GATKSVPipelinePhase1.clustered_num_outlier_samples

    File evidence_metrics = GATKSVPipelinePhase1.evidence_metrics
    File evidence_metrics_common = GATKSVPipelinePhase1.evidence_metrics_common

    File filtered_depth_vcf = select_first([GATKSVPipelinePhase1.filtered_depth_vcf])
    File filtered_pesr_vcf = select_first([GATKSVPipelinePhase1.filtered_pesr_vcf])
    File cohort_pesr_vcf = select_first([GATKSVPipelinePhase1.filtered_pesr_vcf])
    File cohort_depth_vcf = select_first([GATKSVPipelinePhase1.filtered_depth_vcf])
    File? sites_filtered_manta_vcf = GATKSVPipelinePhase1.sites_filtered_manta_vcf
    File? sites_filtered_wham_vcf = GATKSVPipelinePhase1.sites_filtered_wham_vcf
    File? sites_filtered_melt_vcf = GATKSVPipelinePhase1.sites_filtered_melt_vcf
    File? sites_filtered_depth_vcf = GATKSVPipelinePhase1.sites_filtered_depth_vcf
    File cutoffs = GATKSVPipelinePhase1.cutoffs
    File genotyped_pesr_vcf = GenotypeBatch.genotyped_pesr_vcf
    File genotyped_depth_vcf = GenotypeBatch.genotyped_depth_vcf
    File regeno_coverage_medians = GenotypeBatch.regeno_coverage_medians
    File regenotyped_depth_vcf = RegenotypeCNVs.regenotyped_depth_vcfs[0]

    File genotype_pesr_pesr_sepcutoff = select_first([GenotypeBatch.trained_genotype_pesr_pesr_sepcutoff])
    File genotype_pesr_depth_sepcutoff = select_first([GenotypeBatch.trained_genotype_pesr_depth_sepcutoff])
    File genotype_depth_pesr_sepcutoff = select_first([GenotypeBatch.trained_genotype_depth_pesr_sepcutoff])
    File genotype_depth_depth_sepcutoff = select_first([GenotypeBatch.trained_genotype_depth_depth_sepcutoff])
    File depth_gt_rd_sep_file = select_first([GenotypeBatch.trained_genotype_depth_depth_sepcutoff])
    File PE_metrics = select_first([GenotypeBatch.trained_PE_metrics])
    File SR_metrics = select_first([GenotypeBatch.trained_SR_metrics])
    File raw_sr_bothside_pass_file = GenotypeBatch.sr_bothside_pass
    File raw_sr_background_fail_file = GenotypeBatch.sr_background_fail

    # CombineBatches
    Array[File] combined_vcfs = MakeCohortVcf.combined_vcfs
    Array[File] combined_vcf_indexes = MakeCohortVcf.combined_vcf_indexes
    Array[File] cluster_bothside_pass_lists = MakeCohortVcf.cluster_bothside_pass_lists
    Array[File] cluster_background_fail_lists = MakeCohortVcf.cluster_background_fail_lists

    # ResolveComplexVariants
    Array[File] complex_resolve_vcfs = MakeCohortVcf.complex_resolve_vcfs
    Array[File] complex_resolve_vcf_indexes = MakeCohortVcf.complex_resolve_vcf_indexes
    Array[File] complex_resolve_bothside_pass_lists = MakeCohortVcf.complex_resolve_bothside_pass_lists
    Array[File] complex_resolve_background_fail_lists = MakeCohortVcf.complex_resolve_background_fail_lists
    Array[File] breakpoint_overlap_dropped_record_vcfs = MakeCohortVcf.breakpoint_overlap_dropped_record_vcfs
    Array[File] breakpoint_overlap_dropped_record_vcf_indexes = MakeCohortVcf.breakpoint_overlap_dropped_record_vcf_indexes

    # GenotypeComplexVariants
    Array[File] complex_genotype_vcfs = MakeCohortVcf.complex_genotype_vcfs
    Array[File] complex_genotype_vcf_indexes = MakeCohortVcf.complex_genotype_vcfs
  }
}

