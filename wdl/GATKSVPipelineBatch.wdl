version 1.0

import "Module00aBatch.wdl" as m00a
import "Module00b.wdl" as m00b
import "GATKSVPipelinePhase1.wdl" as phase1
import "Module04.wdl" as m04
import "Module04b.wdl" as m04b
import "Module0506.wdl" as m0506
import "Utils.wdl" as utils
import "Structs.wdl"
import "TestUtils.wdl" as tu

# GATK SV Pipeline batch mode
# Runs modules 00abc, 01, 02, 03, 04, 0506

workflow GATKSVPipelineBatch {
  input {
    # Batch data
    String batch
    Array[String] sample_ids

    # Required unless caller and evidence outputs are provided (below)
    Array[File]? bam_or_cram_files
    Array[File]? bam_or_cram_indexes

    # Optionally provide calls and evidence (override caller flags below)
    Array[File]? counts_files
    Array[File]? pe_files
    Array[File]? sr_files
    Array[File?]? baf_files
    Array[File]? delly_vcfs
    Array[File]? manta_vcfs
    Array[File]? melt_vcfs
    Array[File]? wham_vcfs

    # Enable different callers
    Boolean use_delly = false
    Boolean use_manta = true
    Boolean use_melt = true
    Boolean use_wham = true

    # BAF Generation (if baf_files unavailable)
    # BAF Option #1 (provide all)
    # From single-sample gVCFS
    Array[File]? gvcfs

    # BAF Option #2
    # From multi-sample VCFs (sharded by position)
    Array[File]? snp_vcfs
    File? snp_vcf_header # Required only if VCFs are unheadered

    # Merge contig vcfs at each stage of Module 0506 for QC
    Boolean module0506_merge_cluster_vcfs = false
    Boolean module0506_merge_complex_resolve_vcfs = false
    Boolean module0506_merge_complex_genotype_vcfs = false

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

    # Run module metrics - all modules on by default for batch WDL
    Boolean? run_00a_metrics
    Boolean? run_00c_metrics = true  # 00c metrics is off by default standalone but on for batch WDL
    Boolean? run_01_metrics
    Boolean? run_02_metrics
    Boolean? run_03_metrics
    Boolean? run_04_metrics
    Boolean? run_0506_metrics

    File? baseline_00a_metrics
    File? baseline_00c_metrics
    File? baseline_01_metrics
    File? baseline_02_metrics
    File? baseline_03_metrics
    File? baseline_04_metrics
    File? baseline_0506_metrics

    String sv_base_mini_docker
    String sv_base_docker
    String sv_pipeline_docker
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
    String? delly_docker
    String? manta_docker
    String? melt_docker
    String? wham_docker
    String cloud_sdk_docker

    # Do not use
    Array[File]? NONE_ARRAY_
    String? NONE_STRING_
  }

  Boolean collect_coverage_ = !defined(counts_files)
  Boolean collect_pesr_ = !(defined(pe_files) && defined(sr_files))

  String? delly_docker_ = if (!defined(delly_vcfs) && use_delly) then delly_docker else NONE_STRING_
  String? manta_docker_ = if (!defined(manta_vcfs) && use_manta) then manta_docker else NONE_STRING_
  String? melt_docker_ = if (!defined(melt_vcfs) && use_melt) then melt_docker else NONE_STRING_
  String? wham_docker_ = if (!defined(wham_vcfs) && use_wham) then wham_docker else NONE_STRING_

  Boolean run_module00a = collect_coverage_ || collect_pesr_ || defined(delly_docker_) || defined(manta_docker_) || defined(melt_docker_) || defined(wham_docker_)

  if (run_module00a) {
    call m00a.Module00aBatch {
      input:
        bam_or_cram_files=select_first([bam_or_cram_files]),
        bam_or_cram_indexes=bam_or_cram_indexes,
        collect_coverage=collect_coverage_,
        collect_pesr=collect_pesr_,
        sample_ids=sample_ids,
        primary_contigs_list=primary_contigs_list,
        reference_fasta=reference_fasta,
        reference_index=reference_index,
        reference_dict=reference_dict,
        run_module_metrics = run_00a_metrics,
        primary_contigs_fai = primary_contigs_fai,
        batch = batch,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        linux_docker = linux_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        delly_docker=delly_docker,
        manta_docker=manta_docker,
        melt_docker=melt_docker,
        wham_docker=wham_docker,
        gatk_docker=gatk_docker,
        gatk_docker_pesr_override = gatk_docker_pesr_override,
        genomes_in_the_cloud_docker=genomes_in_the_cloud_docker,
        samtools_cloud_docker=samtools_cloud_docker,
        cloud_sdk_docker = cloud_sdk_docker
    }
  }

  Array[File] counts_files_ = if collect_coverage_ then select_all(select_first([Module00aBatch.coverage_counts])) else select_first([counts_files])
  Array[File] pe_files_ = if collect_pesr_ then select_all(select_first([Module00aBatch.pesr_disc])) else select_first([pe_files])
  Array[File] sr_files_ = if collect_pesr_ then select_all(select_first([Module00aBatch.pesr_split])) else select_first([sr_files])

  if (use_delly) {
    Array[File] delly_vcfs_ = if defined(delly_vcfs) then select_first([delly_vcfs]) else select_all(select_first([Module00aBatch.delly_vcf]))
  }
  if (use_manta) {
    Array[File] manta_vcfs_ = if defined(manta_vcfs) then select_first([manta_vcfs]) else select_all(select_first([Module00aBatch.manta_vcf]))
  }
  if (use_melt) {
    Array[File] melt_vcfs_ = if defined(melt_vcfs) then select_first([melt_vcfs]) else select_all(select_first([Module00aBatch.melt_vcf]))
  }
  if (use_wham) {
    Array[File] wham_vcfs_ = if defined(wham_vcfs) then select_first([wham_vcfs]) else select_all(select_first([Module00aBatch.wham_vcf]))
  }

  call m00b.Module00b as Module00b {
    input:
      batch=batch,
      samples=sample_ids,
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
      batch=batch,
      samples=sample_ids,
      ped_file=ped_file,
      genome_file=genome_file,
      contigs=primary_contigs_fai,
      reference_fasta=reference_fasta,
      reference_index=reference_index,
      reference_dict=reference_dict,
      BAF_files=baf_files,
      counts=counts_files_,
      bincov_matrix=Module00b.bincov_matrix,
      bincov_matrix_index=Module00b.bincov_matrix_index,
      PE_files=pe_files_,
      SR_files=sr_files_,
      delly_vcfs=delly_vcfs_,
      manta_vcfs=manta_vcfs_,
      melt_vcfs=melt_vcfs_,
      wham_vcfs=wham_vcfs_,
      gvcfs=gvcfs,
      snp_vcfs=snp_vcfs,
      snp_vcf_header=snp_vcf_header,
      cnmops_chrom_file=autosome_file,
      cnmops_allo_file=allosome_file,
      allosome_contigs=allosome_file,
      autosome_contigs=autosome_file,
      run_00c_metrics = run_00c_metrics,
      run_01_metrics = run_01_metrics,
      run_02_metrics = run_02_metrics,
      run_03_metrics = run_03_metrics,
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

  call m04.Module04 as Module04 {
    input:
      batch_pesr_vcf=GATKSVPipelinePhase1.filtered_pesr_vcf,
      batch_depth_vcf=select_first([GATKSVPipelinePhase1.filtered_depth_vcf]),
      cohort_pesr_vcf=GATKSVPipelinePhase1.filtered_pesr_vcf,
      cohort_depth_vcf=select_first([GATKSVPipelinePhase1.filtered_depth_vcf]),
      batch=batch,
      rf_cutoffs=GATKSVPipelinePhase1.cutoffs,
      medianfile=GATKSVPipelinePhase1.median_cov,
      coveragefile=GATKSVPipelinePhase1.merged_bincov,
      coveragefile_index=GATKSVPipelinePhase1.merged_bincov_index,
      discfile=GATKSVPipelinePhase1.merged_PE,
      discfile_index=GATKSVPipelinePhase1.merged_PE_index,
      splitfile=GATKSVPipelinePhase1.merged_SR,
      splitfile_index=GATKSVPipelinePhase1.merged_SR_index,
      ped_file=ped_file,
      ref_dict=reference_dict,
      run_module_metrics = run_04_metrics,
      primary_contigs_list = primary_contigs_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
      linux_docker=linux_docker
  }

  call m04b.Module04b as Module04b {
    input:
      depth_vcfs=[Module04.genotyped_depth_vcf],
      batch_depth_vcfs=[select_first([GATKSVPipelinePhase1.filtered_depth_vcf])],
      cohort_depth_vcf=select_first([GATKSVPipelinePhase1.filtered_depth_vcf]),
      batches=[batch],
      cohort=batch,
      medianfiles=[GATKSVPipelinePhase1.median_cov],
      coveragefiles=[GATKSVPipelinePhase1.merged_bincov],
      coveragefile_idxs=[GATKSVPipelinePhase1.merged_bincov_index],
      ped_file=ped_file,
      RD_depth_sepcutoffs=[select_first([Module04.trained_genotype_depth_depth_sepcutoff])],
      contig_list=primary_contigs_list,
      regeno_coverage_medians=[Module04.regeno_coverage_medians],
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
      sv_pipeline_base_docker=sv_pipeline_base_docker
  }
  

  call m0506.Module0506 as Module0506 {
    input:
      merge_cluster_vcfs = module0506_merge_cluster_vcfs,
      merge_complex_resolve_vcfs = module0506_merge_complex_resolve_vcfs,
      merge_complex_genotype_vcfs = module0506_merge_complex_genotype_vcfs,
      raw_sr_bothside_pass_files=[Module04.sr_bothside_pass],
      raw_sr_background_fail_files=[Module04.sr_background_fail],
      ped_file=ped_file,
      pesr_vcfs=[Module04.genotyped_pesr_vcf],
      depth_vcfs=Module04b.regenotyped_depth_vcfs,
      contig_list=primary_contigs_fai,
      allosome_fai=allosome_file,
      ref_dict=reference_dict,
      disc_files=[GATKSVPipelinePhase1.merged_PE],
      disc_files_index=[GATKSVPipelinePhase1.merged_PE_index],
      bincov_files=[GATKSVPipelinePhase1.merged_bincov],
      cohort_name=batch,
      rf_cutoff_files=[GATKSVPipelinePhase1.cutoffs],
      batches=[batch],
      depth_gt_rd_sep_files=[select_first([Module04.trained_genotype_depth_depth_sepcutoff])],
      median_coverage_files=[GATKSVPipelinePhase1.median_cov],
      run_module_metrics = run_0506_metrics,
      primary_contigs_list = primary_contigs_list,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      linux_docker=linux_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
      sv_pipeline_qc_docker=sv_pipeline_qc_docker,
      sv_base_mini_docker=sv_base_mini_docker
  }

  call tu.CatMetrics as CatBatchMetrics {
      input:
        prefix = "batch_sv." + batch,
        metric_files = select_all([Module00aBatch.metrics_file_00a, GATKSVPipelinePhase1.metrics_file_00c, GATKSVPipelinePhase1.metrics_file_01, GATKSVPipelinePhase1.metrics_file_02, GATKSVPipelinePhase1.metrics_file_03, Module04.metrics_file_04, Module0506.metrics_file_0506]),
        linux_docker = linux_docker
    }

  Array[File] defined_baseline_metrics = select_all([baseline_00a_metrics, baseline_00c_metrics, baseline_01_metrics, baseline_02_metrics, baseline_03_metrics, baseline_04_metrics, baseline_0506_metrics])
  if (length(defined_baseline_metrics) > 0) {
    call tu.CatMetrics as CatBaselineMetrics {
      input:
        prefix = "baseline." + batch,
        metric_files = defined_baseline_metrics,
        linux_docker = linux_docker
    }
    call tu.PlotMetrics {
      input:
        name = batch,
        samples = sample_ids,
        test_metrics = CatBatchMetrics.out,
        base_metrics = CatBaselineMetrics.out,
        sv_pipeline_base_docker = sv_pipeline_base_docker
    }
  }

  call utils.RunQC as BatchQC {
    input:
      name = batch,
      metrics = CatBatchMetrics.out,
      sv_pipeline_base_docker = sv_pipeline_base_docker
  }

  output {
    File vcf = Module0506.vcf
    File vcf_index = Module0506.vcf_index
    File metrics_file_batch = CatBatchMetrics.out
    File qc_file = BatchQC.out

    # Additional outputs for creating a reference panel
    Array[File] pesr_disc_files = pe_files_
    Array[File] pesr_split_files = sr_files_
    Array[File]? std_delly_vcfs = GATKSVPipelinePhase1.std_delly_vcf
    Array[File]? std_manta_vcfs = GATKSVPipelinePhase1.std_manta_vcf
    Array[File]? std_melt_vcfs = GATKSVPipelinePhase1.std_melt_vcf
    Array[File]? std_wham_vcfs = GATKSVPipelinePhase1.std_wham_vcf
    File bincov_matrix = GATKSVPipelinePhase1.merged_bincov
    File del_bed = GATKSVPipelinePhase1.merged_dels
    File dup_bed = GATKSVPipelinePhase1.merged_dups

    File final_sample_list = GATKSVPipelinePhase1.batch_samples_postOutlierExclusion_file
    File final_sample_outlier_list = GATKSVPipelinePhase1.outlier_samples_excluded_file

    File cutoffs = GATKSVPipelinePhase1.cutoffs
    File genotype_pesr_pesr_sepcutoff = select_first([Module04.trained_genotype_pesr_pesr_sepcutoff])
    File genotype_pesr_depth_sepcutoff = select_first([Module04.trained_genotype_pesr_depth_sepcutoff])
    File genotype_depth_pesr_sepcutoff = select_first([Module04.trained_genotype_depth_pesr_sepcutoff])
    File genotype_depth_depth_sepcutoff = select_first([Module04.trained_genotype_depth_depth_sepcutoff])
    File PE_metrics = select_first([Module04.trained_PE_metrics])
    File SR_metrics = select_first([Module04.trained_SR_metrics])
  }
}

