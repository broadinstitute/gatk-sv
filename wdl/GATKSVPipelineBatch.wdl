version 1.0

import "GatherSampleEvidenceBatch.wdl" as sampleevidence
import "EvidenceQC.wdl" as evidenceqc
import "GATKSVPipelinePhase1.wdl" as phase1
import "GenotypeBatch.wdl" as genotypebatch
import "RegenotypeCNVs.wdl" as regenocnvs
import "MakeCohortVcf.wdl" as makecohortvcf
import "TasksClusterBatch.wdl" as tasks_cluster
import "TasksMakeCohortVcf.wdl" as tasks_makecohortvcf
import "StripyWorkflow.wdl" as stripy
import "AnnotateVcf.wdl" as annotate
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
    Array[File]? stripy_vcfs_input

    # Enable different callers
    Boolean use_manta = true
    Boolean use_melt = false
    Boolean use_scramble = true
    Boolean use_wham = true
    Boolean use_stripy = false

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

    # Annotation resources
    File? protein_coding_gtf
    File? noncoding_bed
    Int? promoter_window
    Int? max_breakend_as_cnv_length
    String? svannotate_additional_args
    File? sample_pop_assignments
    File? sample_keep_list
    File? par_bed
    File? allosomes_list
    Int annotation_sv_per_shard = 5000
    File? external_af_ref_bed
    String? external_af_ref_bed_prefix
    Array[String]? external_af_population

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
    String sv_pipeline_qc_docker
    String linux_docker
    String cnmops_docker
    String gatk_docker
    String? gcnv_gatk_docker
    String genomes_in_the_cloud_docker
    String samtools_cloud_docker
    String? manta_docker
    String? melt_docker
    String? scramble_docker
    String? wham_docker
    String? stripy_docker
    String cloud_sdk_docker

    # Batch metrics
    RuntimeAttr? runtime_attr_cat_metrics
    RuntimeAttr? runtime_attr_plot_metrics

    # Batch ploidy generation
    RuntimeAttr? runtime_attr_create_ploidy

    # AnnotateVcf
    RuntimeAttr? runtime_attr_svannotate
    RuntimeAttr? runtime_attr_scatter_vcf
    RuntimeAttr? runtime_attr_subset_vcf_by_samples_list
    RuntimeAttr? runtime_attr_compute_AFs
    RuntimeAttr? runtime_attr_modify_vcf
    RuntimeAttr? runtime_attr_split_ref_bed
    RuntimeAttr? runtime_attr_split_query_vcf
    RuntimeAttr? runtime_attr_bedtools_closest
    RuntimeAttr? runtime_attr_select_matched_svs
    RuntimeAttr? runtime_attr_concat
    RuntimeAttr? runtime_attr_preconcat
    RuntimeAttr? runtime_attr_fix_header
    RuntimeAttr? runtime_attr_merge_stripy_vcf

    # Do not use
    Array[File]? NONE_ARRAY_
    File? NONE_FILE_
    String? NONE_STRING_
  }

  Boolean collect_coverage_ = !defined(counts_files_input)
  Boolean collect_pesr_ = !(defined(pe_files_input) && defined(sr_files_input) && defined(sd_files_input))

  String? manta_docker_ = if (!defined(manta_vcfs_input) && use_manta) then manta_docker else NONE_STRING_
  String? melt_docker_ = if (!defined(melt_vcfs_input) && use_melt) then melt_docker else NONE_STRING_
  String? scramble_docker_ = if (!defined(scramble_vcfs_input) && use_scramble) then scramble_docker else NONE_STRING_
  String? wham_docker_ = if (!defined(wham_vcfs_input) && use_wham) then wham_docker else NONE_STRING_
  Boolean run_stripy_ = use_stripy && !defined(stripy_vcfs_input)

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
        linux_docker = linux_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        manta_docker=manta_docker_,
        melt_docker=melt_docker_,
        scramble_docker=scramble_docker_,
        wham_docker=wham_docker_,
        gatk_docker=gatk_docker,
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

  if (run_stripy_) {
    scatter (i in range(length(samples))) {
      call stripy.StripyWorkflow {
        input:
          bam_or_cram_file = select_first([bam_or_cram_files])[i],
          bam_or_cram_index = if defined(bam_or_cram_indexes) then select_first([bam_or_cram_indexes])[i] else NONE_FILE_,
          sample_name = samples[i],
          ped_file = ped_file,
          reference_fasta = reference_fasta,
          linux_docker = linux_docker,
          stripy_docker = select_first([stripy_docker])
      }
    }
  }

  Array[File] generated_stripy_vcfs_ = select_all(select_first([StripyWorkflow.stripy_vcf, []]))
  Array[File]? stripy_vcfs_ = if use_stripy then (if defined(stripy_vcfs_input) then select_first([stripy_vcfs_input]) else generated_stripy_vcfs_) else NONE_ARRAY_

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
      stripy_vcfs=stripy_vcfs_,

      cnmops_chrom_file=autosome_file,
      cnmops_allo_file=allosome_file,
      run_batchevidence_metrics = run_batchevidence_metrics,
      run_clusterbatch_metrics = run_clusterbatch_metrics,
      run_batchmetrics_metrics = run_batchmetrics_metrics,
      run_filterbatch_metrics = run_filterbatch_metrics,
      primary_contigs_list = primary_contigs_list,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_base_docker=sv_base_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_qc_docker=sv_pipeline_qc_docker,
      linux_docker=linux_docker,
      cnmops_docker=cnmops_docker,
      gatk_docker=gatk_docker,
      gcnv_gatk_docker=gcnv_gatk_docker,
      runtime_attr_merge_stripy_vcf_cluster_batch=runtime_attr_merge_stripy_vcf
  }

  Array[File] stripy_vcfs_for_annotation_ = select_all([GATKSVPipelinePhase1.merged_stripy_vcf])
  Array[File] merge_vcfs_ = select_all([GATKSVPipelinePhase1.filtered_pesr_vcf, GATKSVPipelinePhase1.filtered_depth_vcf])
  call tasks_makecohortvcf.ConcatVcfs as MergePesrDepthVcfs {
    input:
    vcfs = merge_vcfs_,
    vcfs_idx = [merge_vcfs_[0] + ".tbi", merge_vcfs_[1] + ".tbi"],
    allow_overlaps = true,
    outfile_prefix = "~{name}.merge_pesr_depth",
    sv_base_mini_docker = sv_base_mini_docker
  }

  call tasks_cluster.CreatePloidyTableFromPed {
    input:
      ped_file = ped_file,
      contig_list = primary_contigs_list,
      retain_female_chr_y = false,
      chr_x = chr_x,
      chr_y = chr_y,
      output_prefix = "~{name}.ploidy",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_create_ploidy
  }

  call genotypebatch.GenotypeBatch as GenotypeBatch {
    input:
      vcf=MergePesrDepthVcfs.concat_vcf,
      batch=name,
      rf_cutoffs=GATKSVPipelinePhase1.cutoffs,
      median_coverage=GATKSVPipelinePhase1.median_cov,
      rd_file=GATKSVPipelinePhase1.merged_bincov,
      pe_file=GATKSVPipelinePhase1.merged_PE,
      sr_file=GATKSVPipelinePhase1.merged_SR,
      reference_dict=reference_dict,
        ploidy_table=CreatePloidyTableFromPed.out,
      contig_list = primary_contigs_list,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      gatk_docker=gatk_docker
  }

  call regenocnvs.RegenotypeCNVs as RegenotypeCNVs {
    input:
      depth_vcfs=[GenotypeBatch.genotyped_depth_vcf],
      batch_depth_vcfs=[select_first([GATKSVPipelinePhase1.filtered_depth_vcf])],
      batches=[name],
      cohort=name,
      medianfiles=[GATKSVPipelinePhase1.median_cov],
      coveragefiles=[GATKSVPipelinePhase1.merged_bincov],
      coveragefile_idxs=[GATKSVPipelinePhase1.merged_bincov_index],
      genotyping_rd_table=[select_first([GenotypeBatch.genotyping_rd_table])],
        ploidy_tables=[CreatePloidyTableFromPed.out],
      contig_list=primary_contigs_list,
      regeno_coverage_medians=[GenotypeBatch.regeno_coverage_medians],
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker
  }
  

  call makecohortvcf.MakeCohortVcf as MakeCohortVcf {
    input:
      merge_cluster_vcfs = makecohortvcf_merge_cluster_vcfs,
      merge_complex_resolve_vcfs = makecohortvcf_merge_complex_resolve_vcfs,
      merge_complex_genotype_vcfs = makecohortvcf_merge_complex_genotype_vcfs,
      ped_file=ped_file,
      pesr_vcfs=[GenotypeBatch.genotyped_pesr_vcf],
      depth_vcfs=RegenotypeCNVs.regenotyped_depth_vcfs,
      contig_list=primary_contigs_fai,
      allosome_fai=allosome_file,
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_index,
      reference_dict=reference_dict,
      chr_x=chr_x,
      chr_y=chr_y,
      disc_files=[GATKSVPipelinePhase1.merged_PE],
      bincov_files=[GATKSVPipelinePhase1.merged_bincov],
      cohort_name=name,
      rf_cutoff_files=[GATKSVPipelinePhase1.cutoffs],
      batches=[name],
      genotyping_rd_tables=[select_first([GenotypeBatch.genotyping_rd_table])],
      median_coverage_files=[GATKSVPipelinePhase1.median_cov],
      run_module_metrics = run_makecohortvcf_metrics,
      primary_contigs_list = primary_contigs_list,
      linux_docker=linux_docker,
      gatk_docker=gatk_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_qc_docker=sv_pipeline_qc_docker,
      sv_base_mini_docker=sv_base_mini_docker
  }

  call annotate.AnnotateVcf {
    input:
      vcf = MakeCohortVcf.vcf,
      contig_list = primary_contigs_list,
      prefix = name,
      stripy_vcfs = stripy_vcfs_for_annotation_,
      protein_coding_gtf = protein_coding_gtf,
      noncoding_bed = noncoding_bed,
      promoter_window = promoter_window,
      max_breakend_as_cnv_length = max_breakend_as_cnv_length,
      svannotate_additional_args = svannotate_additional_args,
      sample_pop_assignments = sample_pop_assignments,
      sample_keep_list = sample_keep_list,
      ped_file = ped_file,
      par_bed = par_bed,
      allosomes_list = allosomes_list,
      sv_per_shard = annotation_sv_per_shard,
      external_af_ref_bed = external_af_ref_bed,
      external_af_ref_prefix = external_af_ref_bed_prefix,
      external_af_population = external_af_population,
      sv_pipeline_docker = sv_pipeline_docker,
      sv_base_mini_docker = sv_base_mini_docker,
      gatk_docker = gatk_docker,
      runtime_attr_svannotate = runtime_attr_svannotate,
      runtime_attr_scatter_vcf = runtime_attr_scatter_vcf,
      runtime_attr_subset_vcf_by_samples_list = runtime_attr_subset_vcf_by_samples_list,
      runtime_attr_compute_AFs = runtime_attr_compute_AFs,
      runtime_attr_modify_vcf = runtime_attr_modify_vcf,
      runtime_attr_split_ref_bed = runtime_attr_split_ref_bed,
      runtime_attr_split_query_vcf = runtime_attr_split_query_vcf,
      runtime_attr_bedtools_closest = runtime_attr_bedtools_closest,
      runtime_attr_select_matched_svs = runtime_attr_select_matched_svs,
      runtime_attr_concat = runtime_attr_concat,
      runtime_attr_preconcat = runtime_attr_preconcat,
      runtime_attr_fix_header = runtime_attr_fix_header,
      runtime_attr_merge_stripy_vcf = runtime_attr_merge_stripy_vcf
  }

  call tu.CatMetrics as CatBatchMetrics {
      input:
        prefix = "batch_sv." + name,
        metric_files = select_all([GatherSampleEvidenceBatch.metrics_file_sampleevidence, GATKSVPipelinePhase1.metrics_file_batchevidence, GATKSVPipelinePhase1.metrics_file_clusterbatch, GATKSVPipelinePhase1.metrics_file_batchmetrics, GATKSVPipelinePhase1.metrics_file_filterbatch, MakeCohortVcf.metrics_file_makecohortvcf]),
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
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_plot_metrics
    }
  }

  call utils.RunQC as BatchQC {
    input:
      name = name,
      metrics = CatBatchMetrics.out,
      qc_definitions = qc_definitions,
      sv_pipeline_docker = sv_pipeline_docker
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
    File annotated_vcf = AnnotateVcf.annotated_vcf
    File annotated_vcf_index = AnnotateVcf.annotated_vcf_index
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
    Array[File]? stripy_vcfs = stripy_vcfs_
    Array[File]? stripy_json_outputs = StripyWorkflow.stripy_json
    Array[File]? stripy_tsv_outputs = StripyWorkflow.stripy_tsv
    Array[File]? stripy_html_outputs = StripyWorkflow.stripy_html
    File? merged_stripy_vcf = GATKSVPipelinePhase1.merged_stripy_vcf
    File? merged_stripy_vcf_index = GATKSVPipelinePhase1.merged_stripy_vcf_index

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

    File genotyping_rd_table = GenotypeBatch.genotyping_rd_table
    File genotyping_pe_table = GenotypeBatch.genotyping_pe_table
    File genotyping_sr_table = GenotypeBatch.genotyping_sr_table

    # CombineBatches
    Array[File] combined_vcfs = MakeCohortVcf.combined_vcfs
    Array[File] combined_vcf_indexes = MakeCohortVcf.combined_vcf_indexes
    Array[File] cluster_bothside_pass_lists = MakeCohortVcf.cluster_bothside_pass_lists
    Array[File] cluster_background_fail_lists = MakeCohortVcf.cluster_background_fail_lists

    # ResolveComplexVariants
    Array[File] complex_resolve_vcfs = MakeCohortVcf.complex_resolve_vcfs
    Array[File] complex_resolve_vcf_indexes = MakeCohortVcf.complex_resolve_vcf_indexes
    File complex_resolve_bothside_pass_list = MakeCohortVcf.complex_resolve_bothside_pass_list
    File complex_resolve_background_fail_list = MakeCohortVcf.complex_resolve_background_fail_list
    Array[File] breakpoint_overlap_dropped_record_vcfs = MakeCohortVcf.breakpoint_overlap_dropped_record_vcfs
    Array[File] breakpoint_overlap_dropped_record_vcf_indexes = MakeCohortVcf.breakpoint_overlap_dropped_record_vcf_indexes

    # GenotypeComplexVariants
    Array[File] complex_genotype_vcfs = MakeCohortVcf.complex_genotype_vcfs
    Array[File] complex_genotype_vcf_indexes = MakeCohortVcf.complex_genotype_vcfs
  }
}

