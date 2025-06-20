version 1.0

import "AnnotateVcf.wdl" as annotate
import "FilterGenotypes.wdl" as fg
import "JoinRawCalls.wdl" as jrc
import "RefineComplexVariants.wdl" as rcv
import "SingleSampleFiltering.wdl" as SingleSampleFiltering
import "GATKSVPipelineSingleSampleMetrics.wdl" as SingleSampleMetrics
import "SVConcordance.wdl" as svc
import "Utils.wdl" as utils

import "Structs.wdl"

workflow GATKSVPipelineSingleSamplePart2 {
  input {
    File cleaned_vcf
    String sample_id
    String batch
    File ped_file  # Ped file with sample and panel

    # RefineComplexVariants
    File samples_list
    File PE_metrics
    File PE_metrics_index
    File Depth_DEL_bed
    File Depth_DUP_bed

    Int? min_pe_cpx
    Int? min_pe_ctx

    # JoinRawCalls
    File? clustered_manta_vcf
    File? clustered_manta_vcf_index
    File? clustered_melt_vcf
    File? clustered_melt_vcf_index
    File? clustered_scramble_vcf
    File? clustered_scramble_vcf_index
    File? clustered_wham_vcf
    File? clustered_wham_vcf_index
    File? clustered_depth_vcf
    File? clustered_depth_vcf_index

    # FilterGenotypes
    File gq_recalibrator_model_file

    Array[String] recalibrate_gq_args = []
    Array[File] genome_tracks = []
    Float no_call_rate_cutoff = 0.05  # Set to 1 to disable NCR filtering
    String? sl_filter_args  # Explicitly set SL cutoffs. See apply_sl_filter.py for arguments.

    Boolean run_main_vcf_qc = false

    # Downstream steps
    File wgd_scores
    File case_counts_file
    File case_pe_file
    File case_sr_file

    File convert_cnvs_without_depth_support_vcf
    File genotyped_depth_vcf
    File non_genotyped_unique_depth_calls_vcf

    File ref_samples_list
    File qc_definitions

    File protein_coding_gtf
    File noncoding_bed
    Int? promoter_window
    Int? max_breakend_as_cnv_length
    Int annotation_sv_per_shard

    File? external_af_ref_bed             # bed file with population AFs for annotation
    String? external_af_ref_bed_prefix    # name of external AF bed file call set
    Array[String]? external_af_population # populations to annotate external AFs (required if ref_bed set, use "ALL" for all)

    # Reference genome
    String? chr_x
    String? chr_y
    File contig_list
    File primary_contigs_fai
    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    # Dockers
    String gatk_docker
    String gq_recalibrator_gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    String linux_docker

    # Do not use
    Array[File]? NONE_ARRAY_
  }

  call rcv.RefineComplexVariants {
    input:
      vcf=cleaned_vcf,
      prefix=sample_id,
      batch_name_list=[sample_id],
      batch_sample_lists=[samples_list],
      PE_metrics=[PE_metrics],
      PE_metrics_indexes=[PE_metrics_index],
      Depth_DEL_beds=[Depth_DEL_bed],
      Depth_DUP_beds=[Depth_DUP_bed],
      min_pe_cpx=min_pe_cpx,
      min_pe_ctx=min_pe_ctx,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      linux_docker=linux_docker
  }

  call jrc.JoinRawCalls {
    input:
      prefix=sample_id,
      clustered_manta_vcfs=if defined(clustered_manta_vcf) then select_all([clustered_manta_vcf]) else NONE_ARRAY_,
      clustered_manta_vcf_indexes=if defined(clustered_manta_vcf_index) then select_all([clustered_manta_vcf_index]) else NONE_ARRAY_,
      clustered_melt_vcfs=if defined(clustered_melt_vcf) then select_all([clustered_melt_vcf]) else NONE_ARRAY_,
      clustered_melt_vcf_indexes=if defined(clustered_melt_vcf_index) then select_all([clustered_melt_vcf_index]) else NONE_ARRAY_,
      clustered_scramble_vcfs=if defined(clustered_scramble_vcf) then select_all([clustered_scramble_vcf]) else NONE_ARRAY_,
      clustered_scramble_vcf_indexes=if defined(clustered_scramble_vcf_index) then select_all([clustered_scramble_vcf_index]) else NONE_ARRAY_,
      clustered_wham_vcfs=if defined(clustered_wham_vcf) then select_all([clustered_wham_vcf]) else NONE_ARRAY_,
      clustered_wham_vcf_indexes=if defined(clustered_wham_vcf_index) then select_all([clustered_wham_vcf_index]) else NONE_ARRAY_,
      clustered_depth_vcfs=if defined(clustered_depth_vcf) then select_all([clustered_depth_vcf]) else NONE_ARRAY_,
      clustered_depth_vcf_indexes=if defined(clustered_depth_vcf_index) then select_all([clustered_depth_vcf_index]) else NONE_ARRAY_,
      ped_file=ped_file,
      contig_list=contig_list,
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_fasta_fai,
      reference_dict=reference_dict,
      chr_x=chr_x,
      chr_y=chr_y,
      gatk_docker=gatk_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
  }

  call svc.SVConcordance {
    input:
      eval_vcf=RefineComplexVariants.cpx_refined_vcf,
      truth_vcf=JoinRawCalls.joined_raw_calls_vcf,
      output_prefix=sample_id,
      contig_list=contig_list,
      reference_dict=reference_dict,
      gatk_docker=gatk_docker,
      sv_base_mini_docker=sv_base_mini_docker
  }

  call fg.FilterGenotypes {
    input:
      vcf=SVConcordance.concordance_vcf,
      output_prefix=sample_id,
      ploidy_table=JoinRawCalls.ploidy_table,
      gq_recalibrator_model_file=gq_recalibrator_model_file,
      recalibrate_gq_args=recalibrate_gq_args,
      genome_tracks=genome_tracks,
      no_call_rate_cutoff=no_call_rate_cutoff,
      sl_filter_args=sl_filter_args,
      run_qc=run_main_vcf_qc,
      primary_contigs_fai=primary_contigs_fai,
      linux_docker=linux_docker,
      gatk_docker=gq_recalibrator_gatk_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker
  }

  Array[String] ref_samples = read_lines(ref_samples_list)
  call SingleSampleMetrics.SingleSampleMetrics as SampleFilterMetrics {
    input:
      name = batch,
      ref_samples = ref_samples,
      case_sample = sample_id,
      wgd_scores = wgd_scores,
      sample_counts = case_counts_file,
      contig_list = contig_list,
      linux_docker = linux_docker,
      sv_pipeline_docker = sv_pipeline_docker
  }

  call utils.RunQC as SampleFilterQC {
    input:
      name=batch,
      metrics=SampleFilterMetrics.metrics_file,
      qc_definitions = qc_definitions,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call SingleSampleFiltering.SampleQC as FilterSample {
    input:
      vcf=FilterGenotypes.filtered_vcf,
      sample_filtering_qc_file=SampleFilterQC.out,
      sv_pipeline_docker=sv_pipeline_docker,
  }

  call annotate.AnnotateVcf {
    input:
      vcf = FilterSample.out,
      prefix = batch,
      contig_list = contig_list,
      protein_coding_gtf = protein_coding_gtf,
      noncoding_bed = noncoding_bed,
      promoter_window = promoter_window,
      max_breakend_as_cnv_length = max_breakend_as_cnv_length,
      external_af_ref_bed = external_af_ref_bed,
      external_af_ref_prefix = external_af_ref_bed_prefix,
      external_af_population = external_af_population,
      sv_per_shard = annotation_sv_per_shard,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      gatk_docker = gatk_docker
  }

  call SingleSampleFiltering.VcfToBed as VcfToBed {
    input:
      vcf = AnnotateVcf.annotated_vcf,
      prefix = batch,
      sv_pipeline_docker = sv_pipeline_docker
  }

  call SingleSampleFiltering.UpdateBreakendRepresentation {
    input:
      vcf=AnnotateVcf.annotated_vcf,
      vcf_idx=AnnotateVcf.annotated_vcf_index,
      ref_fasta=reference_fasta,
      ref_fasta_idx=reference_fasta_fai,
      prefix=basename(AnnotateVcf.annotated_vcf, ".vcf.gz") + ".final_cleanup",
      sv_pipeline_docker=sv_pipeline_docker
  }

  call SingleSampleMetrics.SingleSampleMetrics {
    input:
      name = batch,
      ref_samples = ref_samples,
      case_sample = sample_id,
      wgd_scores = wgd_scores,
      sample_pe = case_pe_file,
      sample_sr = case_sr_file,
      sample_counts = case_counts_file,
      cleaned_vcf = cleaned_vcf,
      final_vcf = UpdateBreakendRepresentation.out,
      genotyped_pesr_vcf = convert_cnvs_without_depth_support_vcf,
      genotyped_depth_vcf = genotyped_depth_vcf,
      non_genotyped_unique_depth_calls_vcf = non_genotyped_unique_depth_calls_vcf,
      contig_list = contig_list,
      linux_docker = linux_docker,
      sv_pipeline_docker = sv_pipeline_docker
  }

  call utils.RunQC as SingleSampleQC {
    input:
      name = batch,
      metrics = SingleSampleMetrics.metrics_file,
      qc_definitions = qc_definitions,
      sv_pipeline_docker = sv_pipeline_docker
  }

  output {
      File final_vcf = UpdateBreakendRepresentation.out
      File final_vcf_idx = UpdateBreakendRepresentation.out_idx

      File final_bed = VcfToBed.bed

      # These files contain events reported in the internal VCF representation
      # They are less VCF-spec compliant but may be useful if components of the pipeline need to be re-run
      # on the output.
      File pre_cleanup_vcf = AnnotateVcf.annotated_vcf
      File pre_cleanup_vcf_idx = AnnotateVcf.annotated_vcf_index

      File metrics_file = SingleSampleMetrics.metrics_file
      File qc_file = SingleSampleQC.out
  }
}
