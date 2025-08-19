version 1.0

import "GatherSampleEvidence.wdl" as gse
import "Scramble.wdl" as scramble

workflow ScrambleRealignSR {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    Boolean move_bam_or_cram_files = false

    String sample_id

    # Required to do standalone scramble SR realignment
    File manta_vcf
    File counts
    File scramble_table_input
    Boolean keep_realigned_reads = false

    # Common parameters
    File primary_contigs_list
    File reference_fasta
    File reference_index    # Index (.fai), must be in same dir as fasta
    File reference_dict     # Dictionary (.dict), must be in same dir as fasta
    String? reference_version   # Either "38" or "19"

    # Reference bwa index files, only required for alignments with Dragen 3.7.8
    File? reference_bwa_alt
    File? reference_bwa_amb
    File? reference_bwa_ann
    File? reference_bwa_bwt
    File? reference_bwa_pac
    File? reference_bwa_sa

    # Scramble inputs
    File mei_bed
    Int? scramble_alignment_score_cutoff
    Int? scramble_percent_align_cutoff
    Float? scramble_min_clipped_reads_fraction
    Int? scramble_part2_threads
    File? scramble_vcf_script

    # Docker
    String sv_pipeline_docker
    String sv_base_mini_docker
    String? scramble_docker


    # Runtime configuration overrides
    RuntimeAttr? runtime_attr_localize_reads
    RuntimeAttr? runtime_attr_scramble_part1
    RuntimeAttr? runtime_attr_scramble_part2
    RuntimeAttr? runtime_attr_scramble_make_vcf
    RuntimeAttr? runtime_attr_realign_soft_clips
    RuntimeAttr? runtime_attr_scramble_part1_realigned
    RuntimeAttr? runtime_attr_scramble_part2_realigned
    RuntimeAttr? runtime_attr_scramble_make_vcf_realigned


    # Never assign these values! (workaround until None type is implemented)
    Float? NONE_FLOAT_
    Int? NONE_INT_
    File? NONE_FILE_

  }

  call gse.LocalizeReads {
    input:
      reads_path = bam_or_cram_file,
      reads_index = bam_or_cram_index,
      move_files = move_bam_or_cram_files,
      runtime_attr_override = runtime_attr_localize_reads
  }

  Int scramble_alignment_score_cutoff_ = select_first([scramble_alignment_score_cutoff, 60])
  call gse.RealignSoftClippedReads {
    input:
      reads_path = LocalizeReads.output_file,
      reads_index = LocalizeReads.output_index,
      scramble_table = scramble_table_input,
      is_bam = false,
      sample_id = sample_id,
      reference_fasta = reference_fasta,
      reference_index = reference_index,
      reference_bwa_alt = select_first([reference_bwa_alt]),
      reference_bwa_amb = select_first([reference_bwa_amb]),
      reference_bwa_ann = select_first([reference_bwa_ann]),
      reference_bwa_bwt = select_first([reference_bwa_bwt]),
      reference_bwa_pac = select_first([reference_bwa_pac]),
      reference_bwa_sa = select_first([reference_bwa_sa]),
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_realign_soft_clips
  }
  call scramble.Scramble as ScrambleRealigned {
    input:
      bam_or_cram_file = RealignSoftClippedReads.out,
      bam_or_cram_index = RealignSoftClippedReads.out_index,
      original_bam_or_cram_file = LocalizeReads.output_file,
      original_bam_or_cram_index = LocalizeReads.output_index,
      counts_file = counts,
      input_vcf = manta_vcf,
      sample_name = sample_id,
      reference_fasta = reference_fasta,
      reference_index = reference_index,
      mei_bed = mei_bed,
      regions_list = primary_contigs_list,
      alignment_score_cutoff = scramble_alignment_score_cutoff_,
      percent_align_cutoff = scramble_percent_align_cutoff,
      min_clipped_reads_fraction = scramble_min_clipped_reads_fraction,
      part2_threads = scramble_part2_threads,
      scramble_vcf_script = scramble_vcf_script,
      scramble_docker = select_first([scramble_docker]),
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_scramble_part1 = runtime_attr_scramble_part1_realigned,
      runtime_attr_scramble_part2 = runtime_attr_scramble_part2_realigned,
      runtime_attr_scramble_make_vcf = runtime_attr_scramble_make_vcf_realigned
  }

  output {
    File scramble_vcf = ScrambleRealigned.vcf
    File scramble_index = ScrambleRealigned.index
    File scramble_clusters = ScrambleRealigned.clusters
    File scramble_table = ScrambleRealigned.table

    File? realigned_reads = if (keep_realigned_reads) then RealignSoftClippedReads.out else NONE_FILE_
    File? realigned_reads_index = if (keep_realigned_reads) then RealignSoftClippedReads.out_index else NONE_FILE_

  }

}
