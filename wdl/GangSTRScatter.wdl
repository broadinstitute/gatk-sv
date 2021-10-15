version 1.0

import "Structs.wdl"
import "GangSTR.wdl" as GangSTR

workflow GangSTRScatter {

  input {
    Array[File] bams_or_crams
    Array[File]? bams_or_crams_indexes
    File reference_fasta
    File? reference_fasta_index
    File target_tr_loci_regions_bed
    String str_docker
    RuntimeAttr? runtime_attr
  }

  scatter (i in range(length(bams_or_crams))) {
    File bam_or_cram_ = bams_or_crams[i]
    Boolean is_bam =
      basename(bam_or_cram_, ".bam") + ".bam" == basename(bam_or_cram_)
    File bam_or_cram_index_ =
      if defined(bams_or_crams_indexes) then
        select_first([bams_or_crams_indexes])[i]
      else
        bam_or_cram_ + if is_bam then ".bai" else ".crai"
    File reference_fasta_index_ = select_first([
      reference_fasta_index, reference_fasta + ".fai"])

    call GangSTR.GangSTR as gangSTR {
      input:
        bam_or_cram=bam_or_cram_,
        bam_or_cram_index=bam_or_cram_index_,
        reference_fasta=reference_fasta,
        reference_fasta_index=reference_fasta_index_,
        target_tr_loci_regions_bed=target_tr_loci_regions_bed,
        str_docker=str_docker,
        runtime_attr=runtime_attr
    }
  }

  output {
    Array[File] output_vcfs = gangSTR.output_vcf
    Array[File] samples_stats = gangSTR.sample_stats
    Array[File] insdatas = gangSTR.insdata
  }
}
