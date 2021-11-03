version 1.0

import "Structs.wdl"
import "ExpansionHunter.wdl" as ExpansionHunter

workflow ExpansionHunterScatter {

    input {
        Array[File] bams_or_crams
        Array[File]? bams_or_crams_indexes
        File reference_fasta
        File? reference_fasta_index
        File variant_catalog
        String eh_docker
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

        call ExpansionHunter.ExpansionHunter as expanionHunter {
            input:
                bam_or_cram=bam_or_cram_,
                bam_or_cram_index=bam_or_cram_index_,
                reference_fasta=reference_fasta,
                reference_fasta_index=reference_fasta_index_,
                variant_catalog=variant_catalog,
                docker_file=eh_docker,
                runtime_attr=runtime_attr
        }
    }

    output {
        Array[File] jsons = expanionHunter.json
        Array[File] vcfs = expanionHunter.vcf
        Array[File] overlapping_reads = expanionHunter.overlapping_reads
    }
}
