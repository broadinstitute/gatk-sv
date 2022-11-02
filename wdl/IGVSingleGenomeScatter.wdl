version 1.0

import "Structs.wdl"
import "IGVSingleGenome.wdl" as igv_sg

workflow IGVSingleGenome {
  input{
    # Bed files containing regions to screenshot; 4th column must be SVID
    Array[File] beds
    Array[String] run_names

    # Sample reads
    File bam_or_cram
    File bam_or_cram_index

    # Sample id and prefix for output filenames
    String sample_id

    # Reference corresponding to read alignments
    File ref_fasta
    File ref_fai

    Int? records_per_shard

    String linux_docker
    String igv_docker

    RuntimeAttr? runtime_attr_scatter
    RuntimeAttr? runtime_attr_igv
  }

  scatter (i in range(length(beds))) {
    call igv_sg.IGVSingleGenome {
      input:
        bed=beds[i],
        run_name=run_names[i],
        bam_or_cram=bam_or_cram,
        bam_or_cram_index=bam_or_cram_index,
        sample_id=sample_id,
        ref_fasta=ref_fasta,
        ref_fai=ref_fai,
        records_per_shard=records_per_shard,
        linux_docker=linux_docker,
        igv_docker=igv_docker,
        runtime_attr_scatter=runtime_attr_scatter,
        runtime_attr_igv=runtime_attr_igv
    }
  }

  output{
    Array[Array[File]] igv_plots_tarballs_scattered = IGVSingleGenome.igv_plots_tarballs
    Array[Array[File]] igv_scripts_scattered = IGVSingleGenome.igv_scripts
  }
}
