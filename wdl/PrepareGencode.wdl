version 1.0

import "Structs.wdl"

# Prepare input GTF according to other input dictionaries and output GTF for annotation sub-module
workflow PrepareGencode {

  input {

    File gencode_annotation_gtf         # Gencode annotation GTF
    File gencode_pc_translations_fa     # Gencode protein-coding translation fasta
    File gencode_pc_transcripts_fa      # Gencode protein-coding transcript fasta
    File gencode_transcript_source      # Gencode transcript source metadata
    Int  promoter_window                # Window upstream of TSS to consider as promoter region

    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_get_canonical_transcripts
    RuntimeAttr? runtime_attr_make_canonical_gtf
    RuntimeAttr? runtime_attr_make_promoters
    RuntimeAttr? runtime_attr_subset_gtf
  }

  call GetCanonicalTranscripts {
    input:

      gencode_annotation_gtf     = gencode_annotation_gtf,
      gencode_pc_translations_fa = gencode_pc_translations_fa,
      gencode_pc_transcripts_fa  = gencode_pc_transcripts_fa,
      gencode_transcript_source  = gencode_transcript_source,

      sv_pipeline_docker = sv_pipeline_docker,

      runtime_attr_override      = runtime_attr_get_canonical_transcripts
  }

  call MakeCanonicalGtf {
    input:

      gencode_annotation_gtf = gencode_annotation_gtf,
      canonical_transcripts  = GetCanonicalTranscripts.canon_tx,

      sv_pipeline_docker = sv_pipeline_docker,

      runtime_attr_override  = runtime_attr_make_canonical_gtf
  }

  call MakePromoters {
    input:

      gtf    = MakeCanonicalGtf.canon_gtf,
      window = promoter_window,

      sv_pipeline_docker = sv_pipeline_docker,

      runtime_attr_override = runtime_attr_make_promoters
  }

  call SubsetGtf as SubsetAntisense {
    input:

      annotation_gtf = gencode_annotation_gtf,
      subset         = "antisense",
      
      sv_base_mini_docker = sv_base_mini_docker,

      runtime_attr_override = runtime_attr_subset_gtf
  }

  call SubsetGtf as SubsetLincRNA {
    input:

      annotation_gtf = gencode_annotation_gtf,
      subset         = "lincRNA",

      sv_base_mini_docker = sv_base_mini_docker,

      runtime_attr_override = runtime_attr_subset_gtf
  }

  call SubsetGtf as SubsetProcessedTranscript {
    input:

      annotation_gtf = gencode_annotation_gtf,
      subset         = "processed_transcript",

      sv_base_mini_docker = sv_base_mini_docker,

      runtime_attr_override = runtime_attr_subset_gtf
  }

  call SubsetGtf as SubsetPseudogene {
    input:

      annotation_gtf = gencode_annotation_gtf,
      subset         = "pseudogene",

      sv_base_mini_docker = sv_base_mini_docker,

      runtime_attr_override = runtime_attr_subset_gtf
  }

  output {
    File canonical_gtf            = MakeCanonicalGtf.canon_gtf
    File canonical_promoters      = MakePromoters.promoter_bed
    File antisense_gtf            = SubsetAntisense.gtf
    File lincRNA_gtf              = SubsetLincRNA.gtf
    File processed_transcript_gtf = SubsetProcessedTranscript.gtf
    File pseudogene_gtf           = SubsetPseudogene.gtf
  }
}

task SubsetGtf {

  input {

    File   annotation_gtf
    String subset

    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_override
  }

  output {
    File gtf = "gencode.${subset}.gtf.gz"
  }

  command <<<

    set -euo pipefail

    zcat ~{annotation_gtf} \
      | grep -e "gene_type \"~{subset}\"" - \
      | sort -k1,1V -k4,4n \
      | bgzip -c \
      > gencode.~{subset}.gtf.gz
  
  >>>

  #########################
  RuntimeAttr default_attr = object {
    cpu_cores:          1, 
    mem_gb:             3.75, 
    disk_gb:            10,
    boot_disk_gb:       10,
    preemptible_tries:  3,
    max_retries:        0
    # docker:             "gatksv/sv-base-mini:v0.1",
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
    docker:                 "gatksv/sv-base-mini:v0.1" #select_first([runtime_attr.docker,            default_attr.docker])
  }
}

task GetCanonicalTranscripts {

  input {

    File gencode_annotation_gtf
    File gencode_pc_translations_fa
    File gencode_pc_transcripts_fa
    File gencode_transcript_source

    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_override
  }
  
  output {
    File canon_tx = "gencode.canonical_transcripts.txt"
  }

  command <<<

    set -euo pipefail
    /opt/sv-pipeline/05_annotation/scripts/get_canonical_transcripts.py \
      ~{gencode_annotation_gtf} \
      ~{gencode_pc_translations_fa} \
      ~{gencode_pc_transcripts_fa} \
      ~{gencode_transcript_source} \
      gencode.canonical_transcripts.txt
  
  >>>
  
  #########################
  RuntimeAttr default_attr = object {
    cpu_cores:          1, 
    mem_gb:             3.75, 
    disk_gb:            10,
    boot_disk_gb:       10,
    preemptible_tries:  3,
    max_retries:        0
    # docker:             "gatksv/sv-pipeline:v0.1",
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
    docker:                 sv_pipeline_docker #select_first([runtime_attr.docker,            default_attr.docker])
  }
}

task MakeCanonicalGtf {
  input {
    File gencode_annotation_gtf
    File canonical_transcripts

    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_override
  }
  
  output {
    File canon_gtf = "gencode.canonical_pc.gtf.gz"
  }
  
  command <<<

    set -euo pipefail

    cat \
      <(cut -f2 ~{canonical_transcripts} | sed -e '1d' | fgrep -w -f - <(zcat ~{gencode_annotation_gtf})) \
      <(cut -f1 ~{canonical_transcripts} | sed -e '1d' | fgrep -w -f - <(zcat ~{gencode_annotation_gtf}) | awk '($3=="gene")') \
      | sort -k10,10 \
      | /opt/sv-pipeline/05_annotation/scripts/filter_UTRs.py stdin stdout \
      | sort -k1,1V -k4,4n \
      | bgzip -c \
      > gencode.canonical_pc.gtf.gz
  
  >>>
  
  #########################
  RuntimeAttr default_attr = object {
    cpu_cores:          1, 
    mem_gb:             3.75, 
    disk_gb:            10,
    boot_disk_gb:       10,
    preemptible_tries:  3,
    max_retries:        0
    # docker:             "gatksv/sv-pipeline:v0.1",
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
    docker:                 sv_pipeline_docker #select_first([runtime_attr.docker,            default_attr.docker])
  }
}

task MakePromoters {
  
  input {
    
    File gtf
    Int  window

    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_override
  }
  
  output {
    File promoter_bed = "gencode.canonical_pc.promoter.bed"
  }
  
  command <<<

    /opt/sv-pipeline/05_annotation/scripts/make_promoters.sh ~{gtf} ~{window} > gencode.canonical_pc.promoter.bed
  
  >>>
  
  #########################
  RuntimeAttr default_attr = object {
    cpu_cores:          1, 
    mem_gb:             3.75, 
    disk_gb:            10,
    boot_disk_gb:       10,
    preemptible_tries:  3,
    max_retries:        0
    # docker:             "gatksv/sv-pipeline:v0.1",
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
    docker:                 sv_pipeline_docker #select_first([runtime_attr.docker,            default_attr.docker])
  }
}
