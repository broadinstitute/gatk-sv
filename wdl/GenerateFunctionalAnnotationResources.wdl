version 1.0

import "PrepareGencode.wdl" as pg
# import "PrepareNoncoding.wdl" as pn

# Workflow for preprocessing for functional annotation
workflow GenerateFunctionalAnnotationResources {
  input {
    
    ### args for PrepareGencode
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
    
#     ### args for PrepareNoncoding
#     File noncoding_bed_list
# 
#     RuntimeAttr? runtime_attr_clean_noncoding_bed
#     RuntimeAttr? runtime_attr_make_noncoding_bed
  }
  
  call pg.PrepareGencode as PrepareGencode {
    input:
      
      gencode_annotation_gtf     = gencode_annotation_gtf,
      gencode_pc_translations_fa = gencode_pc_translations_fa,
      gencode_pc_transcripts_fa  = gencode_pc_transcripts_fa,
      gencode_transcript_source  = gencode_transcript_source,
      promoter_window            = promoter_window,

      sv_base_mini_docker = sv_base_mini_docker,
      sv_pipeline_docker  = sv_pipeline_docker,
      
      runtime_attr_get_canonical_transcripts = runtime_attr_get_canonical_transcripts,
      runtime_attr_make_canonical_gtf        = runtime_attr_make_canonical_gtf,
      runtime_attr_make_promoters            = runtime_attr_make_promoters,
      runtime_attr_subset_gtf                = runtime_attr_subset_gtf
  }

#   call pn.PrepareNoncoding as PrepareNoncoding {
#     input:
#       
#       noncoding_bed_list = noncoding_bed_list,
# 
#       sv_base_mini_docker = sv_base_mini_docker,
# 
#       runtime_attr_clean_noncoding_bed = runtime_attr_clean_noncoding_bed,
#       runtime_attr_make_noncoding_bed  = runtime_attr_make_noncoding_bed
#   }

  output {
    
    File canonical_gtf            = PrepareGencode.canonical_gtf
    File canonical_promoters      = PrepareGencode.canonical_promoters
    File antisense_gtf            = PrepareGencode.antisense_gtf
    File lincRNA_gtf              = PrepareGencode.lincRNA_gtf
    File processed_transcript_gtf = PrepareGencode.processed_transcript_gtf
    File pseudogene_gtf           = PrepareGencode.pseudogene_gtf

    # File noncoding_bed            = PrepareNoncoding.noncoding_bed
  }
}
