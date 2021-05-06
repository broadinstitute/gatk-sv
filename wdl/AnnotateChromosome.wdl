version 1.0

import "Structs.wdl"

# After VCF and dictionary are scattered per chromosome, apply annotation to scattered per-chromosome VCF
workflow AnnotateChromosome {

  input {

    String prefix
    File  vcf
    File  protein_coding_gtf
    File  linc_rna_gtf
    File  promoter_bed
    File  noncoding_bed
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_annotate_intervals
    RuntimeAttr? runtime_attr_merge_annotations
  }

  String coding_flag    = "--gencode"
  String noncoding_flag = "--noncoding"

  call AnnotateIntervals as AnnotateProteinCoding {
    input:
      vcf                = vcf,
      intervals_file     = protein_coding_gtf,
      intervals_flag     = coding_flag,
      prefix             = prefix,
      intervals_set      = "protein_coding",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_annotate_intervals
  }

  call AnnotateIntervals as AnnotateLincRNA {
    input:
      vcf                = vcf,
      intervals_file     = linc_rna_gtf,
      intervals_flag     = coding_flag,
      prefix             = prefix,
      intervals_set      = "lincRNA",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_annotate_intervals
  }

  call AnnotateIntervals as AnnotatePromoters {
    input:
      vcf                = vcf,
      intervals_file     = promoter_bed,
      intervals_flag     = noncoding_flag,
      prefix             = prefix,
      intervals_set      = "promoters",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_annotate_intervals
  }

  call AnnotateIntervals as annotate_noncoding_elements {
    input:
      vcf = vcf,
      intervals_file = noncoding_bed,
      intervals_flag = noncoding_flag,
      prefix = prefix,
      intervals_set = "noncoding",
      sv_pipeline_docker  =  sv_pipeline_docker,
      runtime_attr_override = runtime_attr_annotate_intervals
  }
  


  call MergeAnnotations {
    input:
      vcf                = vcf,
      annotated_vcfs     = [AnnotateProteinCoding.annotated_vcf, AnnotateLincRNA.annotated_vcf, AnnotatePromoters.annotated_vcf,
      annotate_noncoding_elements.annotated_vcf],
      prefix             = prefix,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_merge_annotations
  }
  
  output {
    File annotated_vcf = MergeAnnotations.annotated_vcf
    File annotated_vcf_idx = MergeAnnotations.annotated_vcf_idx
  }
}

# Apply annoattion
task AnnotateIntervals {

  input {

    File   vcf
    File   intervals_file   # gtf or bed file
    String intervals_flag   # "--gencode" or "--noncoding"
    String prefix
    String intervals_set
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_override
  }

  output {
    File annotated_vcf = "${prefix}.${intervals_set}.vcf.gz"
  }

  command <<<

    set -euo pipefail

    svtk annotate \
      ~{intervals_flag} ~{intervals_file} \
      ~{vcf} \
      ~{prefix}.~{intervals_set}.vcf
    
    orig=$( zcat ~{vcf} | cut -f1 | grep -cv "^#" || true )
    new=$( cut -f1 ~{prefix}.~{intervals_set}.vcf | grep -cv "^#" || true )
    
    if [ "$new" -ne "$orig" ]; then
      echo "ANNOTATED VCF DOES NOT HAVE THE SAME NUMBER OF RECORDS AS INPUT VCF ($new vs $orig)"
      exit 1
    fi
    
    bgzip -f ~{prefix}.~{intervals_set}.vcf
  
  >>>
  
  #########################
  RuntimeAttr default_attr = object {
    cpu_cores:          1, 
    mem_gb:             3.75, 
    disk_gb:            50,
    boot_disk_gb:       10,
    preemptible_tries:  3,
    max_retries:        0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
    docker:                 sv_pipeline_docker
  }
}

# Merge annotations from different dictionaries (genocode, lincRNA, promoters, potentially more in the future)
task MergeAnnotations {

  input {

    File        vcf
    Array[File] annotated_vcfs
    String      prefix
    String      sv_pipeline_docker

    RuntimeAttr? runtime_attr_override
  }

  output {
    File annotated_vcf = "${prefix}.annotated.vcf.gz"
    File annotated_vcf_idx = "${prefix}.annotated.vcf.gz.tbi"
  }

  command <<<

    set -euo pipefail

    /opt/sv-pipeline/05_annotation/scripts/merge_annotations.py \
      ~{vcf} \
      ~{sep=" " annotated_vcfs} \
      ~{prefix}.annotated.vcf

    bgzip ~{prefix}.annotated.vcf
    tabix ~{prefix}.annotated.vcf.gz
    
    orig=$( zcat ~{vcf} | cut -f1 | grep -cv "^#" || true)
    new=$( zcat ~{prefix}.annotated.vcf.gz | cut -f1 | grep -cv "^#" || true)
    if [ "$new" -ne "$orig" ]; then
      echo "ANNOTATED VCF DOES NOT HAVE THE SAME NUMBER OF RECORDS AS INPUT VCF ($new vs $orig)"
      exit 1
    fi
  
  >>>

  #########################
  RuntimeAttr default_attr = object {
    cpu_cores:          1, 
    mem_gb:             8, 
    disk_gb:            250,
    boot_disk_gb:       10,
    preemptible_tries:  3,
    max_retries:        0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
    docker:                 sv_pipeline_docker
  }
}
