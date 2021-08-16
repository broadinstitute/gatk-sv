# Workflow to parallelize VCF annotation by chromosome

version 1.0

import "TasksMakeCohortVcf.wdl" as MiniTasks
import "AnnotateChromosome.wdl" as annotate_by_chrom

# Scatter VCF and apply prepared annotations
workflow ScatterAnnotateVcfByChrom {

  input {

    File vcf
    String prefix
    File vcf_idx
    File contig_list
    File protein_coding_gtf
    File linc_rna_gtf
    File promoter_bed
    File noncoding_bed

    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_annotate_intervals
    RuntimeAttr? runtime_attr_merge_annotations
    RuntimeAttr? runtime_attr_subset_vcf
    RuntimeAttr? runtime_attr_concat_vcfs
  }

  Array[Array[String]] contigs = read_tsv(contig_list)

  # Annotate, scattered by chromosome
  scatter (contig in contigs) {
    # Remote tabix each chromosome
    call SubsetVcf {
      input:
        vcf     = vcf,
        vcf_idx = vcf_idx,
        contig  = contig[0],
        prefix  = "${prefix}.${contig[0]}",
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_subset_vcf
    }

    # Annotate per chromosome
    call annotate_by_chrom.AnnotateChromosome as AnnotateChromosome {
      input:
        vcf                = SubsetVcf.subsetted_vcf,
        prefix             = "${prefix}.${contig[0]}",
        protein_coding_gtf = protein_coding_gtf,
        linc_rna_gtf       = linc_rna_gtf,
        promoter_bed       = promoter_bed,
        noncoding_bed      = noncoding_bed,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_annotate_intervals = runtime_attr_annotate_intervals,
        runtime_attr_merge_annotations  = runtime_attr_merge_annotations
    }
  }

  # Merge integrated vcfs across chromosomes
  call MiniTasks.ConcatVcfs as ConcatVcfs {
    input:
      vcfs = AnnotateChromosome.annotated_vcf,
      vcfs_idx = AnnotateChromosome.annotated_vcf_idx,
      outfile_prefix = "${prefix}.annotated",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_vcfs
  }

  output {
    File annotated_vcf = ConcatVcfs.concat_vcf
    File annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
  }
}

# Scatter VCF by chromosome
task SubsetVcf {

  input {

    File vcf
    File   vcf_idx
    String contig
    String prefix
    
    String sv_pipeline_docker
    
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    vcf: {
      localization_optional: true
    }
    vcf_idx: {
      localization_optional: true
    }
  }

  output {
    File subsetted_vcf     = "${prefix}.${contig}.vcf.gz"
    File subsetted_vcf_idx = "${prefix}.${contig}.vcf.gz.tbi"
  }

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

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  command <<<

    set -euo pipefail

    java -Xmx~{java_mem_mb}M -jar ${GATK_JAR} SelectVariants \
      -V "~{vcf}" \
      -L "~{contig}" \
      -O ~{prefix}.~{contig}.vcf.gz
  
  >>>

  runtime {
    cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:                 mem_gb + " GiB"
    disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
    docker:                 sv_pipeline_docker
  }
}
