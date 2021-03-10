##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

## Copyright Broad Institute, 2020
## 
## This WDL pipeline implements Duphold 
##
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

version 1.0

import "Structs.wdl"
import "TasksBenchmark.wdl" as mini_tasks
workflow VaPoR{
  input{
    String prefix
    String bam_or_cram_file
    String bam_or_cram_index
    File? vcf_file
    File? bed_file
    File ref_fasta
    File ref_fai
    File ref_dict
    Array[String] contigs
    String vapor_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_Vapor 
    RuntimeAttr? runtime_attr_bcf2vcf
    RuntimeAttr? runtime_attr_vcf2bed
    RuntimeAttr? runtime_attr_SplitVcf
    RuntimeAttr? runtime_attr_ConcatBeds
  }

  scatter ( contig in contigs ) {

    #call mini_tasks.LocalizeCram as LocalizeCram{
    #  input:
    #    contig = contig,
    #    ref_fasta=ref_fasta,
    #    ref_fai=ref_fai,
    #    ref_dict=ref_dict,
    #    bam_or_cram_file=bam_or_cram_file,
    #    bam_or_cram_index=bam_or_cram_index,
    #    sv_pipeline_docker=sv_pipeline_docker,
    #    runtime_attr_override=runtime_attr_LocalizeCram
    #  }

    if (defined(vcf_file)) {
      call mini_tasks.SplitVcf as SplitVcf{
        input:
          contig = contig,
          vcf_file = vcf_file,
          sv_pipeline_docker=sv_pipeline_docker,
          runtime_attr_override=runtime_attr_SplitVcf
      }

      call mini_tasks.vcf2bed as vcf2bed{
        input:
          vcf = SplitVcf.contig_vcf,
          vcf_index = SplitVcf.contig_vcf_index,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_vcf2bed
      }

      File contig_bed = vcf2bed.bed
   }

    if (defined(bed_file)){
      call mini_tasks.SplitBed as SplitBed{
        input:
          contig = contig,
          bed_file = bed_file,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override=runtime_attr_SplitVcf
      }
      File contig_bed = SplitBed.contig_bed
    }

    call RunVaPoRWithCram as RunVaPoR{
      input:
        prefix = prefix,
        contig = contig,
        bam_or_cram_file=bam_or_cram_file,
        bam_or_cram_index=bam_or_cram_index,
        bed = contig_bed,
        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        vapor_docker = vapor_docker,
        runtime_attr_override = runtime_attr_Vapor
    }
  }

  call mini_tasks.ConcatVaPoR as ConcatVaPoR{
    input:
      shard_bed_files=RunVaPoR.vapor,
      prefix=prefix,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_ConcatBeds
      }

  output{
      File bed = ConcatVaPoR.merged_bed_file
      Array[File] plots = RunVaPoR.vapor_plot
    }
  }


task RunVaPoR{
  input{
    String prefix
    String contig
    File bam_or_cram_file
    File bam_or_cram_index
    File bed
    File ref_fasta
    File ref_fai
    File ref_dict
    String vapor_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 5,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File vapor = "~{prefix}.~{contig}.vapor.gz"
    File vapor_plot = "~{prefix}.~{contig}.tar.gz"
  }

  command <<<

    set -Eeuo pipefail

    mkdir ~{prefix}.~{contig}
  
    vapor bed \
    --sv-input ~{bed} \
    --output-path ~{prefix}.~{contig} \
    --output-file ~{prefix}.~{contig}.vapor \
    --reference ~{ref_fasta} \
    --pacbio-input ~{bam_or_cram_file}

    tar -czf ~{prefix}.~{contig}.tar.gz ~{prefix}.~{contig}
    bgzip  ~{prefix}.~{contig}.vapor
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: vapor_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task RunVaPoRWithCram{
  input{
    String prefix
    String contig
    File bam_or_cram_file
    File bam_or_cram_index
    File bed
    File ref_fasta
    File ref_fai
    File ref_dict
    String vapor_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File vapor = "~{prefix}.~{contig}.vapor.gz"
    File vapor_plot = "~{prefix}.~{contig}.tar.gz"
  }

  command <<<

    set -Eeuo pipefail

    #localize cram files
    java -Xmx~{java_mem_mb}M -jar ${GATK_JAR}  PrintReads \
    -I ~{bam_or_cram_file} \
    -L ~{contig} \
    -O ~{contig}.bam \
    -R ~{ref_fasta}

    samtools index ~{contig}.bam
  
    #run vapor
    mkdir ~{prefix}.~{contig}

    vapor bed \
    --sv-input ~{bed} \
    --output-path ~{prefix}.~{contig} \
    --output-file ~{prefix}.~{contig}.vapor \
    --reference ~{ref_fasta} \
    --pacbio-input ~{contig}.bam

    tar -czf ~{prefix}.~{contig}.tar.gz ~{prefix}.~{contig}
    bgzip  ~{prefix}.~{contig}.vapor
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: vapor_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
