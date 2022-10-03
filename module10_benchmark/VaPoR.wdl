##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

## Copyright Broad Institute, 2020
## 
## This WDL pipeline implements VaPoR 
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
import "VaPoRBed.wdl" as vapor_bed
import "VaPoRVcf.wdl" as vapor_vcf
workflow VaPoR{
    input{
        String prefix
        String sample
        String bam_or_cram_file
        String bam_or_cram_index
        File? vcf_file
        File? bed_file
        File ref_fasta
        File ref_fai
        File ref_dict
        Int min_shard_size
        String vapor_docker
        String sv_base_mini_docker
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_vapor 
        RuntimeAttr? runtime_attr_bcf2vcf
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_SplitVcf
        RuntimeAttr? runtime_attr_ConcatBeds
    }

    if (defined(vcf_file)) {
        call ExtractContigList as extract_contig_vcf{
            input:
                input_file = vcf_file,
                sv_base_mini_docker = sv_base_mini_docker

        }

        Array[String] contigs_vcf = transpose(read_tsv(extract_contig_vcf.contig_file))[0]

        call vapor_vcf.VaPoRVcf as VaPoR_vcf{
            input:
                prefix = prefix,
                sample =  sample,
                bam_or_cram_file = bam_or_cram_file,
                bam_or_cram_index = bam_or_cram_index,
                vcf_file = vcf_file,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                contigs = contigs_vcf,
                min_shard_size = min_shard_size,
                vapor_docker = vapor_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_vapor    = runtime_attr_vapor,
                runtime_attr_bcf2vcf = runtime_attr_bcf2vcf,
                runtime_attr_vcf2bed = runtime_attr_vcf2bed,
                runtime_attr_SplitVcf = runtime_attr_SplitVcf,
                runtime_attr_ConcatBeds = runtime_attr_ConcatBeds
        }
    }

    if (defined(bed_file)) {
        call ExtractContigList as extract_contig_bed{
            input:
                input_file = bed_file,
                sv_base_mini_docker = sv_base_mini_docker

        }

        Array[String] contigs_bed = transpose(read_tsv(extract_contig_bed.contig_file))[0]

        call vapor_bed.VaPoRBed as VaPoR_bed{
            input:
                prefix = prefix,
                sample = sample,
                bam_or_cram_file = bam_or_cram_file,
                bam_or_cram_index = bam_or_cram_index,
                bed_file = bed_file,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                contigs = contigs_bed,
                min_shard_size = min_shard_size,
                vapor_docker = vapor_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_vapor    = runtime_attr_vapor,
                runtime_attr_bcf2vcf = runtime_attr_bcf2vcf,
                runtime_attr_vcf2bed = runtime_attr_vcf2bed,
                runtime_attr_SplitVcf = runtime_attr_SplitVcf,
                runtime_attr_ConcatBeds = runtime_attr_ConcatBeds
        }
    }
    
    output{
        File vcf_out = select_first([VaPoR_vcf.bed, VaPoR_bed.bed])
        File bed_out = select_first([VaPoR_vcf.bed, VaPoR_bed.bed])
        File vcf_plots = select_first([VaPoR_vcf.plots, VaPoR_bed.plots])
        File bed_plots = select_first([VaPoR_vcf.plots, VaPoR_bed.plots])
    }
}

task ExtractContigList{
  input{
    File? input_file
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 1, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File contig_file = "contig.tsv"
  }


  command <<<

    set -Eeuo pipefail

    if [[ ~{input_file} == *.gz ]] ;  then
      zcat ~{input_file} | grep -v "#" | cut -f1 | sort | uniq  > "contig.tsv"
    else
      grep -v "#" ~{input_file} | cut -f1 | sort | uniq > "contig.tsv"
    fi

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


