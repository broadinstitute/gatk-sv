version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "TasksBenchmark.wdl" as tasks10

workflow RdPeSrAnno {
  input {
    String prefix
    String bam_or_cram_file
    String bam_or_cram_index
    File vcf_file
    File ref_fasta
    File ref_fai
    File ref_dict
    File contig_list
    File pe_metrics
    File sr_metrics
    File rd_metrics
    String pesrrd_annotation_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_rdpesr
    RuntimeAttr? runtime_attr_bcf2vcf
    RuntimeAttr? runtime_attr_LocalizeCram
    RuntimeAttr? runtime_attr_vcf2bed
    RuntimeAttr? runtime_attr_SplitVcf
    RuntimeAttr? runtime_attr_ConcatBeds
  }

  Array[String] contigs = transpose(read_tsv(contig_list))[0]
  scatter ( contig in contigs ) {

    call tasks10.LocalizeCram as LocalizeCram{
      input:
        contig = contig,
        ref_fasta=ref_fasta,
        ref_fai=ref_fai,
        ref_dict=ref_dict,
        bam_or_cram_file=bam_or_cram_file,
        bam_or_cram_index=bam_or_cram_index,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_LocalizeCram
      }

    call tasks10.SplitVcf as SplitVcf{
      input:
        contig = contig,
        vcf_file = vcf_file,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_SplitVcf
      }

    call tasks10.vcf2bed as vcf2bed{
      input:
        vcf = SplitVcf.contig_vcf,
        vcf_index = SplitVcf.contig_vcf_index,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_vcf2bed
    }

    call RunRdPeSrAnnotation{
      input:
        prefix = prefix,
        contig = contig,
        bam_or_cram_file=LocalizeCram.local_bam,
        bam_or_cram_index=LocalizeCram.local_bai,
        bed = vcf2bed.bed,
        pe_metrics = pe_metrics,
        sr_metrics = sr_metrics,
        rd_metrics = rd_metrics,
        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        ref_dict=ref_dict,
        pesrrd_annotation_docker = pesrrd_annotation_docker,
        runtime_attr_override = runtime_attr_rdpesr
      }
  }

  call MiniTasks.ConcatBeds as ConcatPesrAnno{
    input:
      shard_bed_files=RunRdPeSrAnnotation.pesr_anno,
      prefix=prefix,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_ConcatBeds
      }

  call MiniTasks.ConcatBeds as ConcatRdAnno{
    input:
      shard_bed_files=RunRdPeSrAnnotation.cov,
      prefix=prefix,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_ConcatBeds
      }

  call MiniTasks.ConcatBeds as ConcatRdAnnoLeFlank{
    input:
      shard_bed_files=RunRdPeSrAnnotation.cov_le_flank,
      prefix=prefix,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_ConcatBeds
      }

  call MiniTasks.ConcatBeds as ConcatRdAnnoRiFlank{
    input:
      shard_bed_files=RunRdPeSrAnnotation.cov_ri_flank,
      prefix=prefix,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_ConcatBeds
      }


  output{
      File PesrAnno = ConcatPesrAnno.merged_bed_file
      File RdAnno = ConcatRdAnno.merged_bed_file
      File RdAnnoLeFlank = ConcatRdAnnoLeFlank.merged_bed_file
      File RdAnnoRiFlank = ConcatRdAnnoRiFlank.merged_bed_file
    }
  }


task RunRdPeSrAnnotation{
  input{
    String prefix
    String contig
    File bam_or_cram_file
    File bam_or_cram_index
    File bed
    File pe_metrics
    File sr_metrics
    File rd_metrics
    File ref_fasta
    File ref_fai
    File ref_dict
    String pesrrd_annotation_docker
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
  String filename = basename(bed, '.bed')

  output {
    File pesr_anno = "~{bed}.INS_with_SR_PE"
    File cov = "~{bed}.Seq_Cov"
    File cov_ri_flank = "~{bed}.ri_flank.Seq_Cov"
    File cov_le_flank = "~{bed}.le_flank.Seq_Cov"
  }

  command <<<

    set -Eeuo pipefail
    Rscript /src/modify_bed_for_PE_SR_RD_labeling.R -i ~{bed}

    python3 /src/add_SR_PE_to_PB_INS.V2.py ~{bed} ~{pe_metrics} ~{sr_metrics} 

    zcat ~{rd_metrics} | grep -v '@' | grep -v CONTIG |bgzip >  bincov.tsv.gz
    Rscript /src/bincov_to_normCov.R -i bincov.tsv.gz
    bgzip normCov.tsv
    tabix normCov.tsv.gz

    python3 /src/add_RD_to_SVs.py ~{bed} normCov.tsv.gz
    python3 /src/add_RD_to_SVs.py ~{filename}.ri_flank normCov.tsv.gz
    python3 /src/add_RD_to_SVs.py ~{filename}.le_flank normCov.tsv.gz

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: pesrrd_annotation_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
