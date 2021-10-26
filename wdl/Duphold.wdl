version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "TasksBenchmark.wdl" as tasks10

workflow Duphold {
  input {
    String prefix
    String bam_or_cram_file
    String bam_or_cram_index
    File vcf_file
    File ref_fasta
    File ref_fai
    File ref_dict
    File contig_list
    String duphold_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_duphold
    RuntimeAttr? runtime_attr_bcf2vcf
    RuntimeAttr? runtime_attr_LocalizeCram
    RuntimeAttr? runtime_attr_SplitVcf
    RuntimeAttr? runtime_attr_ConcatVcfs
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

    call RunDupholdPerContig{
      input:
        prefix = prefix,
        contig = contig,
        bam_or_cram_file=LocalizeCram.local_bam,
        bam_or_cram_index=LocalizeCram.local_bai,
        vcf_file = SplitVcf.contig_vcf,
        vcf_index = SplitVcf.contig_vcf_index,
        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        duphold_docker = duphold_docker,
        runtime_attr_override = runtime_attr_duphold
      }

    call Bcf2Vcf{
      input:
        prefix = prefix,
        contig = contig,
        bcf = RunDupholdPerContig.bcf,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_bcf2vcf
      }
  }

  call MiniTasks.ConcatVcfs as ConcatVcfs{
    input:
      vcfs=Bcf2Vcf.vcf,
      outfile_prefix=prefix,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_ConcatVcfs
  }

    output{
      File vcf = ConcatVcfs.concat_vcf
      File vcf_idx = ConcatVcfs.concat_vcf_idx
    }
  }

task RunDupholdPerContig{
  input{
    String prefix
    String contig
    File bam_or_cram_file
    File bam_or_cram_index
    File vcf_file
    File vcf_index
    File ref_fasta
    File ref_fai
    File ref_dict
    String duphold_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 10, 
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  output {
    File bcf = "~{prefix}.~{contig}.bcf"
  }
  command <<<

    set -Eeuo pipefail
    
    duphold -t 4 \
    -v ~{vcf_file} \
    -b ~{bam_or_cram_file} \
    -f ~{ref_fasta} \
    -o ~{prefix}.~{contig}.bcf

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: duphold_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task Bcf2Vcf{
  input{
    String prefix
    String contig
    File bcf
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 5,
    boot_disk_gb: 5,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output{
    File vcf = "~{prefix}.~{contig}.duphold.vcf.gz"
  }
  command <<<
      set -Eeuo pipefail
      bcftools view ~{bcf} | bgzip > ~{prefix}.~{contig}.duphold.vcf.gz
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

task RunDuphold{
  input{
    String prefix
    File bam_or_cram_file
    File bam_or_cram_index
    File vcf_file
    File ref_fasta
    File ref_fai
    File ref_dict
    String duphold_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 10, 
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  output {
    File bcf = "~{prefix}.bcf"
  }
  command <<<

    set -Eeuo pipefail
    
    duphold -t 4 \
    -v ~{vcf_file} \
    -b ~{bam_or_cram_file} \
    -f ~{ref_fasta} \
    -o ~{prefix}.bcf

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: duphold_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}





