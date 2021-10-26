version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "TasksBenchmark.wdl" as tasks10

import "VaPoR.wdl" as vapor
import "Duphold.wdl" as duphold
import "RdPeSrAnno.wdl" as rdpesr

# WARNING: This workflow is potentially very expensive! Start small and scale gradually, or consider running the
# subworkflows separately.

workflow BenchmarkAnnotation {
  input {
    String prefix
    String il_bam
    String il_bam_bai
    String pb_bam
    String pb_bam_bai
    File vcf_file

    File pe_metrics
    File sr_metrics
    File rd_metrics

    File ref_fasta
    File ref_fai
    File ref_dict
    File contig_list

    String pacbio_benchmark_docker
    String vapor_docker
    String duphold_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_Vapor 
    RuntimeAttr? runtime_attr_duphold
    RuntimeAttr? runtime_attr_rdpesr
    RuntimeAttr? runtime_attr_bcf2vcf
    RuntimeAttr? runtime_attr_LocalizeCram
    RuntimeAttr? runtime_attr_vcf2bed
    RuntimeAttr? runtime_attr_SplitVcf
    RuntimeAttr? runtime_attr_ConcatBeds
    RuntimeAttr? runtime_attr_ConcatVcfs

  }

  Array[String] contigs = transpose(read_tsv(contig_list))[0]
  scatter ( contig in contigs ) {

    call tasks10.LocalizeCram as LocalizeCramPB{
      input:
        contig = contig,
        ref_fasta=ref_fasta,
        ref_fai=ref_fai,
        ref_dict=ref_dict,
        bam_or_cram_file=pb_bam,
        bam_or_cram_index=pb_bam_bai,
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

    call RunDupholdPerContig as RunDupholdPerContigPB{
      input:
        prefix = prefix,
        contig = contig,
        bam_or_cram_file=LocalizeCramPB.local_bam,
        bam_or_cram_index=LocalizeCramPB.local_bai,
        vcf_file = SplitVcf.contig_vcf,
        vcf_index = SplitVcf.contig_vcf_index,
        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        pacbio_benchmark_docker = duphold_docker,
        runtime_attr_override = runtime_attr_duphold
      }

    call Bcf2Vcf as Bcf2VcfPB{
      input:
        prefix = prefix,
        contig = contig,
        bcf = RunDupholdPerContigPB.bcf,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_bcf2vcf
      }
 
    call RunVaPoR{
      input:
        prefix = prefix,
        contig = contig,
        bam_or_cram_file=LocalizeCramPB.local_bam,
        bam_or_cram_index=LocalizeCramPB.local_bai,
        bed = vcf2bed.bed,
        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        pacbio_benchmark_docker = vapor_docker,
        runtime_attr_override = runtime_attr_Vapor
      }

    call RunRdPeSrAnnotation{
      input:
        prefix = prefix,
        contig = contig,
        bam_or_cram_file=LocalizeCramPB.local_bam,
        bam_or_cram_index=LocalizeCramPB.local_bai,
        bed = vcf2bed.bed,
        pe_metrics = pe_metrics,
        sr_metrics = sr_metrics,
        rd_metrics = rd_metrics,
        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        ref_dict=ref_dict,
        pacbio_benchmark_docker = pacbio_benchmark_docker,
        runtime_attr_override = runtime_attr_rdpesr
      }

    call tasks10.LocalizeCramRequestPay as LocalizeCramIL{
      input:
        contig = contig,
        ref_fasta=ref_fasta,
        ref_fai=ref_fai,
        ref_dict=ref_dict,
        project_id="talkowski-sv-gnomad",
        bam_or_cram_file=il_bam,
        bam_or_cram_index=il_bam_bai,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_LocalizeCram
      }

    call RunDupholdPerContig as RunDupholdPerContigIL{
      input:
        prefix = prefix,
        contig = contig,
        bam_or_cram_file=LocalizeCramIL.local_bam,
        bam_or_cram_index=LocalizeCramIL.local_bai,
        vcf_file = SplitVcf.contig_vcf,
        vcf_index = SplitVcf.contig_vcf_index,
        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        pacbio_benchmark_docker = duphold_docker,
        runtime_attr_override = runtime_attr_duphold
      }

    call Bcf2Vcf as Bcf2VcfIL{
      input:
        prefix = prefix,
        contig = contig,
        bcf = RunDupholdPerContigIL.bcf,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_bcf2vcf
      }
    }

  call MiniTasks.ConcatVcfs as ConcatVcfsPB{
    input:
      vcfs=Bcf2VcfPB.vcf,
      outfile_prefix="~{prefix}.PB",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_ConcatVcfs
  }

  call MiniTasks.ConcatVcfs as ConcatVcfsIL{
    input:
      vcfs=Bcf2VcfIL.vcf,
      outfile_prefix="~{prefix}.PB",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_ConcatVcfs
  }

  call MiniTasks.ConcatBeds as ConcatBeds{
    input:
      shard_bed_files=RunVaPoR.vapor,
      prefix=prefix,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_ConcatBeds
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
      File duphold_vcf_pb = ConcatVcfsPB.concat_vcf
      File duphold_vcf_il = ConcatVcfsIL.concat_vcf
      File vapor_bed = ConcatBeds.merged_bed_file
      File PesrAnno = ConcatPesrAnno.merged_bed_file
      File RdAnno = ConcatRdAnno.merged_bed_file
      File RdAnnoLeFlank = ConcatRdAnnoLeFlank.merged_bed_file
      File RdAnnoRiFlank = ConcatRdAnnoRiFlank.merged_bed_file
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
    String pacbio_benchmark_docker
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
    docker: pacbio_benchmark_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task vcf2bed{
  input{
    File vcf
    File vcf_index
    String sv_pipeline_docker
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
  String filename = basename(vcf, ".vcf.gz")

  output {
    File bed = "${filename}.bed"
  }

  command <<<

    set -Eeuo pipefail
    
    svtk vcf2bed -i SVTYPE -i SVLEN ~{vcf} tmp1.bed
    
    cat \
      <(awk '{if ($5=="DEL") print}' tmp1.bed | cut -f1-5)  \
      <(awk '{if ($5=="DUP") print}' tmp1.bed | cut -f1-5) \
      <(awk '{if ($5=="INV") print}' tmp1.bed | cut -f1-5) \
      > ${filename}.bed

    paste -d '_' \
    <(awk '{if ($5=="INS") print}' tmp1.bed | cut -f1-5) \
    <(awk '{if ($5=="INS") print}' tmp1.bed | cut -f6) \
    >> ${filename}.bed

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
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
    String pacbio_benchmark_docker
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
    docker: pacbio_benchmark_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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
    String pacbio_benchmark_docker
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
    File vapor = "~{bed}.vapor"
    File vapor_plot = "~{prefix}.~{contig}.tar.gz"
  }

  command <<<

    set -Eeuo pipefail

    mkdir ~{prefix}.~{contig}
  
    vapor bed \
    --sv-input ~{bed} \
    --output-path ~{prefix}.~{contig} \
    --reference ~{ref_fasta} \
    --pacbio-input ~{bam_or_cram_index} \

    tar -czf ~{prefix}.~{contig}.tar.gz ~{prefix}.~{contig}
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: pacbio_benchmark_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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
    String pacbio_benchmark_docker
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
    Rscript /modify_bed_for_PE_SR_RD_labeling.R -i ~{bed}

    python /add_SR_PE_to_PB_INS.V2.py ~{bed} ~{pe_metrics} ~{sr_metrics} 

    zcat ~{rd_metrics} | grep -v '@' | grep -v CONTIG |bgzip >  bincov.tsv.gz
    Rscript /bincov_to_normCov.R -i bincov.tsv.gz
    bgzip normCov.tsv
    tabix normCov.tsv.gz

    python /add_RD_to_SVs.py ~{bed} normCov.tsv.gz
    python /add_RD_to_SVs.py ~{filename}.ri_flank normCov.tsv.gz
    python /add_RD_to_SVs.py ~{filename}.le_flank normCov.tsv.gz

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: pacbio_benchmark_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}



