version 1.0
import "Structs.wdl"

workflow ProcessLrVsSrStat {
  input {
    File full_vcf
    File tp_vcf
    File vcf2bed_py
    File add_GC_R
    File calcu_stat_R
    File SVID_GC
    Array[String] contig_list
    Boolean related = false   # default is false

    String sv_base_mini_docker

  }

  call Vcf2Bed { 
      input: 
          vcf = full_vcf, 
          vcf2bed_py = vcf2bed_py,
          docker_file = sv_base_mini_docker

      }
  
  scatter (contig in contig_list){

    call SplitSvidGc{
        input:
          SVID_GC = SVID_GC,
          contig = contig,
          docker_file = sv_base_mini_docker
    }

    call SplitBed{
        input:
          bed = Vcf2Bed.bed, 
          contig = contig,
          docker_file = sv_base_mini_docker
    }

    call AddGC { 
        input: 
            bed = SplitBed.contig_bed, 
            SVID_GC = SplitSvidGc.contig_SVID_GC,
            add_GC_R = add_GC_R, 
            docker_file = sv_base_mini_docker
        }



  }
  call CalcuStat { 
      input: 
          bed = AddGC.out_bed, 
          calcu_stat_R = calcu_stat_R, 
          related = related,
          appdix = "full",
          docker_file = sv_base_mini_docker
      }

  call Vcf2Bed as TpVcf2Bed { 
      input: 
          vcf = tp_vcf, 
          vcf2bed_py = vcf2bed_py,
          docker_file = sv_base_mini_docker
          }
  
  call CalcuStat as TpCalcuStat { 
      input: 
          bed = TpVcf2Bed.bed, 
          calcu_stat_R = calcu_stat_R, 
          related = related,
          appdix = "TP",
          docker_file = sv_base_mini_docker
     }

  output {
    File full_stat = CalcuStat.stat
    File tp_stat = TpCalcuStat.stat
  }
}

task Vcf2Bed {
  input {
    File vcf
    File vcf2bed_py
    String docker_file
  }

  String prefix = basename(vcf, ".vcf.gz")

  command <<<
    python ~{vcf2bed_py} ~{vcf} ~{prefix}.bed
  >>>

  output {
    File bed = "~{prefix}.bed"
  }

  runtime {
    docker: docker_file
  }
}

task AddGC {
  input {
    File bed
    File add_GC_R
    File SVID_GC
    String docker_file
    RuntimeAttr? runtime_attr_override

  }

  String prefix = basename(bed, ".bed")

  command <<<
    Rscript ~{add_GC_R} ~{bed} ~{SVID_GC} ~{prefix}.with_GC
  >>>

  output {
    File out_bed = "~{prefix}.with_GC"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(bed, "GiB")*2),
    disk_gb: 15 + ceil(size(bed, "GiB")*2),
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_file
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task CalcuStat {
  input {
    File bed
    File calcu_stat_R
    String appdix
    Boolean related = false   # default is false
    String docker_file
  }

  String prefix = basename(bed, ".with_GC")

  command <<<
    if [ "~{related}" == "true" ]; then
      Rscript ~{calcu_stat_R} -i ~{bed} -o ~{prefix}_~{appdix}_stat -r
    else
      Rscript ~{calcu_stat_R} -i ~{bed} -o ~{prefix}_~{appdix}_stat
    fi
  >>>

  output {
    File stat = "~{prefix}_~{appdix}_stat"
  }

  runtime {
    docker: docker_file
  }
}

task SplitSvidGc {
  input {
    File SVID_GC        # input file
    String contig       # chromosome/contig name
    String docker_file
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euo pipefail

    head -1 ~{SVID_GC} > ~{contig}.SVID_GC
    sed -e 's/_/\t/' ~{SVID_GC} \
      | awk -v c="~{contig}" '{ if ($1 == c) print }' \
      | sed -e 's/\t/_/' >> ~{contig}.SVID_GC
  >>>

  output {
    File contig_SVID_GC = "~{contig}.SVID_GC"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3 + ceil(size(SVID_GC, "GiB")*2),
    disk_gb: 5 + ceil(size(SVID_GC, "GiB")*2),
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_file
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task SplitBed {
  input {
    File bed   # input BED file
    String contig
    String docker_file
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euo pipefail

    # Get header (assume first line is header)
    head -1 ~{bed} > ~{contig}.bed

    set -euo pipefail

    awk '{if ($1==~{contig}) print}' ~{bed} >> ~{contig}.bed
  >>>

  output {
    File contig_bed = "~{contig}.bed"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 5 + ceil(size(bed, "GiB")*2),
    disk_gb: 10 + ceil(size(bed, "GiB")*2),
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_file
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

