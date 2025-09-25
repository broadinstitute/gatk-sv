version 1.0

workflow ProcessLrVsSrStat {
  input {
    File full_vcf
    File tp_vcf
    File vcf2bed_py
    File add_GC_R
    File calcu_stat_R
    File SVID_GC
    Boolean related = false   # default is false

    String sv_base_mini_docker

  }

  call Vcf2Bed { 
      input: 
          vcf = full_vcf, 
          vcf2bed_py = vcf2bed_py,
          docker_file = sv_base_mini_docker

      }
  
  call AddGC { 
      input: 
          bed = Vcf2Bed.bed, 
          add_GC_R = add_GC_R, 
          SVID_GC = SVID_GC,
          docker_file = sv_base_mini_docker
      }

  call CalcuStat { 
      input: 
          bed = AddGC.out_bed, 
          calcu_stat_R = calcu_stat_R, 
          related = related,
          String appdix = "full",
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
          String appdix = "TP",
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
  }

  String prefix = basename(bed, ".bed")

  command <<<
    Rscript ~{add_GC_R} ~{bed} ~{SVID_GC} ~{prefix}.with_GC
  >>>

  output {
    File out_bed = "~{prefix}.with_GC"
  }

  runtime {
    docker: docker_file
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

}







