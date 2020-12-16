version 1.0

import "Structs.wdl"

workflow CleanVcf1 {
  input {
    File vcf
    File background_list
    File ped_file
    String sv_pipeline_docker
    String linux_docker
    File bothsides_pass_list
    File allosome_fai
    RuntimeAttr? runtime_attr_override # TODO
  }

  call CreateEmptyFile {
    input:
      linux_docker=linux_docker
  }

  call CleanVcf1_1 {
    input:
      vcf=vcf,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf1_2 {
    input:
      EV_update_vcf=CleanVcf1_1.EV_update_vcf,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf1_3 {
    input:
      EV_update_vcf=CleanVcf1_1.EV_update_vcf,
      vcf_convert_svtype=CleanVcf1_2.vcf_convert_svtype,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf1_4 {
    input:
      convertsvtype_vcf=CleanVcf1_3.convertsvtype_vcf,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf1_5 {
    input:
      convertsvtype_vcf=CleanVcf1_3.convertsvtype_vcf,
      vargq_persample=CleanVcf1_4.vargq_persample,
      sv_pipeline_docker=sv_pipeline_docker
  }

  if (CleanVcf1_5.count_xy > 0) {
    call CleanVcf1_6 {
      input:
        cleaninfo_vcf=CleanVcf1_5.cleaninfo_vcf,
        cleaninfo_vcf_tbi=CleanVcf1_5.cleaninfo_vcf_tbi,
        ped_file=ped_file,
        sv_pipeline_docker=sv_pipeline_docker
    }

    if (CleanVcf1_6.clean_bed_ids_count > 0) {
      call CleanVcf1_7 {
        input:
          allosome_fai=allosome_fai,
          cleaninfo_vcf=CleanVcf1_5.cleaninfo_vcf,
          cleaninfo_vcf_tbi=CleanVcf1_5.cleaninfo_vcf_tbi,
          clean_bed_ids=CleanVcf1_6.clean_bed_ids,
          male=CleanVcf1_6.male,
          sv_pipeline_docker=sv_pipeline_docker
      }
      call CleanVcf1_8 {
        input:
          allosome_fai=allosome_fai,
          cleaninfo_vcf=CleanVcf1_5.cleaninfo_vcf,
          cleaninfo_vcf_tbi=CleanVcf1_5.cleaninfo_vcf_tbi,
          clean_bed_ids=CleanVcf1_6.clean_bed_ids,
          female=CleanVcf1_6.female,
          sv_pipeline_docker=sv_pipeline_docker
      }
      call CleanVcf1_9 {
        input:
          RD_CN_sexcheck_FORMAT_male=CleanVcf1_7.RD_CN_sexcheck_FORMAT_male,
          sv_pipeline_docker=sv_pipeline_docker
      }
      call CleanVcf1_10 {
        input:
          RD_CN_sexcheck_FORMAT_female=CleanVcf1_8.RD_CN_sexcheck_FORMAT_female,
          sv_pipeline_docker=sv_pipeline_docker
      }
    }

    call CleanVcf1_11 {
      input:
        cleaninfo_vcf=CleanVcf1_5.cleaninfo_vcf,
        cleaninfo_vcf_tbi=CleanVcf1_5.cleaninfo_vcf_tbi,
        clean_bed=CleanVcf1_6.clean_bed,
        male_median_value_pervar=select_first([CleanVcf1_9.male_median_value_pervar, CreateEmptyFile.empty]),
        female_median_value_pervar=select_first([CleanVcf1_10.female_median_value_pervar, CreateEmptyFile.empty]),
        sv_pipeline_docker=sv_pipeline_docker
    }

    call CleanVcf1_12 {
      input:
        sexchr_revise_1=CleanVcf1_11.sexchr_revise_1,
        cleaninfo_vcf=CleanVcf1_5.cleaninfo_vcf,
        cleaninfo_vcf_tbi=CleanVcf1_5.cleaninfo_vcf_tbi,
        clean_bed=CleanVcf1_6.clean_bed,
        male_median_value_pervar=select_first([CleanVcf1_9.male_median_value_pervar, CreateEmptyFile.empty]),
        female_median_value_pervar=select_first([CleanVcf1_10.female_median_value_pervar, CreateEmptyFile.empty]),
        sv_pipeline_docker=sv_pipeline_docker
    }

    call CleanVcf1_13 {
      input:
        cleaninfo_vcf=CleanVcf1_5.cleaninfo_vcf,
        cleaninfo_vcf_csi=CleanVcf1_5.cleaninfo_vcf_csi,
        male=CleanVcf1_6.male,
        sv_pipeline_docker=sv_pipeline_docker
    }

    call CleanVcf1_14 {
      input:
        cleaninfo_vcf=CleanVcf1_5.cleaninfo_vcf,
        cleaninfo_vcf_csi=CleanVcf1_5.cleaninfo_vcf_csi,
        female=CleanVcf1_6.female,
        sv_pipeline_docker=sv_pipeline_docker
    }

    call CleanVcf1_15 {
      input:
        male_vcf=CleanVcf1_13.male_vcf,
        sexchr_revise_2=CleanVcf1_12.sexchr_revise_2,
        sv_pipeline_docker=sv_pipeline_docker
    }

    call CleanVcf1_16 {
      input:
        male_vcf=CleanVcf1_13.male_vcf,
        sexchr_revise_2=CleanVcf1_12.sexchr_revise_2,
        sv_pipeline_docker=sv_pipeline_docker
    }

    if ((CleanVcf1_15.count + CleanVcf1_16.count) > 0) {
      call CleanVcf1_17 {
        input:
          male_vcf=CleanVcf1_13.male_vcf,
          male_dup_revise_txt=CleanVcf1_16.male_dup_revise_txt,
          male_del_revise_txt=CleanVcf1_15.male_del_revise_txt,
          sv_pipeline_docker=sv_pipeline_docker
      }
    }

    if (CleanVcf1_5.count_y > 0) {
      call CleanVcf1_18 {
        input:
          female_vcf=CleanVcf1_14.female_vcf,
          sv_pipeline_docker=sv_pipeline_docker
      }
      call CleanVcf1_19 {
        input:
          female_vcf=CleanVcf1_14.female_vcf,
          female_y_revise_txt=CleanVcf1_18.female_y_revise_txt,
          sv_pipeline_docker=sv_pipeline_docker
      }
    }

    if (CleanVcf1_6.ped_file_count > 0) {
      call CleanVcf1_20 {
        input:
          cleaninfo_vcf=CleanVcf1_5.cleaninfo_vcf,
          cleaninfo_vcf_csi=CleanVcf1_5.cleaninfo_vcf_csi,
          ped_file=ped_file,
          sv_pipeline_docker=sv_pipeline_docker
      }
      call CleanVcf1_21 {
        input:
          other_vcf=CleanVcf1_20.other_vcf,
          sv_pipeline_docker=sv_pipeline_docker
      }
      call CleanVcf1_22 {
        input:
          other_vcf=CleanVcf1_20.other_vcf,
          other_revise_txt=CleanVcf1_21.other_revise_txt,
          sv_pipeline_docker=sv_pipeline_docker
      }
      call CleanVcf1_23 {
        input:
          cleanmale_vcf=select_first([CleanVcf1_17.cleanmale_vcf, CleanVcf1_13.male_vcf]),
          cleanfemale_vcf=select_first([CleanVcf1_19.cleanfemale_vcf, CleanVcf1_14.female_vcf]),
          cleanother_vcf=CleanVcf1_22.cleanother_vcf,
          sv_pipeline_docker=sv_pipeline_docker
      }
    }
    if (CleanVcf1_6.ped_file_count == 0) {
      call CleanVcf1_24 {
        input:
          cleanmale_vcf=select_first([CleanVcf1_17.cleanmale_vcf, CleanVcf1_13.male_vcf]),
          cleanfemale_vcf=select_first([CleanVcf1_19.cleanfemale_vcf, CleanVcf1_14.female_vcf]),
          sv_pipeline_docker=sv_pipeline_docker
      }
    }

    call CleanVcf1_25 {
      input:
        combinedsex_vcf=select_first([CleanVcf1_23.combinedsex_vcf, CleanVcf1_24.combinedsex_vcf]),
        combinedsex_vcf_tbi=select_first([CleanVcf1_23.combinedsex_vcf_tbi, CleanVcf1_24.combinedsex_vcf_tbi]),
        cleaninfo_vcf=CleanVcf1_5.cleaninfo_vcf,
        cleaninfo_vcf_tbi=CleanVcf1_5.cleaninfo_vcf_tbi,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  call CleanVcf1_26 {
    input:
      background_list=background_list,
      cleanallo_vcf=select_first([CleanVcf1_25.cleanallo_vcf, CleanVcf1_5.cleaninfo_vcf]),
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf1_27 {
    input:
      int_vcf=CleanVcf1_26.int_vcf,
      bothsides_pass_list=bothsides_pass_list,
      sv_pipeline_docker=sv_pipeline_docker
  }

  output {
    File include_list=CleanVcf1_1.include_list
    File sex=select_first([CleanVcf1_12.sexchr_revise_2, CreateEmptyFile.empty])
    File intermediate_vcf=CleanVcf1_27.intermediate_vcf
    File intermediate_vcf_idx=CleanVcf1_27.intermediate_vcf_idx
  }
}


task CleanVcf1_1 {
  input {
    File vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    ##get sampleids from VCF##
    zcat ~{vcf} \
      |sed -n '1,1000p' \
      |egrep "^#" \
      |tail -n -1 \
      |cut -f10- \
      |tr '\t' '\n' \
      > includelist.txt

    ##convert EV integer back into string##
    /opt/sv-pipeline/04_variant_resolution/scripts/replace_ev_numeric_code_with_string.py ~{vcf} - | bgzip -c > EV.update.vcf.gz
  >>>

  output {
    File include_list="includelist.txt"
    File EV_update_vcf="EV.update.vcf.gz"
  }
}

task CleanVcf1_2 {
  input {
    File EV_update_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(EV_update_vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    ##convert all alt to svtype and alt to N##
    svtk vcf2bed ~{EV_update_vcf} stdout -i SVTYPE  \
      |awk -F"\t" '{ if ($5!~"ME")$5=$7; print $4"\t" "<"$5 ">"}' \
      |gzip \
      >vcf.convert.svtype.bed.gz
  >>>

  output {
    File vcf_convert_svtype="vcf.convert.svtype.bed.gz"
  }
}
task CleanVcf1_3 {
  input {
    File EV_update_vcf
    File vcf_convert_svtype
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([EV_update_vcf, vcf_convert_svtype], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    zcat ~{EV_update_vcf} \
      |awk -F'\t' -v OFS='\t' 'ARGIND==1{inFileA[$1]=$2; next} {if ($3 in inFileA && $1!~"#") $5=inFileA[$3]; print }'  \
        <(zcat ~{vcf_convert_svtype}) - \
      |awk '{if ($1!~"#") $4="N"; print}' OFS='\t' \
      |bgzip \
      >convertsvtype.vcf.gz
  >>>

  output {
    File convertsvtype_vcf="convertsvtype.vcf.gz"
  }
}

task CleanVcf1_4 {
  input {
    File convertsvtype_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(convertsvtype_vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail
    ##get rid of multiallelic tage in INFO field and add varGQ to QUAL column and Members field##
    svtk vcf2bed ~{convertsvtype_vcf} stdout -i varGQ \
      |awk -F"\t" '{print $4 "\t" $7}' \
      >vargq.persample
  >>>

  output {
    File vargq_persample="vargq.persample"
  }
}

task CleanVcf1_5 {
  input {
    File convertsvtype_vcf
    File vargq_persample
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([convertsvtype_vcf, vargq_persample], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    # use BCFTOOLS 1.9, installed in /usr/local/bin/bcftools in our docker
    BCFTOOLS=/usr/local/bin/bcftools

    zcat ~{convertsvtype_vcf} \
      |sed 's/;MULTIALLELIC//g' \
      |sed 's/UNRESOLVED;//g' \
      |sed 's/;varGQ=[0-9]*//g' \
      |awk -F'\t' -v OFS='\t' 'ARGIND==1{inFileA[$1]=$2; next} {if ($3 in inFileA && $1!~"#") $6=inFileA[$3]; print }' ~{vargq_persample} - \
      |bgzip \
      >cleaninfo.vcf.gz

    tabix -p vcf cleaninfo.vcf.gz
    ${BCFTOOLS} index cleaninfo.vcf.gz

    zcat cleaninfo.vcf.gz|awk '{if (($1~"X" || $1~"Y") && $1!~"#" ) print}'|wc -l > count_xy.txt
    zcat cleaninfo.vcf.gz|awk '{if ($1~"Y" && $1!~"#") print}'|wc -l > count_y.txt
  >>>

  output {
    File cleaninfo_vcf="cleaninfo.vcf.gz"
    File cleaninfo_vcf_tbi="cleaninfo.vcf.gz.tbi"
    File cleaninfo_vcf_csi="cleaninfo.vcf.gz.csi"
    Int count_xy = read_int("count_xy.txt")
    Int count_y= read_int("count_y.txt")
  }
}

task CleanVcf1_6 {
  input {
    File cleaninfo_vcf
    File cleaninfo_vcf_tbi
    File ped_file
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([cleaninfo_vcf, ped_file], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail
    svtk vcf2bed ~{cleaninfo_vcf} stdout \
      |awk '{if (($5=="DEL" || $5=="DUP") && $3-$2>=5000 && ($1~"X" || $1~"Y") && $1!~"#") print}' \
      >clean.bed

    awk '{print $4}' clean.bed>clean.bed.ids.txt

    ##male##
    awk '{if ($5==1) print $2}' ~{ped_file} \
      |fgrep -wf <(zcat ~{cleaninfo_vcf}|head -n 1000|fgrep "CHROM"|fgrep POS|cut -f10-|tr '\t' '\n') >male.txt

    ##female##
    awk '{if ($5==2) print $2}' ~{ped_file} \
      |fgrep -wf <(zcat ~{cleaninfo_vcf}|head -n 1000|fgrep "CHROM"|fgrep POS|cut -f10-|tr '\t' '\n') >female.txt

    cat clean.bed.ids.txt|wc -l > clean_bed_ids_count.txt
    awk '{if ($5!=2 && $5!=1) print $2}' ~{ped_file}|wc -l > ped_file_count.txt
  >>>

  output {
    File clean_bed="clean.bed"
    File clean_bed_ids="clean.bed.ids.txt"
    File male="male.txt"
    File female="female.txt"
    Int clean_bed_ids_count=read_int("clean_bed_ids_count.txt")
    Int ped_file_count=read_int("ped_file_count.txt")
  }
}

task CleanVcf1_7 {
  input {
    File allosome_fai
    File cleaninfo_vcf
    File cleaninfo_vcf_tbi
    File clean_bed_ids
    File male
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([cleaninfo_vcf, clean_bed_ids], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    # use BCFTOOLS 1.9, installed in /usr/local/bin/bcftools in our docker
    BCFTOOLS=/usr/local/bin/bcftools

    awk '{print $1"\t0\t"$2}' < ~{allosome_fai} > allosomes.list
    ${BCFTOOLS} query -R allosomes.list -S ~{male} -i 'ID=@~{clean_bed_ids}' -f '[%ID\t%SAMPLE\t%RD_CN\n]' ~{cleaninfo_vcf} \
      | awk '{if ($3!=".") print}' \
      | gzip > RD_CN.sexcheck.FORMAT.male.gz
  >>>

  output {
    File RD_CN_sexcheck_FORMAT_male="RD_CN.sexcheck.FORMAT.male.gz"
  }
}

task CleanVcf1_8 {
  input {
    File cleaninfo_vcf
    File cleaninfo_vcf_tbi
    File allosome_fai
    File clean_bed_ids
    File female
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(clean_bed_ids, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    # use BCFTOOLS 1.9, installed in /usr/local/bin/bcftools in our docker
    BCFTOOLS=/usr/local/bin/bcftools

    awk '{print $1"\t0\t"$2}' < ~{allosome_fai} > allosomes.list
    ${BCFTOOLS} query -R allosomes.list -S ~{female} -i 'ID=@~{clean_bed_ids}' -f '[%ID\t%SAMPLE\t%RD_CN\n]' ~{cleaninfo_vcf} \
      | awk '{if ($3!=".") print}' \
      | gzip > RD_CN.sexcheck.FORMAT.female.gz
  >>>

  output {
    File RD_CN_sexcheck_FORMAT_female="RD_CN.sexcheck.FORMAT.female.gz"
  }
}

task CleanVcf1_9 {
  input {
    File RD_CN_sexcheck_FORMAT_male
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(RD_CN_sexcheck_FORMAT_male, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    zcat ~{RD_CN_sexcheck_FORMAT_male}| Rscript -e 'd<-read.table("stdin")' \
      -e 'x<-tapply(d[,3],d[,1],median)' \
      -e 'write.table(x,"male.median.value.pervar.txt",col.names=FALSE,quote=FALSE,sep = "\t")'
  >>>

  output {
    File male_median_value_pervar="male.median.value.pervar.txt"
  }
}

task CleanVcf1_10 {
  input {
    File RD_CN_sexcheck_FORMAT_female
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(RD_CN_sexcheck_FORMAT_female, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    zcat ~{RD_CN_sexcheck_FORMAT_female}| Rscript -e 'd<-read.table("stdin")' \
      -e 'x<-tapply(d[,3],d[,1],median)' \
      -e 'write.table(x,"female.median.value.pervar.txt",col.names=FALSE,quote=FALSE,sep = "\t")'
  >>>

  output {
    File female_median_value_pervar="female.median.value.pervar.txt"
  }
}

task CreateEmptyFile {
  input {
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr runtime_default = object {
                                  mem_gb: 1.0,
                                  disk_gb: 10,
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: linux_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    touch empty.txt
  >>>

  output {
    File empty="empty.txt"
  }
}

task CleanVcf1_11 {
  input {
    File cleaninfo_vcf
    File cleaninfo_vcf_tbi
    File clean_bed
    File? male_median_value_pervar
    File? female_median_value_pervar
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([cleaninfo_vcf, clean_bed, male_median_value_pervar, female_median_value_pervar], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail
    ##Pull out ids where male copy state 1 to normal when female normal and on X##
    echo "">sexchr.revise.1.txt

    if [ $(awk '{if (($5=="DEL" || $5=="DUP") && $3-$2>=5000) print }' ~{clean_bed}|awk '{if (($1~"X") && $1!~"#" ) print}'|wc -l) -gt 0 ]
    then
      awk '{if ($2==1) print $1}' ~{male_median_value_pervar} \
        |{ fgrep -wf <(awk '{if ($2==2) print $1}' ~{female_median_value_pervar}) || true; } \
        |{ fgrep -wf  - <(zcat ~{cleaninfo_vcf}|awk '{if ($1~"X" && $1!~"#") print $3}') || true; } \
        >sexchr.revise.1.txt
    fi
  >>>

  output {
    File sexchr_revise_1="sexchr.revise.1.txt"
  }
}

task CleanVcf1_12 {
  input {
    File sexchr_revise_1
    File cleaninfo_vcf
    File cleaninfo_vcf_tbi
    File clean_bed
    File male_median_value_pervar
    File female_median_value_pervar
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([sexchr_revise_1, cleaninfo_vcf, clean_bed, male_median_value_pervar, female_median_value_pervar], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail
    mv ~{sexchr_revise_1} sexchr.revise.2.txt
    if [ $(awk '{if (($5=="DEL" || $5=="DUP") && $3-$2>=5000) print }' ~{clean_bed}|awk '{if (($1~"Y") && $1!~"#" ) print}'|wc -l) -gt 0 ]
    then
      awk '{if ($2==1) print $1}' ~{male_median_value_pervar} \
        |{ fgrep -wf <(awk '{if ($2==0) print $1}' ~{female_median_value_pervar}) || true; } \
        |{ fgrep -wf - <(zcat ~{cleaninfo_vcf}|awk '{if ($1~"Y" && $1!~"#") print $3}') || true; } \
        >>sexchr.revise.2.txt
    fi
  >>>

  output {
    File sexchr_revise_2="sexchr.revise.2.txt"
  }
}


task CleanVcf1_13 {
  input {
    File cleaninfo_vcf
    File cleaninfo_vcf_csi
    File male
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([cleaninfo_vcf, male], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    # use BCFTOOLS 1.9, installed in /usr/local/bin/bcftools in our docker
    BCFTOOLS=/usr/local/bin/bcftools

    ##Pull out male sex chr##
    ${BCFTOOLS} view ~{cleaninfo_vcf} -S ~{male} -r chrX:1-1000000000,chrY:1-1000000000,X:1-1000000000,Y:1-1000000000 --no-update|bgzip>male.vcf.gz
    ${BCFTOOLS} index male.vcf.gz
  >>>

  output {
    File male_vcf="male.vcf.gz"
    File male_vcf_csi="male.vcf.gz.csi"
  }
}


task CleanVcf1_14 {
  input {
    File cleaninfo_vcf
    File cleaninfo_vcf_csi
    File female
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([cleaninfo_vcf, female], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    # use BCFTOOLS 1.9, installed in /usr/local/bin/bcftools in our docker
    BCFTOOLS=/usr/local/bin/bcftools

    ##Pull out female sex chr##
    ${BCFTOOLS} view ~{cleaninfo_vcf} -S ~{female} -r chrX:1-1000000000,chrY:1-1000000000,X:1-1000000000,Y:1-1000000000 --no-update|bgzip>female.vcf.gz
    ${BCFTOOLS} index female.vcf.gz
  >>>

  output {
    File female_vcf="female.vcf.gz"
    File female_vcf_csi="female.vcf.gz.csi"
  }
}


task CleanVcf1_15 {
  input {
    File male_vcf
    File sexchr_revise_2
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([male_vcf, sexchr_revise_2], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    zcat ~{male_vcf}\
      |awk -F'\t' '{if ($5~"DEL" && $1!~"#") print $0 "\t" "ENDOFLINE"}' \
      |{ fgrep -wf ~{sexchr_revise_2} || true; } \
      |tr '\t' '\n' \
      |awk -F':' '{if ($3>=1 && NF>4 && $1!="GT") $1="0/0";else if ($3==0 && NF>4 && $1!="GT" ) $1="0/1"; if (NF>4 && $1!="GT") $3=$3+1;print}' OFS=":" \
      |tr '\n' '\t' \
      |sed 's/ENDOFLINE/\n/g' \
      |sed -e 's/^[ \t]*//' \
      |sed -e 's/[\t]$//g' \
      |bgzip \
      >male_del.revise.txt.gz
    zcat male_del.revise.txt.gz|wc -l > count.txt
  >>>

  output {
    File male_del_revise_txt="male_del.revise.txt.gz"
    Int count=read_int("count.txt")
  }
}


task CleanVcf1_16 {
  input {
    File male_vcf
    File sexchr_revise_2
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([male_vcf, sexchr_revise_2], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail
    zcat ~{male_vcf}\
      |awk -F'\t' '{if ($5~"DUP" && $1!~"#") print $0 "\t" "ENDOFLINE"}' \
      |{ fgrep -wf ~{sexchr_revise_2} || true; } \
      |tr '\t' '\n' \
      |awk -F':' '{if ($3<=1 && NF>4 && $1!="GT") $1="0/0";else if ($3==2 && NF>4 && $1!="GT" ) $1="0/1";else if (NF>4 && $1!="GT" ) $1="1/1"; if (NF>4 && $1!="GT" ) $3=$3+1;print}' OFS=":" \
      |tr '\n' '\t' \
      |sed 's/ENDOFLINE/\n/g' \
      |sed -e 's/^[ \t]*//' \
      |sed -e 's/[\t]$//g' \
      |bgzip \
      >male_dup.revise.txt.gz
    zcat male_dup.revise.txt.gz|wc -l > count.txt
  >>>

  output {
    File male_dup_revise_txt="male_dup.revise.txt.gz"
    Int count=read_int("count.txt")
  }
}

task CleanVcf1_17 {
  input {
    File male_vcf
    File male_dup_revise_txt
    File male_del_revise_txt
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([male_vcf, male_dup_revise_txt, male_del_revise_txt], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 10.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    # use BCFTOOLS 1.9, installed in /usr/local/bin/bcftools in our docker
    BCFTOOLS=/usr/local/bin/bcftools

    cat <(zcat ~{male_vcf}|fgrep -wvf <(zcat ~{male_dup_revise_txt} ~{male_del_revise_txt}|awk '{print $3}' )) \
      <(zcat ~{male_del_revise_txt} ~{male_dup_revise_txt}|awk '{if ($1!="") print}'|tr ' ' '\t') \
      |vcf-sort \
      |bgzip \
      >cleanmale.vcf.gz
    ${BCFTOOLS} index cleanmale.vcf.gz
  >>>

  output {
    File cleanmale_vcf="cleanmale.vcf.gz"
    File cleanmale_vcf_csi="cleanmale.vcf.gz.csi"
  }
}



task CleanVcf1_18 {
  input {
    File female_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(female_vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail
    zcat ~{female_vcf}\
      |awk -F'\t' '{if ($1!~"#" && $1~"Y") print $0 "\t" "ENDOFLINE"}' \
      |tr '\t' '\n' \
      |awk -F':' '{ if (NF>4 && $1!="GT" ) $1="./."  \
      ;if (NF>4 && $1!="GT" ) $2=$3=$4=$5=$6=$7=$8=$9=".";print}' OFS=":" \
      |tr '\n' '\t' \
      |sed 's/ENDOFLINE/\n/g' \
      |sed -e 's/^[ \t]*//' \
      |sed -e 's/[\t]$//g' \
      |bgzip \
      >female.y.revise.txt.gz
  >>>

  output {
    File female_y_revise_txt="female.y.revise.txt.gz"
  }
}



task CleanVcf1_19 {
  input {
    File female_vcf
    File female_y_revise_txt
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([female_vcf, female_y_revise_txt], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    # use BCFTOOLS 1.9, installed in /usr/local/bin/bcftools in our docker
    BCFTOOLS=/usr/local/bin/bcftools

    cat <(zcat ~{female_vcf} \
      |fgrep -wvf <(zcat ~{female_y_revise_txt}|awk '{print $3}' )) \
      <(zcat ~{female_y_revise_txt}) \
      |vcf-sort \
      |bgzip \
      >cleanfemale.vcf.gz
    ${BCFTOOLS} index cleanfemale.vcf.gz
  >>>

  output {
    File cleanfemale_vcf="cleanfemale.vcf.gz"
    File cleanfemale_vcf_csi="cleanfemale.vcf.gz.csi"
  }
}



task CleanVcf1_20 {
  input {
    File cleaninfo_vcf
    File cleaninfo_vcf_csi
    File ped_file
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([cleaninfo_vcf, ped_file], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    # use BCFTOOLS 1.9, installed in /usr/local/bin/bcftools in our docker
    BCFTOOLS=/usr/local/bin/bcftools

    awk '{if ($5!=2 && $5!=1) print $2}' ~{ped_file}>other.txt
    ${BCFTOOLS} view ~{cleaninfo_vcf} -S other.txt -r chrX:1-1000000000,chrY:1-1000000000,X:1-1000000000,Y:1-1000000000 --no-update|bgzip>other.vcf.gz
    ${BCFTOOLS} index other.vcf.gz
  >>>

  output {
    File other_vcf="other.vcf.gz"
    File other_vcf_csi="other.vcf.gz.csi"
  }
}

task CleanVcf1_21 {
  input {
    File other_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(other_vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    zcat ~{other_vcf}\
      |awk -F'\t' '{if ($1!~"#") print $0 "\t" "ENDOFLINE"}' \
      |tr '\t' '\n' \
      |awk -F':' '{ if (NF>4 && $1!="GT" ) $1="./.";print}' OFS=":" \
      |tr '\n' '\t' \
      |sed 's/ENDOFLINE/\n/g' \
      |sed -e 's/^[ \t]*//' \
      |sed -e 's/[\t]$//g' \
      |bgzip \
      >other.revise.txt.gz
  >>>

  output {
    File other_revise_txt="other.revise.txt.gz"
  }
}



task CleanVcf1_22 {
  input {
    File other_vcf
    File other_revise_txt
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([other_vcf, other_revise_txt], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    # use BCFTOOLS 1.9, installed in /usr/local/bin/bcftools in our docker
    BCFTOOLS=/usr/local/bin/bcftools

    cat <(zcat ~{other_vcf} \
      |fgrep -wvf <(zcat ~{other_revise_txt}|awk '{print $3}' )) \
      <(zcat ~{other_revise_txt}) \
      |vcf-sort \
      |bgzip \
      >cleanother.vcf.gz
    ${BCFTOOLS} index cleanother.vcf.gz
  >>>

  output {
    File cleanother_vcf="cleanother.vcf.gz"
    File cleanother_vcf_csi="cleanother.vcf.gz.csi"
  }
}



task CleanVcf1_23 {
  input {
    File cleanmale_vcf
    File cleanfemale_vcf
    File cleanother_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([cleanmale_vcf, cleanfemale_vcf, cleanother_vcf], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail
    cat <(zcat ~{cleanmale_vcf}|egrep "##") \
      <(paste <(zcat ~{cleanmale_vcf}|egrep -v "##") <(zcat ~{cleanfemale_vcf}|cut -f10-|egrep -v "##") <(zcat ~{cleanother_vcf}|cut -f10-|egrep -v "##") ) \å
      |bgzip \
      >combinedsex.vcf.gz
    tabix -p vcf combinedsex.vcf.gz
  >>>

  output {
    File combinedsex_vcf="combinedsex.vcf.gz"
    File combinedsex_vcf_tbi="combinedsex.vcf.gz.tbi"
  }
}

task CleanVcf1_24 {
  input {
    File cleanmale_vcf
    File cleanfemale_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([cleanmale_vcf, cleanfemale_vcf], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail
    cat <(zcat ~{cleanmale_vcf}|egrep "##") \
      <(paste <(zcat ~{cleanmale_vcf}|egrep -v "##") <(zcat ~{cleanfemale_vcf}|cut -f10-|egrep -v "##"))  \
      |bgzip \
      >combinedsex.vcf.gz
    tabix -p vcf combinedsex.vcf.gz
  >>>

  output {
    File combinedsex_vcf="combinedsex.vcf.gz"
    File combinedsex_vcf_tbi="combinedsex.vcf.gz.tbi"
  }
}

task CleanVcf1_25 {
  input {
    File combinedsex_vcf
    File combinedsex_vcf_tbi
    File cleaninfo_vcf
    File cleaninfo_vcf_tbi
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([combinedsex_vcf, cleaninfo_vcf], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    zcat ~{combinedsex_vcf}|awk '{if ($1!~"#") print $3}'>modified.ids.txt

    ##shuffle sex ids backinto place to match original vcf and back to initial vcf##
    vcf-shuffle-cols -t ~{cleaninfo_vcf} ~{combinedsex_vcf} \
      |awk '{if ($1!~"#") print}' \
      |cat <(zcat ~{cleaninfo_vcf}|fgrep -wvf modified.ids.txt ) - \
      |vcf-sort \
      |bgzip \
      >cleanallo.vcf.gz
  >>>

  output {
    File cleanallo_vcf="cleanallo.vcf.gz"
  }
}

task CleanVcf1_26 {
  input {
    File background_list
    File cleanallo_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([background_list, cleanallo_vcf], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    # the code below will not print any lines if the background list file is empty, so add a dummy sentinel record at the end
    cat ~{background_list} <(echo "XXX_SENTINEL_XXX") > background_list_with_sentinel.list

    ##change tag for SR background failures and Unresolved##
    zcat ~{cleanallo_vcf} \
      |awk -F'\t' -v OFS='\t' 'ARGIND==1{inFileA[$1]; next} {if ($3 in inFileA && $1!~"#") $7=$7";HIGH_SR_BACKGROUND"; print }' <(awk '{print $NF}' background_list_with_sentinel.list) - \
      |awk '{if (NR==2) print $0 "\n" "##FILTER=<ID=HIGH_SR_BACKGROUND,Description=\"High number of SR splits in background samples indicating messy region\">" ;else print}' \
      |awk '{if ($8~"UNRESOLVED") $7=$7";UNRESOLVED";print}' OFS='\t' \
      |awk '{if (NR==2) print $0 "\n" "##FILTER=<ID=UNRESOLVED,Description=\"Variant is unresolved\">" ;else print}' \
      |bgzip \
      >int.vcf.gz
  >>>

  output {
    File int_vcf="int.vcf.gz"
  }
}

task CleanVcf1_27 {
  input {
    File int_vcf
    File bothsides_pass_list
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([int_vcf, bothsides_pass_list], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail
    /opt/sv-pipeline/04_variant_resolution/scripts/add_bothsides_support_filter.py \
      --bgzip \
      --outfile int.w_bothsides.vcf.gz \
      ~{int_vcf} \
      ~{bothsides_pass_list}
    tabix int.w_bothsides.vcf.gz
  >>>

  output {
    File intermediate_vcf="int.w_bothsides.vcf.gz"
    File intermediate_vcf_idx="int.w_bothsides.vcf.gz.tbi"
  }
}
