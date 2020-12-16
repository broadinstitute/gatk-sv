version 1.0

import "Structs.wdl"

workflow CleanVcf5 {
  input {
    File revise_vcf_lines
    File normal_revise_vcf
    File ped_file
    File sex_chr_revise
    File multi_ids
    File? outlier_samples_list
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override # TODO
  }

  call CleanVcf5_1 {
    input:
      normal_revise_vcf=normal_revise_vcf,
      revise_vcf_lines=revise_vcf_lines,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_2 {
    input:
      overlap_revise_vcf=CleanVcf5_1.overlap_revise_vcf,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_3 {
    input:
      overlap_revise_vcf=CleanVcf5_1.overlap_revise_vcf,
      outlier_samples_list=outlier_samples_list,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_4 {
    input:
      overlap_revise_vcf=CleanVcf5_1.overlap_revise_vcf,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_5 {
    input:
      copystate_rd_cn_format=CleanVcf5_3.copystate_rd_cn_format,
      overlap_revise_bed=CleanVcf5_2.overlap_revise_bed,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_6 {
    input:
      copystate_rd_cn_format=CleanVcf5_3.copystate_rd_cn_format,
      overlap_revise_bed=CleanVcf5_2.overlap_revise_bed,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_7 {
    input:
      copystate_rd_cn_format=CleanVcf5_3.copystate_rd_cn_format,
      overlap_revise_bed=CleanVcf5_2.overlap_revise_bed,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_8 {
    input:
      copystate_rd_cn_format=CleanVcf5_3.copystate_rd_cn_format,
      overlap_revise_bed=CleanVcf5_2.overlap_revise_bed,
      gt4copystate=CleanVcf5_7.gt4copystate,
      multi_dup_ids_1=CleanVcf5_6.multi_dup_ids_1,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_9 {
    input:
      overlap_revise_vcf=CleanVcf5_1.overlap_revise_vcf,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_10 {
    input:
      genotype_gt_format=CleanVcf5_4.genotype_gt_format,
      multi_dup_ids=CleanVcf5_8.multi_dup_ids,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_11 {
    input:
      genotype_gt_format=CleanVcf5_4.genotype_gt_format,
      multi_del_ids=CleanVcf5_5.multi_del_ids,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_12 {
    input:
      multi_dup_ids=CleanVcf5_8.multi_dup_ids,
      regeno_bed=CleanVcf5_9.regeno_bed,
      gt5kb_dup_ids_1=CleanVcf5_10.gt5kb_dup_ids_1,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_13 {
    input:
      multi_del_ids=CleanVcf5_5.multi_del_ids,
      gt5kb_del_ids_1=CleanVcf5_11.gt5kb_del_ids_1,
      regeno_bed=CleanVcf5_9.regeno_bed,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_14 {
    input:
      overlap_revise_vcf=CleanVcf5_1.overlap_revise_vcf,
      gt5kb_dup_ids=CleanVcf5_12.gt5kb_dup_ids,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_15 {
    input:
      overlap_revise_vcf=CleanVcf5_1.overlap_revise_vcf,
      gt5kb_del_ids=CleanVcf5_13.gt5kb_del_ids,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_16 {
    input:
      del_int=CleanVcf5_15.del_int,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_17 {
    input:
      dup_int=CleanVcf5_14.dup_int,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_18 {
    input:
      overlap_revise_vcf=CleanVcf5_1.overlap_revise_vcf,
      gt5kb_dup_ids=CleanVcf5_12.gt5kb_dup_ids,
      gt5kb_del_ids=CleanVcf5_13.gt5kb_del_ids,
      dup_revise=CleanVcf5_17.dup_revise,
      del_revise=CleanVcf5_16.del_revise,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_19 {
    input:
      multi_geno_ids_txt=multi_ids,
      multi_del_ids=CleanVcf5_5.multi_del_ids,
      multi_dup_ids=CleanVcf5_8.multi_dup_ids,
      newdepth_geno_vcf=CleanVcf5_18.newdepth_geno_vcf,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_20 {
    input:
      multi_dup_ids=CleanVcf5_8.multi_dup_ids,
      multitagged_vcf=CleanVcf5_19.multitagged_vcf,
      multitagged_vcf_tbi=CleanVcf5_19.multitagged_vcf_tbi,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_21 {
    input:
      multi_del_ids=CleanVcf5_5.multi_del_ids,
      multitagged_vcf=CleanVcf5_19.multitagged_vcf,
      multitagged_vcf_tbi=CleanVcf5_19.multitagged_vcf_tbi,
      dup_multi_revise_vcf=CleanVcf5_20.dup_multi_revise_vcf,
      all_multi_revised_list_1=CleanVcf5_20.all_multi_revised_list_1,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_22 {
    input:
      multitagged_vcf=CleanVcf5_19.multitagged_vcf,
      multitagged_vcf_tbi=CleanVcf5_19.multitagged_vcf_tbi,
      dup_multi_revise_vcf=CleanVcf5_20.dup_multi_revise_vcf,
      del_multi_revise_vcf=CleanVcf5_21.del_multi_revise_vcf,
      all_multi_revised_list_2=CleanVcf5_21.all_multi_revised_list_2,
      new_header=CleanVcf5_21.new_header,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_23 {
    input:
      multitagged_vcf=CleanVcf5_19.multitagged_vcf,
      multitagged_vcf_tbi=CleanVcf5_19.multitagged_vcf_tbi,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_24 {
    input:
      multi_bed=CleanVcf5_23.multi_bed,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_25 {
    input:
      multi_bed_overlap=CleanVcf5_24.multi_bed_overlap,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_26 {
    input:
      multitagged_geno_vcf=CleanVcf5_22.multitagged_geno_vcf,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_27 {
    input:
      multitagged_geno_vcf=CleanVcf5_22.multitagged_geno_vcf,
      multi_remove=CleanVcf5_25.multi_remove,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_28 {
    input:
      cleantagandmulti_vcf=CleanVcf5_27.cleantagandmulti_vcf,
      sv_pipeline_docker=sv_pipeline_docker
  }

  if (CleanVcf5_28.out > 0) {
    call CleanVcf5_29_TRUE_1 {
      input:
        famfile=ped_file,
        cleantagandmulti_vcf=CleanVcf5_27.cleantagandmulti_vcf,
        sv_pipeline_docker=sv_pipeline_docker
    }

    call CleanVcf5_29_TRUE_2 {
      input:
        malecols=CleanVcf5_29_TRUE_1.malecols,
        sex_chr_revise=sex_chr_revise,
        cleantagandmulti_vcf=CleanVcf5_27.cleantagandmulti_vcf,
        sv_pipeline_docker=sv_pipeline_docker
    }

    call CleanVcf5_29_TRUE_3 {
      input:
        cleantagandmulti_vcf=CleanVcf5_27.cleantagandmulti_vcf,
        sexchr_backtoorig=CleanVcf5_29_TRUE_2.sexchr_backtoorig,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  File cleansexcn_vcf_ = select_first([CleanVcf5_29_TRUE_3.cleansexcn_vcf, CleanVcf5_27.cleantagandmulti_vcf])
  call CleanVcf5_30 {
    input:
      cleansexcn_vcf=cleansexcn_vcf_,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf5_31 {
    input:
      cleansexcn_vcf=cleansexcn_vcf_,
      blankcheck_ids=CleanVcf5_30.blankcheck_ids,
      sv_pipeline_docker=sv_pipeline_docker
  }

  output {
    File polished = CleanVcf5_31.polished
  }
}

task CleanVcf5_1 {
  input {
    File normal_revise_vcf
    File revise_vcf_lines
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([normal_revise_vcf, revise_vcf_lines], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 20.0),
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
    cat <(zcat ~{normal_revise_vcf}|fgrep -wvf <(zcat ~{revise_vcf_lines}|awk '{if ($1!="") print $3}'|sort -u)) \
    <(zcat ~{revise_vcf_lines}|awk '{if ($1!="") print}' |tr ' ' '\t') \
    |vcf-sort \
    |bgzip \
    >overlap.revise.vcf.gz
  >>>

  output {
    File overlap_revise_vcf="overlap.revise.vcf.gz"
  }
}

task CleanVcf5_2 {
  input {
    File overlap_revise_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(overlap_revise_vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    ##create bed of VCF##
    svtk vcf2bed ~{overlap_revise_vcf} overlap.revise.bed
    gzip overlap.revise.bed
  >>>

  output {
    File overlap_revise_bed="overlap.revise.bed.gz"
  }
}

task CleanVcf5_3 {
  input {
    File overlap_revise_vcf
    File? outlier_samples_list
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([overlap_revise_vcf, overlap_revise_vcf], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    ~{if defined(outlier_samples_list) then "ln ~{outlier_samples_list} outliers.txt" else "touch outliers.txt"}
    ##multi check##
    zcat ~{overlap_revise_vcf} \
    |awk '{if ($1!~"#") $1=$3;print}' OFS="\t" \
    |vcftools --vcf - --remove outliers.txt --stdout --extract-FORMAT-info RD_CN \
    |gzip \
    >copystate.RD_CN.FORMAT.gz
  >>>

  output {
    File copystate_rd_cn_format="copystate.RD_CN.FORMAT.gz"
  }
}

task CleanVcf5_4 {
  input {
    File overlap_revise_vcf
    File? outlier_samples_list
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([overlap_revise_vcf, outlier_samples_list], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    ~{if defined(outlier_samples_list) then "ln ~{outlier_samples_list} outliers.txt" else "touch outliers.txt"}
    zcat ~{overlap_revise_vcf} \
    |awk '{if ($1!~"#") $1=$3;print}' OFS="\t" \
    |vcftools --vcf - --remove outliers.txt --stdout --extract-FORMAT-info GT \
    |gzip \
    >genotype.gt.FORMAT.gz
  >>>

  output {
    File genotype_gt_format="genotype.gt.FORMAT.gz"
  }
}

task CleanVcf5_5 {
  input {
    File copystate_rd_cn_format
    File overlap_revise_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([copystate_rd_cn_format, overlap_revise_bed], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    ##New method for determining copy state based on >1% of people having an multi-allelic copy state as define above##
    vf_1=$(zcat ~{copystate_rd_cn_format}|awk 'NR==1{print (NF-2) * 0.01}'|awk '{if ($1<=1) print 2; else print }' )

    zcat ~{copystate_rd_cn_format} \
    |{ fgrep -wf <(zcat ~{overlap_revise_bed} |awk -F"\t" '{if ($5=="DEL" && $3-$2>=1000) print $4}' ) || [[ $? == 1 ]]; } \
    |awk 'NR>1{for(i=3;i<=NF;i++) if ($i!="." && $i>3) print  $1 }' \
    |sort \
    |uniq -c \
    |awk -v vf_1=$vf_1 '{if ($1>vf_1)print $2}' \
    |gzip \
    >multi.del.ids.txt.gz
  >>>

  output {
    File multi_del_ids="multi.del.ids.txt.gz"
  }
}

task CleanVcf5_6 {
  input {
    File copystate_rd_cn_format
    File overlap_revise_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([copystate_rd_cn_format, overlap_revise_bed], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    vf_1=$(zcat ~{copystate_rd_cn_format}|awk 'NR==1{print (NF-2) * 0.01}'|awk '{if ($1<=1) print 2; else print }' )

    zcat ~{copystate_rd_cn_format} \
    |{ fgrep -wf <(zcat ~{overlap_revise_bed}|awk -F"\t" '{if ($5=="DUP" && $3-$2>=1000) print $4}' ) || [[ $? == 1 ]]; } \
    |awk 'NR>1{for(i=3;i<=NF;i++) if ($i!="." && $i>4) print  $1 }' \
    |sort \
    |uniq -c \
    |awk -v vf_1=$vf_1 '{if ($1>vf_1)print $2}' \
    >multi.dup.ids.txt
  >>>

  output {
    File multi_dup_ids_1="multi.dup.ids.txt"
  }
}

task CleanVcf5_7 {
  input {
    File copystate_rd_cn_format
    File overlap_revise_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([copystate_rd_cn_format, overlap_revise_bed], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    ##Case with CN 0,1,2,3,4##
    zcat ~{copystate_rd_cn_format} \
    |{ fgrep -wf <(zcat ~{overlap_revise_bed} | awk -F"\t" '{if ($5=="DUP" && $3-$2>=1000) print $4}') || [[ $? == 1 ]]; } \
    |awk 'NR>1{for(i=3;i<=NF;i++) if ($i!="." && ($i<1 || $i>4)) print  $1 "\t" $i }'\
    |sort -u \
    |awk  '{print $1}' \
    |sort \
    |uniq -c \
    |awk '{if ($1>4) print $2}'>gt4copystate.txt
  >>>

  output {
    File gt4copystate="gt4copystate.txt"
  }
}

task CleanVcf5_8 {
  input {
    File copystate_rd_cn_format
    File overlap_revise_bed
    File gt4copystate
    File multi_dup_ids_1
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([copystate_rd_cn_format, overlap_revise_bed, gt4copystate, multi_dup_ids_1], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    vf_1=$(zcat ~{copystate_rd_cn_format}|awk 'NR==1{print (NF-2) * 0.01}'|awk '{if ($1<=1) print 2; else print }' )

    mv ~{multi_dup_ids_1} multi.dup.ids.txt
    zcat ~{copystate_rd_cn_format} \
      | awk -F'\t' -v OFS='\t' 'ARGIND==1{inFileA[$1]; next} {if ($1 in inFileA) print }' <(zcat ~{overlap_revise_bed}|awk -F"\t" '{if ($5=="DUP" && $3-$2>=1000) print $4}') - \
      | awk 'NR>1{for(i=3;i<=NF;i++) if ($i!="." && ($i<1 || $i>4)) print  $1 }' \
      | sort \
      | uniq -c \
      | awk -F'\t' -v OFS='\t' 'ARGIND==1{inFileA[$1]; next} {if ($1 in inFileA) print }' ~{gt4copystate} - \
      | awk -v vf_1=$vf_1 '{if ($1>vf_1)print $2}' \
      >> multi.dup.ids.txt

    sort -u multi.dup.ids.txt |gzip >multi.dup.ids.txt.gz
  >>>

  output {
    File multi_dup_ids="multi.dup.ids.txt.gz"
  }
}

task CleanVcf5_9 {
  input {
    File overlap_revise_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(overlap_revise_vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    ##Regenotype to determine multiallelic; we just change copy state for some nested variants and we need to make sure we get proper genotype for these; also previous stages have different notaion for multiallelic and we need to make this uniform; this is a CN based regenotyping so restricted to >5kb ##
    ##Genotype big dup##
    svtk vcf2bed ~{overlap_revise_vcf} regeno.bed
    gzip regeno.bed
  >>>

  output {
    File regeno_bed="regeno.bed.gz"
  }
}

task CleanVcf5_10 {
  input {
    File genotype_gt_format
    File multi_dup_ids
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([genotype_gt_format, multi_dup_ids], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    ##add variants that are <5kb because clustering but have a mutliallelic genotype from before##
    zcat ~{genotype_gt_format} \
      |awk '{if ($1~"DUP") print}' \
      |awk '{for (i = 3; i <= NF; ++i) print $1 "\t" $i}' \
      |awk '{if ($2!="1/1" && $2!="0/0" && $2!="0/1" && $2!="./.") print $1}' \
      |{ fgrep -wvf <(zcat ~{multi_dup_ids}) || [[ $? == 1 ]]; } \
      |sort -u \
      >gt5kb.dup.ids.txt
  >>>

  output {
    File gt5kb_dup_ids_1="gt5kb.dup.ids.txt"
  }
}

task CleanVcf5_11 {
  input {
    File genotype_gt_format
    File multi_del_ids
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([genotype_gt_format, multi_del_ids], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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

    zcat ~{genotype_gt_format} \
      |awk '{if ($1~"DEL") print}' \
      |awk '{for (i = 3; i <= NF; ++i) print $1 "\t" $i}' \
      |awk '{if ($2!="1/1" && $2!="0/0" && $2!="0/1" && $2!="./.") print $1}' \
      |{ fgrep -wvf <(zcat ~{multi_del_ids}) || [[ $? == 1 ]]; } \
      |sort -u \
      >gt5kb.del.ids.txt
  >>>

  output {
    File gt5kb_del_ids_1="gt5kb.del.ids.txt"
  }
}

task CleanVcf5_12 {
  input {
    File multi_dup_ids
    File regeno_bed
    File gt5kb_dup_ids_1
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([multi_dup_ids, regeno_bed, gt5kb_dup_ids_1], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    mv ~{gt5kb_dup_ids_1} gt5kb.dup.ids.txt
    ##generate list##
    ##CNV >5kb, split del and dup ##
    if [ -f ~{multi_dup_ids} ]
    then
      zcat ~{regeno_bed}  \
      |awk '{if ($3-$2>=5000 && $5=="DUP")print $4}' \
      |{ fgrep -wvf <(zcat ~{multi_dup_ids}) || [[ $? == 1 ]]; } \
      >>gt5kb.dup.ids.txt
    else
      zcat ~{regeno_bed} \
      |awk '{if ($3-$2>=5000 && $5=="DUP")print $4}' \
      >>gt5kb.dup.ids.txt
    fi
  >>>

  output {
    File gt5kb_dup_ids="gt5kb.dup.ids.txt"
  }
}

task CleanVcf5_13 {
  input {
    File multi_del_ids
    File gt5kb_del_ids_1
    File regeno_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([multi_del_ids, gt5kb_del_ids_1, regeno_bed], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    mv ~{gt5kb_del_ids_1} gt5kb.del.ids.txt
    if [ -f ~{multi_del_ids} ]
    then
      zcat ~{regeno_bed} \
        |awk '{if ($3-$2>=5000 && $5=="DEL")print $4}' \
        |{ fgrep -wvf <(zcat ~{multi_del_ids}) || [[ $? == 1 ]]; } \
        >>gt5kb.del.ids.txt
    else
      zcat ~{regeno_bed} \
        |awk '{if ($3-$2>=5000 && $5=="DEL")print $4}' \
        >>gt5kb.del.ids.txt
    fi

  >>>

  output {
    File gt5kb_del_ids="gt5kb.del.ids.txt"
  }
}

task CleanVcf5_14 {
  input {
    File overlap_revise_vcf
    File gt5kb_dup_ids
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([overlap_revise_vcf, gt5kb_dup_ids], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    zcat ~{overlap_revise_vcf} \
      |{ fgrep -wf ~{gt5kb_dup_ids} || [[ $? == 1 ]]; } \
      >dup.int.txt
  >>>

  output {
    File dup_int="dup.int.txt"
  }
}

task CleanVcf5_15 {
  input {
    File overlap_revise_vcf
    File gt5kb_del_ids
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([overlap_revise_vcf, gt5kb_del_ids], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    zcat ~{overlap_revise_vcf} \
      |{ fgrep -wf ~{gt5kb_del_ids} || [[ $? == 1 ]]; } \
      >>del.int.txt
  >>>

  output {
    File del_int="del.int.txt"
  }
}

task CleanVcf5_16 {
  input {
    File del_int
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(del_int, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    ##regenotype VCF##
    dellen=$(cat ~{del_int}|wc -l)
    columnlen=$(less ~{del_int}|cut -f10-|tr '\t' '\n' |wc -l)
    dellenchange=$(echo $dellen $columnlen|awk '{if ($1 == 0) { print "0" } else { print $2/$1}}')

    paste <(less ~{del_int}|cut -f1-9) <(less ~{del_int}|cut -f10-|tr '\t' '\n' \
    |awk -F':' '{if ($3>=2 && $1!="./.") $1="0/0"; \
    else if ($3==1 && $1!="./.") $1="0/1"; \
    else if ($1!="./.")$1="1/1";print}' OFS=":" \
    |awk -v lenchange=$dellenchange 'NR%lenchange {printf("%s\t", $0); next} \
    {print $0}')>del.revise.txt
  >>>

  output {
    File del_revise="del.revise.txt"
  }
}

task CleanVcf5_17 {
  input {
    File dup_int
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(dup_int, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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

    duplen=$(cat ~{dup_int}|wc -l)
    columnlen=$(less ~{dup_int}|cut -f10-|tr '\t' '\n' |wc -l)
    duplenchange=$(echo $duplen $columnlen|awk '{if ($1 == 0) { print "0" } else { print $2/$1}}')

    paste <(less ~{dup_int}|cut -f1-9) <(less ~{dup_int}|cut -f10-|tr '\t' '\n' \
    |awk -F':' '{if ($3<=2 && $1!="./.") $1="0/0"; \
    else if ($3==3 && $1!="./.") $1="0/1"; \
    else if ($1!="./.") $1="1/1";print}' OFS=":" \
    |awk -v lenchange=$duplenchange 'NR%lenchange {printf("%s\t", $0); next} \
    {print $0}') >dup.revise.txt
  >>>

  output {
    File dup_revise="dup.revise.txt"
  }
}

task CleanVcf5_18 {
  input {
    File overlap_revise_vcf
    File gt5kb_dup_ids
    File gt5kb_del_ids
    File dup_revise
    File del_revise
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([overlap_revise_vcf, gt5kb_dup_ids, gt5kb_del_ids, dup_revise, del_revise], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
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
    cat <(zcat ~{overlap_revise_vcf}|fgrep -wvf <(cat ~{gt5kb_dup_ids} ~{gt5kb_del_ids})) \
    <(cat ~{dup_revise} ~{del_revise}) \
    |vcf-sort \
    |bgzip \
    >newdepth.geno.vcf.gz
  >>>

  output {
    File newdepth_geno_vcf="newdepth.geno.vcf.gz"
  }
}

task CleanVcf5_19 {
  input {
    File multi_geno_ids_txt
    File multi_del_ids
    File multi_dup_ids
    File newdepth_geno_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([multi_geno_ids_txt, multi_del_ids, multi_dup_ids, newdepth_geno_vcf, newdepth_geno_vcf], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    ##Tag multi##
    ##Add filters to header##
    zcat ~{newdepth_geno_vcf} \
      |awk -F'\t' -v OFS='\t' 'ARGIND==1{inFileA[$1]; next} {if ($3 in inFileA && $1!~"#" && $7!~"PESR_GT_OVERDISPERSION") $7=$7";PESR_GT_OVERDISPERSION"; print }' \
        <(cat <(zcat ~{multi_geno_ids_txt}) <(printf "\n")) - \
      |awk -F'\t' -v OFS='\t' 'ARGIND==1{inFileA[$1]; next} {if ($3 in inFileA && $1!~"#") $7=$7";MULTIALLELIC"; print }' \
        <(cat <(zcat ~{multi_del_ids} ~{multi_dup_ids} |sort -u) <(printf "\n")) - \
      |sed 's\PASS;\\g' \
      |awk '{if (NR==2) print $0 "\n" "##FILTER=<ID=PESR_GT_OVERDISPERSION,Description=\"High PESR dispersion count\">" ;else print}' \
      |awk '{if (NR==2) print $0 "\n" "##FILTER=<ID=MULTIALLELIC,Description=\"Multiallelic copy number variant>" ;else print}' \
      |bgzip \
      >multitagged.vcf.gz
    tabix multitagged.vcf.gz
  >>>

  output {
    File multitagged_vcf="multitagged.vcf.gz"
    File multitagged_vcf_tbi="multitagged.vcf.gz.tbi"
  }
}

task CleanVcf5_20 {
  input {
    File multi_dup_ids
    File multitagged_vcf
    File multitagged_vcf_tbi
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([multi_dup_ids, multitagged_vcf], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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

    # use BCFTOOLS 1.9
    BCFTOOLS=/usr/local/bin/bcftools

    touch all.multi.revised.list
    touch dup.multi.revise.vcf
    if [ $(zcat ~{multi_dup_ids}|wc -l) -ge 1  ]
    then
      /opt/sv-pipeline/04_variant_resolution/scripts/reset_multiallelic_format_fields.py ~{multitagged_vcf} <(zcat ~{multi_dup_ids}) > dup.multi.revise.vcf
      ${BCFTOOLS} query -f '%ID\n' dup.multi.revise.vcf >> all.multi.revised.list
    fi
  >>>

  output {
    File dup_multi_revise_vcf="dup.multi.revise.vcf"
    File all_multi_revised_list_1="all.multi.revised.list"
  }
}

task CleanVcf5_21 {
  input {
    File multi_del_ids
    File multitagged_vcf
    File multitagged_vcf_tbi
    File dup_multi_revise_vcf
    File all_multi_revised_list_1
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([multi_del_ids, multitagged_vcf, dup_multi_revise_vcf, all_multi_revised_list_1], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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

    # use BCFTOOLS 1.9
    BCFTOOLS=/usr/local/bin/bcftools

    mv ~{all_multi_revised_list_1} all.multi.revised.list
    touch del.multi.revise.vcf
    if [ $(zcat ~{multi_del_ids}|wc -l) -ge 1 ]
    then
      /opt/sv-pipeline/04_variant_resolution/scripts/reset_multiallelic_format_fields.py ~{multitagged_vcf} <(zcat ~{multi_del_ids}) > del.multi.revise.vcf
      ${BCFTOOLS} query -f '%ID\n' del.multi.revise.vcf >> all.multi.revised.list
    fi

    # make sure that the new header includes CN and CNQ format fields if we set any
    if [ -s ~{dup_multi_revise_vcf} ]
    then
      grep '^#' ~{dup_multi_revise_vcf} > new_header.vcf
    elif [ -s  del.multi.revise.vcf ]
    then
      grep '^#' del.multi.revise.vcf > new_header.vcf
    else
      zcat ~{multitagged_vcf} | grep '^#' > new_header.vcf
    fi
  >>>

  output {
    File del_multi_revise_vcf="del.multi.revise.vcf"
    File all_multi_revised_list_2="all.multi.revised.list"
    File new_header="new_header.vcf"
  }
}


task CleanVcf5_22 {
  input {
    File multitagged_vcf
    File multitagged_vcf_tbi
    File dup_multi_revise_vcf
    File del_multi_revise_vcf
    File all_multi_revised_list_2
    File new_header
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([multitagged_vcf, dup_multi_revise_vcf, del_multi_revise_vcf, all_multi_revised_list_2, new_header], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 20.0),
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

    # use BCFTOOLS 1.9
    BCFTOOLS=/usr/local/bin/bcftools

    # combine the revised variants with the unrevised variants, reheader, resort, and compress
    cat <(zcat ~{multitagged_vcf} | fgrep -wvf ~{all_multi_revised_list_2}) \
      <(cat ~{del_multi_revise_vcf} ~{dup_multi_revise_vcf} | grep -v '^#' | awk '!seen[$3]++') \
      |${BCFTOOLS} reheader -h ~{new_header} \
      |vcf-sort \
      |bgzip \
      >multitagged.geno.vcf.gz
  >>>

  output {
    File multitagged_geno_vcf="multitagged.geno.vcf.gz"
  }
}

task CleanVcf5_23 {
  input {
    File multitagged_vcf
    File multitagged_vcf_tbi
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(multitagged_vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    ##remove overlapping multi###
    zcat ~{multitagged_vcf} \
    |awk -F'\t' '{if ($1~"#" || ($7~"MULTIALLELIC" &&  ($5=="<DEL>" || $5=="<DUP>"))) print}' \
    |svtk vcf2bed stdin tmp.bed
    cut -f1-5 tmp.bed \
    |gzip \
    >multi.bed.gz
  >>>

  output {
    File multi_bed="multi.bed.gz"
  }
}

task CleanVcf5_24 {
  input {
    File multi_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(multi_bed, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    ##strip out overlapping multiallelics##
    bedtools intersect -wa -wb -a  ~{multi_bed} -b ~{multi_bed} \
    |awk -F'\t' '{if ($4!=$9 && $3-$2>=$8-$7) print $0; \
    else if ($4!=$9) print $6,$7,$8,$9,$10,$1,$2,$3,$4,$5}' OFS="\t" \
    |sort -u \
    |awk '{print $3-$2,$8-$7,$0}' OFS="\t"  \
    |sort -nrk1,1 -k2,2nr \
    |cut -f3- \
    >multi.bed.overlap.txt
  >>>

  output {
    File multi_bed_overlap="multi.bed.overlap.txt"
  }
}

task CleanVcf5_25 {
  input {
    File multi_bed_overlap
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(multi_bed_overlap, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    set -euo pipefail
    echo "">multi.remove.txt

    while read bed
    do
      echo "$bed"|cut -d$'\t' -f1-5 >large.bed
      echo "$bed"|cut -d$'\t' -f6-10>small.bed
      overlap=$(bedtools coverage -a small.bed -b large.bed|awk '{if ($NF>0.50) print "YES";else print "NO"}')
      echo $bed|awk '{print $4}'
      if [ "$overlap" == "YES" ] && [ $(awk '{print $4}' large.bed|fgrep -wf - multi.remove.txt|wc -l) -eq 0 ]
      then
        awk '{print $4}' small.bed >>multi.remove.txt
      fi
    done< ~{multi_bed_overlap}
  >>>

  output {
    File multi_remove="multi.remove.txt"
  }
}

task CleanVcf5_26 {
  input {
    File multitagged_geno_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(multitagged_geno_vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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

    # use BCFTOOLS 1.9
    BCFTOOLS=/usr/local/bin/bcftools

    ##get alt tag for multiallelics##
    ## produces a file with a row for each distinct multialllic variant ID and copy number combination
    ${BCFTOOLS} query -i 'FILTER = "MULTIALLELIC"' -f '[%ID\t%CN\n]' ~{multitagged_geno_vcf} \
    |sort -u >multi.cn.txt
  >>>

  output {
    File multi_cn="multi.cn.txt"
  }
}

task CleanVcf5_27 {
  input {
    File multitagged_geno_vcf
    File multi_remove
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([multitagged_geno_vcf, multi_remove], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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

    # use BCFTOOLS 1.9
    BCFTOOLS=/usr/local/bin/bcftools

    ##strip out variants with no genotypes and overlapping multiallelics##
    ### Find missing genotype and then add multiallelics that need to be removed###
    ##change multiallelics svtype into mCNV##
    ##add CN information to ALT column##
    zcat ~{multitagged_geno_vcf} \
    |${BCFTOOLS} view -e 'FILTER == "MULTIALLELIC"'  \
    |svtk vcf2bed stdin tmp.bed

    awk -F'\t' '{if ($6=="") print $4}' tmp.bed \
    |cat - ~{multi_remove} \
    |sed '/^$/d' \
    |{ fgrep -wvf - <(zcat ~{multitagged_geno_vcf} ) || [[ $? == 1 ]]; } \
    |awk -F';' '{if ($1~"MULTIALLELIC" && ( $2~"DEL" || $2~"DUP")) $2="SVTYPE=CNV"; print}' OFS=';' \
    |awk '{OFS="\t"; if ($8~"SVTYPE=CNV;") $5="<CNV>"; print}' \
    |bgzip \
    >cleantagandmulti.vcf.gz
  >>>

  output {
    File cleantagandmulti_vcf="cleantagandmulti.vcf.gz"
  }
}

task CleanVcf5_28 {
  input {
    File cleantagandmulti_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(cleantagandmulti_vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    zcat ~{cleantagandmulti_vcf}|awk '{if (($1~"X" || $1~"Y") && $1!~"#") print}'|wc -l > out
  >>>

  output {
    Int out=read_int("out")
  }
}

task CleanVcf5_29_TRUE_1 {
  input {
    File famfile
    File cleantagandmulti_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([famfile, cleantagandmulti_vcf], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    ##Determine columns male columns##
    zcat ~{cleantagandmulti_vcf}\
    |egrep ^# \
    |tail -n 1 \
    |tr '\t' '\n' \
    |cat -n - \
    >col.txt

    awk '{if ($5==1) print $2}' ~{famfile} \
    |{ fgrep -wf - col.txt || [[ $? == 1 ]]; } \
    >malecols.txt
  >>>

  output {
    File malecols="malecols.txt"
  }
}

task CleanVcf5_29_TRUE_2 {
  input {
    File sex_chr_revise
    File cleantagandmulti_vcf
    File malecols
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([sex_chr_revise, cleantagandmulti_vcf, malecols], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    ##regenotype male calls on sex chr and add 1 to copy state for multialleic check##

    python3 <<CODE | bgzip > sexchr.backtoorig.txt.gz
    import pysam
    import sys

    with open("~{malecols}") as f:
      samples = [x.strip().split('\t')[1] for x in f.readlines() if x]

    with open("~{sex_chr_revise}") as f:
      vids = set([x.strip() for x in f.readlines() if x])

    vcf = pysam.VariantFile("~{cleantagandmulti_vcf}")

    for record in vcf:
      if record.id not in vids:
        continue
      for i in samples:
        g = record.samples[i]
        if g['RD_CN'] is not None and g['RD_CN'] >= 1:
          g['RD_CN'] = g['RD_CN'] - 1
      sys.stdout.write(str(record))

    vcf.close()
    CODE

  >>>

  output {
    File sexchr_backtoorig="sexchr.backtoorig.txt.gz"
  }
}

task CleanVcf5_29_TRUE_3 {
  input {
    File cleantagandmulti_vcf
    File sexchr_backtoorig
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([cleantagandmulti_vcf, sexchr_backtoorig], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 50.0),
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
    cat <(zcat ~{cleantagandmulti_vcf}|fgrep -wvf <(zcat ~{sexchr_backtoorig}|awk '{print $3}'  )) \
    <(zcat ~{sexchr_backtoorig} |awk '{if ($1!="") print}' |tr ' ' '\t') \
    |vcf-sort \
    |bgzip \
    >cleansexCN.vcf.gz
  >>>

  output {
    File cleansexcn_vcf="cleansexCN.vcf.gz"
  }
}

task CleanVcf5_30 {
  input {
    File cleansexcn_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(cleansexcn_vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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

    mv ~{cleansexcn_vcf} cleanGQ.vcf.gz
    ##find blank variants with no samples##
    svtk vcf2bed cleanGQ.vcf.gz tmp.bed

    awk -F'\t' '{if ($5!~"CN" && $6=="") print $4}' tmp.bed \
    >blankcheck.ids.txt
  >>>

  output {
    File blankcheck_ids="blankcheck.ids.txt"
  }
}

task CleanVcf5_31 {
  input {
    File cleansexcn_vcf
    File blankcheck_ids
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([cleansexcn_vcf, blankcheck_ids], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
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
    mv ~{cleansexcn_vcf} cleanGQ.vcf.gz
    ##Fix header##
    ##get header to clean##
    ##add new filters##
    zcat cleanGQ.vcf.gz \
    |awk '{if ($1~"##" && NR>1)  print}' \
    |{ fgrep -v "MULTIALLELIC" || [[ $? == 1 ]]; } \
    |awk '{if (NR==2) print $0 "\n" "##FILTER=<ID=MULTIALLELIC,Description=\"Multiallelic site\">" ;else print}' \
    |awk '{if (NR==2) print $0 "\n" "##ALT=<ID=CNV,Description=\"Copy Number Polymorphism\">" ;else print}' \
    |sort -k1,1 \
    |{ egrep -v "CIPOS|CIEND|RMSSTD|EVENT|INFO=<ID=UNRESOLVED,|source|varGQ|bcftools|ALT=<ID=UNR" || [[ $? == 1 ]]; } \
    |cat <(zcat cleanGQ.vcf.gz|head -n 1) - <(zcat cleanGQ.vcf.gz|fgrep -wvf ~{blankcheck_ids} |awk '{if ($1!~"##")  print}') \
    |bgzip >polished.vcf.gz
  >>>

  output {
    File polished="polished.vcf.gz"
  }
}
