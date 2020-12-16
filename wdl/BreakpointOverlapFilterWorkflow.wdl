version 1.0

import "Structs.wdl"

workflow BreakpointOverlapFilterWorkflow {
  input {
    File vcf
    String prefix
    File bothside_pass
    File background_fail
    String sv_pipeline_docker
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override_1a
    RuntimeAttr? runtime_attr_override_1b
    RuntimeAttr? runtime_attr_override_2
    RuntimeAttr? runtime_attr_override_3
    RuntimeAttr? runtime_attr_override_4
    RuntimeAttr? runtime_attr_override_5
    RuntimeAttr? runtime_attr_override_6
    RuntimeAttr? runtime_attr_override_7
  }

  call BreakpointOverlapFilter1a {
    input:
      vcf=vcf,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_override_1a
  }

  call BreakpointOverlapFilter1b {
    input:
      variants_bed=BreakpointOverlapFilter1a.variants_bed,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_override_1b
  }

  call BreakpointOverlapFilter2 {
    input:
      dupside1=BreakpointOverlapFilter1b.dupside1,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_override_2
  }

  call BreakpointOverlapFilter3 {
    input:
      dupside1_freq50=BreakpointOverlapFilter2.dupside1_freq50,
      bothside_pass=bothside_pass,
      background_fail=background_fail,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_override_3
  }

  call BreakpointOverlapFilter4 {
    input:
      dupside1=BreakpointOverlapFilter1b.dupside1,
      dupside1_bothpassfilter=BreakpointOverlapFilter3.dupside1_bothpassfilter,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_override_4
  }

  call BreakpointOverlapFilter5 {
    input:
      dupside1_samplecountfilter=BreakpointOverlapFilter4.dupside1_samplecountfilter,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_override_5
  }

  call BreakpointOverlapFilter6 {
    input:
      dupside1_allfilter=BreakpointOverlapFilter5.dupside1_allfilter,
      dupside1_freq50=BreakpointOverlapFilter2.dupside1_freq50,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_override_6
  }

  call BreakpointOverlapFilter7 {
    input:
      vcf=vcf,
      prefix=prefix,
      remove_side1_var=BreakpointOverlapFilter6.remove_side1_var,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_override_7
  }

  output {
    File out_vcf = BreakpointOverlapFilter7.bp_filtered_vcf
    File out_vcf_index = BreakpointOverlapFilter7.bp_filtered_vcf_idx
  }
}

task BreakpointOverlapFilter1a {
  input {
    File vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GiB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 1.5,
                                  disk_gb: ceil(10 + input_size * 2),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    ##clean out variants that overlap at one site##
    ##pull out variants with duplicate bp that are not driven by depth which will be integrated in the clean vcf##
    ##make sure to flip bed as well so second bp location can be compared with first from other variants##
    svtk vcf2bed ~{vcf} variants.bed.gz -i CHR2 -i STRANDS -i SVLEN -i varGQ -i END -i EVIDENCE -i SVTYPE --split-bnd
  >>>

  output {
    File variants_bed = "variants.bed.gz"
  }
}

task BreakpointOverlapFilter1b {
  input {
    File variants_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(variants_bed, "GiB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 60,
                                  disk_gb: ceil(10 + input_size * 2),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    ##clean out variants that overlap at one site##
    ##pull out variants with duplicate bp that are not driven by depth which will be integrated in the clean vcf##
    ##make sure to flip bed as well so second bp location can be compared with first from other variants##
    zcat ~{variants_bed}  \
      | sed "s/+-/+ "$'\t -/g' \
      | sed "s/-+/- "$'\t +/g' \
      | sed "s/++/+ "$'\t +/g' \
      | sed "s/--/- "$'\t -/g' | \
      ##Convert back to 1-based positions##
      awk -v OFS='\t' '{$2=$2+1; print $0}' \
      | awk -v OFS='\t' \
        '{if (!(($NF=="DEL" || $NF=="DUP") && $10>=5000)) print $0 "\n" $7,$12,$2,$4,$5,$6,$1,$9,$8,$10,$11,$2,$13,$14 }' | \
      ###Find duplicated variants that overlap at same bp one side##
      awk 'cnt[$1"_"$2"_"$8]++{if (cnt[$1"_"$2"_"$8]==2) print prev[$1"_"$2"_"$8] "\t" $1"_"$2"_"$8 \
        ; print $0 "\t" $1"_"$2"_"$8} {prev[$1"_"$2"_"$8]=$0}' \
      | awk '!seen[$4"_"$NF]++' \
      | awk 'cnt[$NF]++{if (cnt[$NF]==2) print prev[$NF] \
        ; print $0 } {prev[$NF]=$0}' \
      > dupside1.bed
  >>>

  output {
    File dupside1 = "dupside1.bed"
  }
}

task BreakpointOverlapFilter2 {
  input {
    File dupside1
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(dupside1, "GiB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 1.5,
                                  disk_gb: ceil(10 + input_size * 15),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    mkdir temp
    awk '{print $NF "\t" $6}' ~{dupside1} \
      | awk -F'[,\t]' '{for (i=2;i<=NF;i++) print $1 "\t" $i}' \
      | sort -T temp \
      | uniq -D \
      | awk '{print $1}'|sort|uniq -c \
      > dupside1.A.txt

    awk '{print $NF "\t" $6}' ~{dupside1} \
      | awk -F'[,\t]' '{for (i=2;i<=NF;i++) print $1 "\t" $i}' \
      | awk '{print $1}' \
      | sort -T temp \
      | uniq -c \
      > dupside1.B.txt

    ##Find 50% overlap between samples for overlaps##
    join -j 2 dupside1.A.txt dupside1.B.txt \
      | awk '{if ($2 >= 0.5 * $3) print $1}' \
      | (fgrep -wf - ~{dupside1} || printf "") \
      > dupside1.freq50.txt

  >>>

  output {
    File dupside1_freq50 = "dupside1.freq50.txt"
  }
}

task BreakpointOverlapFilter3 {
  input {
    File dupside1_freq50
    File background_fail
    File bothside_pass
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([dupside1_freq50, background_fail, bothside_pass], "GiB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 1.0 + input_size * 2.5,
                                  disk_gb: ceil(10 + input_size * 30),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    ##Add SRfail###
    { fgrep -wf <(awk '{print $NF}' ~{background_fail}) ~{dupside1_freq50} || true; } \
      | awk '{print $0 "\t" 0}' \
      > dupside1.passSR.txt

    { fgrep -wvf <(awk '{print $NF}' ~{background_fail}) ~{dupside1_freq50} || true; } \
      | awk '{print $0 "\t" 1}' \
      >> dupside1.passSR.txt
    rm ~{background_fail}

    mkdir temp
    sort -k4,4 -T temp dupside1.passSR.txt > dupside1.passSR.sorted.txt
    rm dupside1.passSR.txt

    awk '{print $NF "\t" $1}' ~{bothside_pass} | sort -T temp -k1,1 > bothside_pass.sorted.txt
    rm ~{bothside_pass}

    ##Attach the % of variants that show SR support at bothends##
    join -1 4 -2 1 -e "0" -a 1 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 2.2 \
      dupside1.passSR.sorted.txt bothside_pass.sorted.txt \
      | tr ' ' '\t' \
      | bgzip \
      > dupside1.bothpassfilter.txt.gz
  >>>

  output {
    File dupside1_bothpassfilter = "dupside1.bothpassfilter.txt.gz"
  }
}


task BreakpointOverlapFilter4 {
  input {
    File dupside1
    File dupside1_bothpassfilter
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([dupside1, dupside1_bothpassfilter], "GiB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 60,
                                  disk_gb: ceil(10 + input_size * 5),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    mkdir temp
    awk '{print $4 "\t" $6}' ~{dupside1} \
      | awk -F'[,\t]' '{print $1 "\t" NF-1}' \
      | sort -T temp -k1,1 \
      > dupside1.A.txt
    rm -r temp ~{dupside1}

    ##count number of samples and indiciate if size gt 50bp##
    zcat ~{dupside1_bothpassfilter} > dupside1_bothpassfilter.txt
    rm ~{dupside1_bothpassfilter}
    join -1 4 -2 1 dupside1_bothpassfilter.txt dupside1.A.txt \
      | tr ' ' '\t' \
      | awk '{if ($10>=50) print $0 "\t" 1;else print $0 "\t" 0}' \
      | bgzip \
      > dupside1.samplecountfilter.txt.gz
  >>>

  output {
    File dupside1_samplecountfilter = "dupside1.samplecountfilter.txt.gz"
  }
}

task BreakpointOverlapFilter5 {
  input {
    File dupside1_samplecountfilter
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(dupside1_samplecountfilter, "GiB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10 + input_size * 5),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    ##Convert Evidence column into Integers for scoring and ##
    ##RD,PE,SR-1,RD,PE-2,PE,SR-3,RD,SR-4,PE-5,RD-6,SR-7##
    zcat ~{dupside1_samplecountfilter} \
      | sed 's/BAF,//g' \
      | awk -v OFS='\t' '
      {
        if ($13=="PE,RD,SR") print $0 "\t" 1
        else if ($13=="PE,RD") print $0 "\t" 2
        else if ($13=="PE,SR") print $0 "\t" 3
        else if ($13=="RD,SR") print $0 "\t" 4
        else if ($13=="PE") print $0 "\t" 5
        else if ($13=="RD") print $0 "\t" 6
        else if ($13=="SR") print $0 "\t" 7
      }' | \
      ##assign BND to bottom
      awk '{if ($14=="BND") print $0 "\t" 0;else print $0 "\t" 1}' \
      | bgzip \
      > dupside1.allfilter.txt.gz
  >>>

  output {
    File dupside1_allfilter = "dupside1.allfilter.txt.gz"
  }
}

task BreakpointOverlapFilter6 {
  input {
    File dupside1_allfilter
    File dupside1_freq50
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([dupside1_allfilter, dupside1_freq50], "GiB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(50 + input_size * 5),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    ##sort file with overlapping samples LevelofSupport->BothEndsupport->SRfail-> Not BND->Higher varq-> Higher Freq -> Smallest size if gt 5kb##
    mkdir temp
    zcat ~{dupside1_allfilter} \
      | sort -T temp -k20,20n -k17,17nr -nrk16,16 -k21,21nr -k11,11nr -k18,18nr  -k19,19nr -k10,10n \
      | awk '!seen[$15]++' \
      | awk '{print $1}' \
      | (fgrep -wvf - ~{dupside1_freq50} || printf "") \
      | awk '{print $4}' \
      | bgzip \
      > remove.side1.var.txt.gz
  >>>

  output {
    File remove_side1_var = "remove.side1.var.txt.gz"
  }
}

task BreakpointOverlapFilter7 {
  input {
    File vcf
    File remove_side1_var
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_file = prefix + ".non_redundant.vcf.gz"

  Float input_size = size([vcf, remove_side1_var], "GiB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10 + input_size * 3),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    ##remove variants with samebp##
    zcat ~{remove_side1_var} > remove.side1.var.txt
    rm ~{remove_side1_var}
    (zgrep -wvf remove.side1.var.txt ~{vcf} || printf "") \
      | bgzip \
      > ~{output_file}

    tabix ~{output_file}
  >>>

  output {
    File bp_filtered_vcf = output_file
    File bp_filtered_vcf_idx = output_file + ".tbi"
  }
}
