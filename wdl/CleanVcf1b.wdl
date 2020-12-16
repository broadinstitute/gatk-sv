version 1.0

import "Structs.wdl"

workflow CleanVcf1b {
  input {
    File intermediate_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override # TODO
  }

  call CleanVcf1b_1 {
    input:
      intermediate_vcf=intermediate_vcf,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf1b_2 {
    input:
      int_bed=CleanVcf1b_1.int_bed,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf1b_3 {
    input:
      int_vcf=intermediate_vcf,
      normoverlap=CleanVcf1b_2.normoverlap,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf1b_4 {
    input:
      int_vcf=intermediate_vcf,
      normoverlap=CleanVcf1b_2.normoverlap,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf1b_5 {
    input:
      normoverlap=CleanVcf1b_2.normoverlap,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf1b_6 {
    input:
      overlap_test=CleanVcf1b_5.overlap_test,
      rd_cn_normcheck=CleanVcf1b_3.rd_cn_normcheck,
      ev_normcheck=CleanVcf1b_4.ev_normcheck,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf1b_7 {
    input:
      geno_normal_revise=CleanVcf1b_6.geno_normal_revise,
      int_vcf=intermediate_vcf,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf1b_8 {
    input:
      subset_vcf=CleanVcf1b_7.subset_vcf,
      geno_normal_revise=CleanVcf1b_6.geno_normal_revise,
      col=CleanVcf1b_1.col,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf1b_9 {
    input:
      normal_revise_vcf_lines=CleanVcf1b_8.normal_revise_vcf_lines,
      int_vcf=intermediate_vcf,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf1b_10 {
    input:
      normal_revise_vcf=CleanVcf1b_9.normal_revise_vcf,
      normal_revise_vcf_csi=CleanVcf1b_9.normal_revise_vcf_csi,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf1b_11 {
    input:
      copystate_rd_cn=CleanVcf1b_10.copystate_rd_cn,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf1b_12 {
    input:
      int_bed=CleanVcf1b_1.int_bed,
      copystate_per_variant=CleanVcf1b_11.copystate_per_variant,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call CleanVcf1b_13 {
    input:
      int_bed=CleanVcf1b_1.int_bed,
      multi_del=CleanVcf1b_12.multi_del,
      copystate_per_variant=CleanVcf1b_11.copystate_per_variant,
      sv_pipeline_docker=sv_pipeline_docker
  }

  output {
    File multi = CleanVcf1b_13.multi
    File normal = CleanVcf1b_9.normal_revise_vcf
    File vcftools_idx = CleanVcf1b_9.normal_revise_vcf_csi
  }
}


task CleanVcf1b_1 {
  input {
    File intermediate_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(intermediate_vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 2.0 + input_size * 1.5,
                                  disk_gb: ceil(10.0 + input_size * 10.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    ##gzipped vcf from clean vcf part1.sh##
    int_vcf_gz=~{intermediate_vcf}

    ##Remove CNVs that are improperly genotyped by depth because they are nested within a real CNV##

    ##Determine columns of VCF after header##
    zcat $int_vcf_gz\
    |sed -n '1,1000p'\
    |egrep ^# \
    |tail -n 1 \
    |tr '\t' '\n' \
    |cat -n - \
    >col.txt

    ##Only affects CNV so pull those out##
    zcat $int_vcf_gz \
      |awk '{if ($5~"DEL" || $5~"DUP" || $1~"#") print}' \
      |svtk vcf2bed stdin tmp.bed
    awk -F"\t" '{if ($6=="") print $6="blanksample";print $0}' OFS='\t' tmp.bed \
      |gzip>int.bed.gz
  >>>

  output {
    File col="col.txt"
    File int_bed="int.bed.gz"
  }

}


task CleanVcf1b_2 {
  input {
    File int_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(int_bed, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 2.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    ##list of potenital overlaps with a normal copy state variant (>5kb variants require depth but nested events could be missed; i.e a duplication with a nest deletion will have a normal copy state for the deletion)##
    ##flip bed intersect so largest is CNV is always first##
    bedtools intersect -wa -wb -a <(zcat ~{int_bed}|awk '{if ($3-$2>=5000 ) print}') \
    -b <(zcat ~{int_bed}|awk '{if ($3-$2>=5000) print}') \
    |awk -F'\t' '{if ($4!=$10 && $3-$2>=$9-$8 && $5!=$11) print ;\
    else if ($4!=$10 && $5!=$11) print $7,$8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6}' OFS='\t' \
    |awk -F'\t' '{if ($6!="blanksample") print}' \
    |sort -u \
    >normaloverlap.txt
  >>>

  output {
    File normoverlap="normaloverlap.txt"
  }

}


task CleanVcf1b_3 {
  input {
    File int_vcf
    File normoverlap
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([int_vcf, normoverlap], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 7.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    ##pull out the depth based copy number variant for each normal overlapping variant##
    int_vcf_gz=~{int_vcf}
    cat <(zcat $int_vcf_gz|awk -F"\t" '{if ($1~"#") print}') \
        <(awk '{print $4 "\n" $10}' ~{normoverlap}|sort -u|fgrep -wf - <(zcat $int_vcf_gz)) \
      |awk '{if ($1!~"#") $1=$3;print}' OFS="\t" \
      |awk '{if ($1~"#" || $5=="<DEL>" || $5=="<DUP>") print}' \
      |vcftools --vcf - --stdout --extract-FORMAT-info RD_CN \
      |awk -F"\t" 'NR==1{for (i=3;i<=NF;i++) header[i]=$i} NR>1{for(j=3;j<=NF;j++) print $1"@"header[j] "\t" $j }' \
      |sort -k1,1 \
      |gzip \
      >RD_CN.normalcheck.FORMAT.gz
  >>>

  output {
    File rd_cn_normcheck="RD_CN.normalcheck.FORMAT.gz"
  }

}



task CleanVcf1b_4 {
  input {
    File int_vcf
    File normoverlap
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([int_vcf, normoverlap], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 7.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    ##pull out evidence supporting each normal overlapping variant##
    int_vcf_gz=~{int_vcf}
    cat <(zcat $int_vcf_gz|awk -F"\t" '{if ($1~"#") print}') \
      <(awk '{print $4 "\n" $10}' ~{normoverlap}|sort -u|fgrep -wf - <(zcat $int_vcf_gz)) \
      |awk '{if ($1!~"#") $1=$3;print}' OFS="\t"\
      |vcftools --vcf - --stdout --extract-FORMAT-info EV \
      |awk -F"\t" 'NR==1{for (i=3;i<=NF;i++) header[i]=$i} NR>1{for(j=3;j<=NF;j++) print $1"@"header[j] "\t" $j }' \
      |sort -k1,1 \
      |gzip \
      >EV.normalcheck.FORMAT.gz
  >>>

  output {
    File ev_normcheck="EV.normalcheck.FORMAT.gz"
  }

}


task CleanVcf1b_5 {
  input {
    File normoverlap
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(normoverlap, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 3.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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

    ##check if nested is incorrectly classified as normal##
    touch overlap.test.txt
    while read bed
    do
      echo $bed|tr ' ' '\t'|cut -f1-6 >large.bed
      echo $bed|tr ' ' '\t'|cut -f7-12>small.bed
      ##require at least 50% coverage to consider a variant overlapping##
      overlap=$(bedtools coverage -a small.bed -b large.bed|awk '{if ($NF>=0.50) print "YES";else print "NO"}')

      if [ "$overlap" == "YES" ]
      then
        smallid=$(awk '{print $4}' small.bed)

        ##pull out variants that are called a variants for both the smaller and larger CNVs (don't have normal copy state to check for)##
        if [ $(awk '{print $NF}' small.bed \
        |tr ',' '\n' \
        |fgrep -wvf - <(awk -F"[,\t]" -v var=$smallid '{for(i=6;i<=NF;i++) print var"@"$i "\t" $4"@"$i "\t" $5}' large.bed)|wc -l) -gt 0 ]
        then
          awk '{print $NF}' small.bed \
          |tr ',' '\n' \
          |fgrep -wvf - <(awk -F"[,\t]" -v var=$smallid '{for(i=6;i<=NF;i++) print var"@"$i "\t" $4"@"$i "\t" $5}' large.bed) \
          >>overlap.test.txt
        fi
      fi
    done<~{normoverlap}
  >>>

  output {
    File overlap_test="overlap.test.txt"
  }

}


task CleanVcf1b_6 {
  input {
    File overlap_test
    File rd_cn_normcheck
    File ev_normcheck
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([overlap_test, rd_cn_normcheck, ev_normcheck], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 4.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    ##determine variants that need to be revised from a normal copy state into a CNV##
    cat ~{overlap_test} \
    |sort -k1,1 \
    |join -j 1 - <(zcat ~{rd_cn_normcheck}) \
    |join -j 1 - <(zcat ~{ev_normcheck}) \
    |tr ' ' '\t' \
    |sort -k2,2 \
    |join -1 2 -2 1 - <(zcat ~{rd_cn_normcheck}) \
    |awk '{if ($3=="DUP" && $4==2 && $6==3) print $2 "\t" 1; else if ($3=="DEL" && $4==2 && $6==1)  print $2 "\t" 3 }' \
    |tr '@' '\t'\
    >geno.normal.revise.txt

  >>>

  output {
    File geno_normal_revise="geno.normal.revise.txt"
  }

}


task CleanVcf1b_7 {
  input {
    File int_vcf
    File geno_normal_revise
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([int_vcf, geno_normal_revise], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 3.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    ##Update genotypes##
    { zfgrep -wf <(awk '{print $1}' ~{geno_normal_revise}|sort -u) ~{int_vcf} || [[ $? == 1 ]]; }\
    |bgzip \
    >subset.vcf.gz
  >>>

  output {
    File subset_vcf="subset.vcf.gz"
  }

}


task CleanVcf1b_8 {
  input {
    File subset_vcf
    File geno_normal_revise
    File col
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([subset_vcf, geno_normal_revise, col], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 10.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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

    python3 <<CODE > normal.revise.vcf.lines.txt
    import gzip
    import sys

    VCF='~{subset_vcf}'
    REVISE='~{geno_normal_revise}'
    COL='~{col}'

    # Grab regenotyped samples of interest
    sys.stderr.write("Reading {}...\n".format(REVISE))
    geno_dict = {}
    with open(REVISE) as f:
      for line in f:
        tokens = line.strip().split('\t')
        vid = tokens[0]
        if vid not in geno_dict:
          geno_dict[vid] = []
        geno_dict[vid].append(tokens[1]) # id.txt but only sample id

    # Column definitions
    sys.stderr.write("Reading {}...\n".format(COL))
    sample_columns_dict = {}
    with open(COL) as f:
      for line in f:
        tokens = line.strip().split('\t')
        sample_columns_dict[tokens[1]] = int(tokens[0]) - 1

    # Assign GT/GQ
    sys.stderr.write("Reassigning genotypes...\n")
    with gzip.open(VCF, 'rb') as f:
      for lineb in f:
        line = lineb.decode('utf-8').strip()
        vid = line.split('\t', 3)[2]
        if vid in geno_dict:
          sample_ids = geno_dict[vid]
          tokens = line.split('\t')
          sample_indexes = [sample_columns_dict[s] for s in sample_ids]
          for i in sample_indexes:
            entry = tokens[i].split(':', 4)
            entry[0] = "0/1"
            entry[1] = entry[3]
            tokens[i] = ":".join(entry)
          sys.stdout.write("{}\t\n".format("\t".join(tokens)))
        else:
          sys.stdout.write("{}\t\n".format(line))
    CODE
  >>>

  output {
    File normal_revise_vcf_lines="normal.revise.vcf.lines.txt"
  }

}



task CleanVcf1b_9 {
  input {
    File int_vcf
    File normal_revise_vcf_lines
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([int_vcf, normal_revise_vcf_lines], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 15,
                                  disk_gb: ceil(10.0 + input_size * 50.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    ##rewrite vcf with updated genotypes##
    awk '{print $3}' ~{normal_revise_vcf_lines}|sort -u > vids.list
    cat <(zcat ~{int_vcf} | fgrep -wvf vids.list) \
      <(sed 's/\t$//' ~{normal_revise_vcf_lines}) \
      |vcf-sort \
      |bgzip \
      >normal.revise.vcf.gz

    bcftools index normal.revise.vcf.gz
  >>>

  output {
    File normal_revise_vcf="normal.revise.vcf.gz"
    File normal_revise_vcf_csi="normal.revise.vcf.gz.csi"
  }

}


task CleanVcf1b_10 {
  input {
    File normal_revise_vcf
    File normal_revise_vcf_csi
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(normal_revise_vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 15,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    ##get copy state per variant##
    zcat ~{normal_revise_vcf} \
      |awk '{if ($1!~"#") $1=$3;print}' OFS="\t" \
      |vcftools --vcf - --stdout --extract-FORMAT-info RD_CN \
      |gzip \
      >copystate.RD_CN.FORMAT.gz
  >>>

  output {
    File copystate_rd_cn="copystate.RD_CN.FORMAT.gz"
  }

}


task CleanVcf1b_11 {
  input {
    File copystate_rd_cn
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(copystate_rd_cn, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 15,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    ##get copy state per variant##
    zcat ~{copystate_rd_cn} \
      |awk 'NR>1{for(i=3;i<=NF;i++) lines[$1 "\t" $i]++ } END{for (x in lines) print x}' \
      |gzip \
      >copystate.per.variant.txt.gz
  >>>

  output {
    File copystate_per_variant="copystate.per.variant.txt.gz"
  }

}


task CleanVcf1b_12 {
  input {
    File int_bed
    File copystate_per_variant
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([int_bed, copystate_per_variant], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 15,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    ##Find multi-allelic for del or dup ; CNV >1kb we trust depth ##
    ##del##
    zcat ~{copystate_per_variant} \
      |awk '{if ($2!="." && $2>3) print $1}' \
      |sort -u \
      |{ fgrep -wf <(zcat ~{int_bed}|awk -F"\t" '{if ($5=="DEL" && $3-$2>=1000) print $4}' ) || [[ $? == 1 ]]; } \
      >multi.cnvs.del.txt
  >>>

  output {
    File multi_del="multi.cnvs.del.txt"
  }

}


task CleanVcf1b_13 {
  input {
    File int_bed
    File multi_del
    File copystate_per_variant
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([int_bed, multi_del, copystate_per_variant], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 15,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    ##dup##
    mv ~{multi_del} multi.cnvs.txt
    zcat ~{copystate_per_variant} \
      |awk '{if ($2!="." && ($2<1 || $2>4)) print $1}' \
      |sort -u \
      |{ fgrep -wf <(zcat ~{int_bed}|awk -F"\t" '{if ($5=="DUP" && $3-$2>=1000) print $4}' ) || [[ $? == 1 ]]; } \
      >>multi.cnvs.txt
  >>>

  output {
    File multi="multi.cnvs.txt"
  }

}
