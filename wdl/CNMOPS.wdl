version 1.0

##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/cnMOPS_hg38/1/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

import "Structs.wdl"

## cnMOPS worflow definition ##
## mode = sex (1=male, 2=female)

workflow CNMOPS {
  input {
    String r1
    String r2
    String batch
    File bincov_matrix
    File chrom_file
    File bincov_matrix_index
    File ped_file
    File blacklist
    File allo_file
    Array[String] samples
    String prefix
    Int? min_size
    Boolean? stitch_and_clean_large_events = false
    Float mem_gb_override_sample10 = 7.5
    Float mem_gb_override_sample3 = 16
    String linux_docker
    String sv_pipeline_docker
    String cnmops_docker
    RuntimeAttr? runtime_attr_sample10
    RuntimeAttr? runtime_attr_sample3
    RuntimeAttr? runtime_attr_ped
    RuntimeAttr? runtime_attr_clean
  }
  Array[Array[String]] Allos = read_tsv(allo_file)
  Array[Array[String]] Chroms = read_tsv(chrom_file)

  scatter (Allo in Allos) {
    call CNSampleNormal as MaleR2 {
      input:
        chr = Allo[0],
        black = blacklist,
        ped = ped_file,
        mode = "1",
        r = r2,
        bincov_matrix = bincov_matrix,
        bincov_matrix_index = bincov_matrix_index,
        mem_gb_override = mem_gb_override_sample10,
        cnmops_docker = cnmops_docker,
        runtime_attr_override = runtime_attr_sample10
    }

    call CNSampleNormal as MaleR1 {
      input:
        chr = Allo[0],
        black = blacklist,
        ped = ped_file,
        mode = "1",
        r = r1,
        bincov_matrix = bincov_matrix,
        bincov_matrix_index = bincov_matrix_index,
        mem_gb_override = mem_gb_override_sample3,
        cnmops_docker = cnmops_docker,
        runtime_attr_override = runtime_attr_sample3
    }
  }

  scatter (Chrom in Chroms) {
    call CNSampleNormal as NormalR2 {
      input:
        chr = Chrom[0],
        black = blacklist,
        ped = ped_file,
        mode = "normal",
        r = r2,
        bincov_matrix = bincov_matrix,
        bincov_matrix_index = bincov_matrix_index,
        mem_gb_override = mem_gb_override_sample10,
        cnmops_docker = cnmops_docker,
        runtime_attr_override = runtime_attr_sample10
    }

    call CNSampleNormal as NormalR1 {
      input:
        chr = Chrom[0],
        black = blacklist,
        ped = ped_file,
        mode = "normal",
        r = r1,
        bincov_matrix = bincov_matrix,
        bincov_matrix_index = bincov_matrix_index,
        mem_gb_override = mem_gb_override_sample3,
        cnmops_docker = cnmops_docker,
        runtime_attr_override = runtime_attr_sample3
    }
  }

  call CNSampleNormal as FemaleR2 {
    input:
      chr = "chrX",
      black = blacklist,
      ped = ped_file,
      mode = "2",
      r = r2,
      bincov_matrix = bincov_matrix,
      bincov_matrix_index = bincov_matrix_index,
      mem_gb_override = mem_gb_override_sample10,
      cnmops_docker = cnmops_docker,
      runtime_attr_override = runtime_attr_sample10
  }

  call CNSampleNormal as FemaleR1 {
    input:
      chr = "chrX",
      black = blacklist,
      ped = ped_file,
      mode = "2",
      r = r1,
      bincov_matrix = bincov_matrix,
      bincov_matrix_index = bincov_matrix_index,
      mem_gb_override = mem_gb_override_sample3,
      cnmops_docker = cnmops_docker,
      runtime_attr_override = runtime_attr_sample3
  }

  call GetPed {
    input:
      ped = ped_file,
      samples = samples,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_ped
  }

  call CleanCNMops {
    input:
      chrom_file=chrom_file,
      allo_file=allo_file,
      samplelist = GetPed.batchped,
      black = blacklist,
      batch = batch,
      NR1 = NormalR1.Gff,
      NR2 = NormalR2.Gff,
      MR1 = MaleR1.Gff,
      MR2 = MaleR2.Gff,
      FR1 = FemaleR1.Gff,
      FR2 = FemaleR2.Gff,
      prefix = prefix,
      min_size = min_size,
      stitch_and_clean_large_events = stitch_and_clean_large_events,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_ped
  }

  output {
    File Del = CleanCNMops.Del
    File Dup = CleanCNMops.Dup
    File Del_idx = CleanCNMops.Del_idx
    File Dup_idx = CleanCNMops.Dup_idx
  }
}

task GetPed {
  input {
    File ped
    Array[String] samples
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File batchped = "batch.ped"
  }
  command <<<
    egrep '~{sep="|" samples}' ~{ped} > batch.ped
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task CleanCNMops {
  input {
    File samplelist
    File black
    String batch
    Array[File] NR1
    Array[File] NR2
    Array[File] MR1
    Array[File] MR2
    File FR1
    File FR2
    String prefix
    Int? min_size = 1000000
    Boolean? stitch_and_clean_large_events = false
    File? chrom_file
    File? allo_file
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3, 
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File Del = "~{batch}.DEL.~{prefix}.bed.gz"
    File Del_idx = "~{batch}.DEL.~{prefix}.bed.gz.tbi"
    File Dup = "~{batch}.DUP.~{prefix}.bed.gz"
    File Dup_idx = "~{batch}.DUP.~{prefix}.bed.gz.tbi"
  }
  command <<<

    set -euo pipefail
    cut -f2 ~{samplelist} > sample.list
    cat ~{sep=" "  NR1} ~{sep=" "  NR2} ~{sep=" "  MR1} ~{sep=" "  MR2} ~{FR1} ~{FR2} > cnmops.gff

    mkdir calls
    grep -v "#" cnmops.gff > cnmops.gff1
    echo "./cnmops.gff1">GFF.list
    /opt/WGD/bin/cleancnMOPS.sh -z -o calls/ -S ~{black} sample.list GFF.list

    zcat calls/*/*.cnMOPS.DEL.bed.gz > DELS.bed 
    awk -v batch=~{batch}_DEL_ 'BEGIN{OFS="\t"} {print $1,$2,$3,batch,$4,"cnmops"}' DELS.bed | cat -n |\
    awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$5$1,$6,"DEL",$7}' | sort -k1,1V -k2,2n > ~{batch}.DEL.bed

    cat <(echo -e "#chr\tstart\tend\tname\tsample\tsvtype\tsources") ~{batch}.DEL.bed  > ~{batch}.DEL.~{prefix}.bed

    zcat calls/*/*.cnMOPS.DUP.bed.gz > DUPS.bed 
    awk -v batch=~{batch}_DUP_ 'BEGIN{OFS="\t"} {print $1,$2,$3,batch,$4,"cnmops"}' DUPS.bed | cat -n |\
    awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$5$1,$6,"DUP",$7}' | sort -k1,1V -k2,2n > ~{batch}.DUP.bed

    cat <(echo -e "#chr\tstart\tend\tname\tsample\tsvtype\tsources") ~{batch}.DUP.bed  > ~{batch}.DUP.~{prefix}.bed

    if ~{stitch_and_clean_large_events}; then
      mv ~{batch}.DUP.~{prefix}.bed ~{batch}.DUP.~{prefix}.prestitch.bed
      mv ~{batch}.DEL.~{prefix}.bed ~{batch}.DEL.~{prefix}.prestitch.bed
      cat ~{chrom_file} ~{allo_file} > contig.fai
      svtk rdtest2vcf --contigs contig.fai ~{batch}.DUP.~{prefix}.prestitch.bed sample.list dup.vcf.gz
      svtk rdtest2vcf --contigs contig.fai ~{batch}.DEL.~{prefix}.prestitch.bed sample.list del.vcf.gz
      bash /opt/sv-pipeline/04_variant_resolution/scripts/stitch_fragmented_CNVs.sh -d dup.vcf.gz dup1.vcf.gz
      bash /opt/sv-pipeline/04_variant_resolution/scripts/stitch_fragmented_CNVs.sh -d del.vcf.gz del1.vcf.gz
      svtk vcf2bed dup1.vcf.gz dup1.bed
      cat <(echo -e "#chr\tstart\tend\tname\tsample\tsvtype\tsources") \
          <(awk -v OFS="\t" -v minsize=~{min_size} '{if($3-$2>minsize)print $1,$2,$3,$4,$6,$5,"cnmops_large"}' dup1.bed) >~{batch}.DUP.~{prefix}.bed
      svtk vcf2bed del1.vcf.gz del1.bed
      cat <(echo -e "#chr\tstart\tend\tname\tsample\tsvtype\tsources") \
          <(awk -v OFS="\t" -v minsize=~{min_size} '{if($3-$2>minsize)print $1,$2,$3,$4,$6,$5,"cnmops_large"}' del1.bed) >~{batch}.DEL.~{prefix}.bed
    fi
    bgzip -f ~{batch}.DEL.~{prefix}.bed
    tabix -f ~{batch}.DEL.~{prefix}.bed.gz

    bgzip -f ~{batch}.DUP.~{prefix}.bed
    tabix -f ~{batch}.DUP.~{prefix}.bed.gz
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

task CNSampleNormal {
  input {
    String chr
    File black
    File ped
    String mode
    String r
    File bincov_matrix
    File bincov_matrix_index
    Float? mem_gb_override
    String cnmops_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    bincov_matrix: {
      localization_optional: true
    }
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 16,
    disk_gb: 20,
    boot_disk_gb: 10,
    preemptible_tries: 3, 
    max_retries : 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File Gff = "calls/cnMOPS.cnMOPS.gff"
  }
  command <<<

    set -euo pipefail
    GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` \
      tabix -h ~{bincov_matrix} "~{chr}" | sed 's/Start/start/' | sed 's/Chr/chr/' | sed 's/End/end/' > bincov_~{chr}.bed

    if [ ~{mode} == "normal" ]; then  
      mv bincov_~{chr}.bed bincov_~{chr}_~{mode}.bed
    else 
      awk -v sex="~{mode}" '$5==sex' ~{ped} | cut -f2 > ids.to.include
      col=$(head -n 1 bincov_~{chr}.bed | tr '\t' '\n'|cat -n| grep -wf ids.to.include | awk -v ORS="," '{print $1}' | sed 's/,$//g' | sed 's:\([0-9]\+\):$&:g')
      col_a="{print \$1,\$2,\$3,$col}"
      awk -f <(echo "$col_a") bincov_~{chr}.bed | tr ' ' '\t' > bincov_~{chr}_~{mode}.bed
    fi

    EMPTY_OUTPUT_ERROR="No CNV regions in result object. Rerun cn.mops with different parameters!"
    set +e
    bash /opt/WGD/bin/cnMOPS_workflow.sh -S ~{black} -x ~{black} -r ~{r} -o . -M bincov_~{chr}_~{mode}.bed &> cnmops.out
    RC=$?
    set -e
    if [ ! $RC -eq 0 ]; then
      if grep -q "$EMPTY_OUTPUT_ERROR" "cnmops.out"; then
        touch calls/cnMOPS.cnMOPS.gff
      else
        echo "cnMOPS_workflow.sh returned a non-zero code that was not due to an empty call file."
        exit $RC
      fi
    fi
    
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([mem_gb_override, runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: cnmops_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
