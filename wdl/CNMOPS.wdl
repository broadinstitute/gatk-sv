version 1.0

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
    File exclude_list
    File allo_file
    File ref_dict
    Array[String] samples
    String prefix
    Int? min_size
    Boolean? stitch_and_clean_large_events = false
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
        exclude = exclude_list,
        ped = ped_file,
        mode = "1",
        r = r2,
        ref_dict = ref_dict,
        bincov_matrix = bincov_matrix,
        bincov_matrix_index = bincov_matrix_index,
        cnmops_docker = cnmops_docker,
        runtime_attr_override = runtime_attr_sample10
    }

    call CNSampleNormal as MaleR1 {
      input:
        chr = Allo[0],
        exclude = exclude_list,
        ped = ped_file,
        mode = "1",
        r = r1,
        ref_dict = ref_dict,
        bincov_matrix = bincov_matrix,
        bincov_matrix_index = bincov_matrix_index,
        cnmops_docker = cnmops_docker,
        runtime_attr_override = runtime_attr_sample3
    }
  }

  scatter (Chrom in Chroms) {
    call CNSampleNormal as NormalR2 {
      input:
        chr = Chrom[0],
        exclude = exclude_list,
        ped = ped_file,
        mode = "normal",
        r = r2,
        ref_dict = ref_dict,
        bincov_matrix = bincov_matrix,
        bincov_matrix_index = bincov_matrix_index,
        cnmops_docker = cnmops_docker,
        runtime_attr_override = runtime_attr_sample10
    }

    call CNSampleNormal as NormalR1 {
      input:
        chr = Chrom[0],
        exclude = exclude_list,
        ped = ped_file,
        mode = "normal",
        r = r1,
        ref_dict = ref_dict,
        bincov_matrix = bincov_matrix,
        bincov_matrix_index = bincov_matrix_index,
        cnmops_docker = cnmops_docker,
        runtime_attr_override = runtime_attr_sample3
    }
  }

  call CNSampleNormal as FemaleR2 {
    input:
      chr = "chrX",
      exclude = exclude_list,
      ped = ped_file,
      mode = "2",
      r = r2,
      ref_dict = ref_dict,
      bincov_matrix = bincov_matrix,
      bincov_matrix_index = bincov_matrix_index,
      cnmops_docker = cnmops_docker,
      runtime_attr_override = runtime_attr_sample10
  }

  call CNSampleNormal as FemaleR1 {
    input:
      chr = "chrX",
      exclude = exclude_list,
      ped = ped_file,
      mode = "2",
      r = r1,
      ref_dict = ref_dict,
      bincov_matrix = bincov_matrix,
      bincov_matrix_index = bincov_matrix_index,
      cnmops_docker = cnmops_docker,
      runtime_attr_override = runtime_attr_sample3
  }

  call CleanCNMops {
    input:
      chrom_file=chrom_file,
      allo_file=allo_file,
      samplelist = ped_file,
      exclude = exclude_list,
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

task CleanCNMops {
  input {
    File samplelist
    File exclude
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
    /opt/WGD/bin/cleancnMOPS.sh -z -o calls/ -S ~{exclude} sample.list GFF.list

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
    File exclude
    File ped
    String mode
    String r
    File ref_dict
    File bincov_matrix
    File bincov_matrix_index
    String cnmops_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    bincov_matrix: {
      localization_optional: true
    }
  }

  Float mem_gb_base = 4.0
  Float mem_gb_scale = 4.0
  Float mem_gb = mem_gb_base + mem_gb_scale * size(bincov_matrix, "GiB")
  Int disk_gb_base = 10
  Float disk_gb_scale = 2.0
  Int disk_gb = disk_gb_base + ceil(size([exclude, ped, bincov_matrix_index], "GiB") + disk_gb_scale * size(bincov_matrix, "GiB"))
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: mem_gb,
    disk_gb: disk_gb,
    boot_disk_gb: 10,
    preemptible_tries: 3, 
    max_retries : 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb_used = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb_used * 1000 * 0.8)

  output {
    File Gff = "calls/cnMOPS.cnMOPS.gff"
  }
  command <<<

    set -euo pipefail

    java -Xmx~{java_mem_mb}M -jar ${GATK_JAR} PrintSVEvidence \
      --sequence-dictionary ~{ref_dict} \
      --evidence-file ~{bincov_matrix} \
      -L ~{chr} \
      -O ~{chr}.RD.txt

    if [ ~{mode} == "normal" ]; then  
      mv ~{chr}.RD.txt ~{chr}.~{mode}.RD.txt
    else 
      awk -v sex="~{mode}" '$5==sex' ~{ped} | cut -f2 > ids.to.include
      col=$(head -n 1 ~{chr}.RD.txt | tr '\t' '\n'|cat -n| grep -wf ids.to.include | awk -v ORS="," '{print $1}' | sed 's/,$//g' | sed 's:\([0-9]\+\):$&:g')
      col_a="{print \$1,\$2,\$3,$col}"
      awk -f <(echo "$col_a") ~{chr}.RD.txt | tr ' ' '\t' > ~{chr}.~{mode}.RD.txt
    fi

    # redirect stdout and stderr to cnmops.out so that EMPTY_OUTPUT_ERROR can be detected, but use tee to also output them to
    # terminal so that errors can be debugged
    EMPTY_OUTPUT_ERROR="No CNV regions in result object. Rerun cn.mops with different parameters!"
    set +e
    bash /opt/WGD/bin/cnMOPS_workflow.sh -S ~{exclude} -x ~{exclude} -r ~{r} -o . -M ~{chr}.~{mode}.RD.txt 2>&1 | tee cnmops.out
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
    memory: mem_gb_used + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: cnmops_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
