version 1.0

import "Structs.wdl"

workflow ClusterDepth {
  input {
    File contigs
    Float frac
    File del_bed
    String flags
    File dup_bed
    String batch
    File? exclude_list
    Float? exclude_list_frac_max


    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_bed_cluster
    RuntimeAttr? runtime_attr_concat
    RuntimeAttr? runtime_attr_depth_vcf
    RuntimeAttr? runtime_attr_rdtest_bed
  }

  Array[Array[String]] contiglist = read_tsv(contigs)

  scatter (contig in contiglist) {
    call BedCluster as BedCluster_del {
      input:
        batch = batch,
        svtype = "DEL",
        chrom = contig[0],
        bed = del_bed,
        frac = frac,
        exclude_list=exclude_list,
        exclude_list_frac_max=exclude_list_frac_max,
        flags = flags,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_bed_cluster
    }

    call BedCluster as BedCluster_dup {
      input:
        batch = batch,
        svtype = "DUP",
        chrom = contig[0],
        bed = dup_bed,
        frac = frac,
        exclude_list=exclude_list,
        exclude_list_frac_max=exclude_list_frac_max,
        flags = flags,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_bed_cluster
    }
  }

  call ConcatBeds as ConcatBeds_del {
    input:
      batch = batch,
      svtype = "DEL",
      beds = BedCluster_del.clustered_bed,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat
  }

  call ConcatBeds as ConcatBeds_dup {
    input:
      batch = batch,
      svtype = "DUP",
      beds = BedCluster_dup.clustered_bed,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat
  }

  call MakeRDTestBed {
    input:
      dels = ConcatBeds_del.merged_bed,
      dups = ConcatBeds_dup.merged_bed,
      batch = batch,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_rdtest_bed
  }

  call MakeDepthVCF {
    input:
      bed = MakeRDTestBed.bed,
      contigs = contigs,
      batch = batch,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_depth_vcf
  }

  output {
    File clustered_vcf = MakeDepthVCF.vcf
    File clustered_vcf_index = MakeDepthVCF.vcf_index
  }
}

task MakeRDTestBed {
  input {
    File dels
    File dups
    String batch
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
    File bed = "${batch}.depth.bed"
  }
  command <<<

    set -euo pipefail
    python3 /opt/sv-pipeline/scripts/make_depth_rdtest_bed.py ~{dels} | sed '1d' > del.bed
    python3 /opt/sv-pipeline/scripts/make_depth_rdtest_bed.py ~{dups} | sed '1d' > dup.bed
    echo -e "#chrom start end name samples svtype" | sed -e 's/ /\t/g' > ~{batch}.depth.bed
    cat del.bed dup.bed | sort -k1,1V -k2,2n >> ~{batch}.depth.bed
  
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

task MakeDepthVCF {
  input {
    File bed
    File contigs
    String batch
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
    File vcf = "${batch}.depth.vcf.gz"
    File vcf_index = "${batch}.depth.vcf.gz.tbi"
  }
  command <<<

    set -euo pipefail
    cut -f5 ~{bed} | sed -e '1d' -e 's/,/\n/g' | sort -u > samples.list
    svtk rdtest2vcf --contigs ~{contigs} ~{bed} samples.list ~{batch}.depth.vcf.gz
  
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

task BedCluster {
  input {
    String batch
    String svtype
    String chrom
    File bed
    File? exclude_list
    Float? exclude_list_frac_max = 0.5
    Float frac
    String flags
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
    File clustered_bed = "${batch}.${svtype}.${chrom}.bed"
  }
  command <<<

    set -euo pipefail
    tabix -p bed ~{bed};
    svtk bedcluster ~{bed} -r ~{chrom} \
      -p ~{batch}_depth_~{svtype}_~{chrom} \
      -f ~{frac} \
      ~{flags} \
      > ~{batch}.~{svtype}.~{chrom}.preexcludelist.bed

    ~{if defined(exclude_list) then
      "bedtools coverage -a ~{batch}.~{svtype}.~{chrom}.preexcludelist.bed -b ~{exclude_list} | awk '$NF < ~{exclude_list_frac_max}' | rev | cut -f5- | rev > excluded.filtered.bed"
      else
      ""}
    ~{if defined(exclude_list) then
       "cat <(head -1 ~{batch}.~{svtype}.~{chrom}.preexcludelist.bed) excluded.filtered.bed > ~{batch}.~{svtype}.~{chrom}.bed"
      else
       "mv ~{batch}.~{svtype}.~{chrom}.preexcludelist.bed ~{batch}.~{svtype}.~{chrom}.bed"}
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

task ConcatBeds {
  input {
    String batch
    String svtype
    Array[File] beds
    String sv_base_mini_docker
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
    File merged_bed = "${batch}.${svtype}.bed"
  }
  command <<<

    awk 'FNR==1 && NR!=1 { while (/^#chrom/) getline; } 1 {print}' ~{sep=" "  beds} > ~{batch}.~{svtype}.bed

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

