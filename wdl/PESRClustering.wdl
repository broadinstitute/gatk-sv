version 1.0

import "Structs.wdl"

workflow ClusterPESR {
  input {
    Array[File] vcfs
    String batch
    Float frac
    String svtypes
    Int svsize
    String algorithm
    File exclude_list
    File contigs
    Int dist
    String flags

    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_cluster
    RuntimeAttr? runtime_attr_concat
  }

  Array[Array[String]] contiglist = read_tsv(contigs)

  scatter (contig in contiglist) {
    call VCFCluster {
      input:
        vcfs = vcfs,
        batch = batch,
        algorithm = algorithm,
        chrom = contig[0],
        dist = dist,
        frac = frac,
        exclude_list = exclude_list,
        svsize = svsize,
        svtypes = svtypes,
        flags = flags,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_cluster
    }
  }
  call ConcatVCFs {
    input:
      vcfs = VCFCluster.clustered_vcf,
      batch = batch,
      algorithm = algorithm,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat
  }

  output {
    File clustered_vcf = ConcatVCFs.vcf
    File clustered_vcf_index = ConcatVCFs.vcf_index
  }
}

task VCFCluster {
  input {
    Array[File] vcfs
    String batch
    String algorithm
    String chrom
    Int dist
    Float frac
    File exclude_list
    Int svsize
    String svtypes
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
    max_retries: 3
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File clustered_vcf = "${batch}.${algorithm}.${chrom}.vcf.gz"
  }
  command <<<

    set -euo pipefail
    for f in ~{sep=" "  vcfs}; do tabix -p vcf -f $f; done;
    tabix -p bed ~{exclude_list};

    svtk vcfcluster ~{write_lines(vcfs)} stdout \
      -r ~{chrom} \
      -p ~{batch}_~{algorithm}_~{chrom} \
      -d ~{dist} \
      -f ~{frac} \
      -x ~{exclude_list} \
      -z ~{svsize} \
      -t ~{svtypes} \
      ~{flags} \
      | vcf-sort -c \
      | bgzip -c > ~{batch}.~{algorithm}.~{chrom}.vcf.gz
  
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

task ConcatVCFs {
  input {
    Array[File] vcfs
    String batch
    String algorithm
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
    File vcf = "${batch}.${algorithm}.vcf.gz"
    File vcf_index = "${batch}.${algorithm}.vcf.gz.tbi"
  }
  command <<<

    set -euo pipefail
    vcf-concat ~{sep=" "  vcfs} | vcf-sort -c | bgzip -c > ~{batch}.~{algorithm}.vcf.gz;
    tabix -p vcf ~{batch}.~{algorithm}.vcf.gz;
  
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

