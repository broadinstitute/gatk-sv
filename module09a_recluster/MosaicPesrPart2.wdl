version 1.0

import "Structs.wdl"

workflow Mosaic{
  input{
    String name
    Int rare_cutoff
    File outlier
    File lookup
    File potential
    File coverage_file
    File coverage_file_idx
    File fam_file
    File median_file
    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker
  }
  call GetPotential{
    input:
        outlier=outlier,
        rare_cutoff=rare_cutoff,
        potential=potential,
        name=name,
        lookup=lookup,
        sv_pipeline_docker=sv_pipeline_docker
  }
  call rdtest{
    input:
      bed=GetPotential.rare,
      coverage_file=coverage_file,
      coverage_file_idx=coverage_file_idx,
      median_file=median_file,
      fam_file=fam_file,
      prefix=name,
      sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker
  }
  output{
    File potentialmosaic = GetPotential.rare
    File igvplots = rdtest.plots
  }
}

task GetPotential{
  input{
    String name
    Int rare_cutoff
    File outlier
    File potential
    File lookup
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 8,
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command<<<
    set -euox pipefail
    fgrep -v -f ~{outlier} ~{potential} >potential.txt
    while read chr start end id type sample;do
      n=$(zfgrep "$id:" ~{lookup}|cut -f 5) ||true
      if [ "$n" -eq "$n" ] ;then
        if [ "$n" -lt ~{rare_cutoff} ]; then
          printf "$chr\t$start\t$end\t$id\t$type\t$sample\n"
        fi
      fi
    done<potential.txt > ~{name}.potentialmosaic.rare.bed

    echo -e "#chr\tstart\tend\tid\ttype\tsample" > header.bed
    cat header.bed ~{name}.potentialmosaic.rare.bed | bgzip > ~{name}.potentialmosaic.rare.bed.gz
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
  output{
    File rare="~{name}.potentialmosaic.rare.bed.gz"
  }
}
# Run rdtest
task rdtest {
  input{
    File bed
    String coverage_file
    File coverage_file_idx
    File median_file
    File fam_file
    String prefix
    String sv_pipeline_rdtest_docker
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
  command <<<
    set -euo pipefail

    zcat ~{bed} | tail -n+2 > rdtest.bed
    /opt/RdTest/localize_bincov.sh rdtest.bed ~{coverage_file}
    awk -v OFS="\t" '{print $1,$2,$3,$4,$6,$5}' rdtest.bed > test.bed

    Rscript /opt/RdTest/RdTest.R \
      -b test.bed \
      -n ~{prefix} \
      -c local_coverage.bed.gz \
      -m ~{median_file} \
      -f ~{fam_file} \
      -p TRUE 
    mkdir plots
    mv *jpg plots
    tar -czvf mosaic.tar.gz plots/
  >>>
  
  output {
    File stats = "~{prefix}.metrics"
    File local_coverage = "local_coverage.bed.gz"
    File plots= "mosaic.tar.gz"
  }  
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_rdtest_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
