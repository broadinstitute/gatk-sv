version 1.0

import "Structs.wdl"

workflow GenotypeGenomicDisorderRegions {
  input {
    Array[String] batch_names
    File genomic_disorder_regions_bed
    Array[File] rd_files
    Array[File] median_files
    Array[File] depth_sepcutoff_files
    String sv_pipeline_docker
    RuntimeAttr? runtime_generate_median_geno
  }
  scatter (i in range(length(batch_names))) {
    call RunRdTest {
      input:
        batch_name = batch_names[i],
        rdtest_bed = rd_files[i],
        rd_file = rd_files[i],
        rd_index = rd_files[i] + ".tbi",
        median_file = median_files[i],
        depth_sepcutoff = median_files[i],
        sv_pipeline_rdtest_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_generate_median_geno
    }
  }
  output{
    Array[File] rdtest_out = RunRdTest.out
  }
}

task RunRdTest {
  input{
    String batch_name
    File rdtest_bed
    File rd_file
    File rd_index
    File median_file
    File depth_sepcutoff
    Int large_size_cutoff = 1000000
    String sv_pipeline_rdtest_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 30,
                               disk_gb: ceil(40.0 + size(rd_file, "GiB") * 4),
                               boot_disk_gb: 30,
                               preemptible_tries: 1,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    mkdir rdtest_~{batch_name}/
    Rscript /opt/RdTest/RdTest.R \
      -v TRUE -g TRUE -p TRUE \
      -r ~{depth_sepcutoff} \
      -b ~{rdtest_bed} \
      -c ~{rd_file} \
      -m ~{median_file} \
      -n ~{batch_name} \
      -s ~{large_size_cutoff} -o rdtest_~{batch_name}
    tar czvf rdtest_~{batch_name}.tar.gz rdtest_~{batch_name}/
  >>>
  output{
    File out = "rdtest_~{batch_name}.tar.gz"
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
