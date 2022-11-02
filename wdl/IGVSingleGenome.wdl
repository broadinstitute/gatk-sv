version 1.0

import "Structs.wdl"

workflow IGVSingleGenome {
  input{
    # Bed file containing regions to screenshot; 4th column must be SVID
    File bed

    # Sample reads
    File bam_or_cram
    File bam_or_cram_index

    # Sample id and prefix for output filenames
    String sample_id
    String run_name

    # Reference corresponding to read alignments
    File ref_fasta
    File ref_fai

    Int? records_per_shard

    String linux_docker
    String igv_docker

    RuntimeAttr? runtime_attr_scatter
    RuntimeAttr? runtime_attr_igv
  }

  Int records_per_shard_ = select_first([records_per_shard, 300])

  call ScatterBed {
    input:
      bed=bed,
      size=records_per_shard_,
      prefix="~{run_name}-~{sample_id}.shard_",
      linux_docker=linux_docker,
      runtime_attr_override=runtime_attr_scatter
  }

  scatter (i in range(length(ScatterBed.out))) {
    call RunIGV {
      input:
        bed=ScatterBed.out[i],
        ref_fasta=ref_fasta,
        ref_fai=ref_fai,
        bam_or_cram=bam_or_cram,
        bam_or_cram_index=bam_or_cram_index,
        prefix="~{run_name}-~{sample_id}-shard-~{i}",
        igv_docker=igv_docker,
        runtime_attr_override=runtime_attr_igv
    }
  }

  output{
    Array[File] igv_plots_tarballs = RunIGV.out
    Array[File] igv_scripts = RunIGV.script
  }
}

task ScatterBed {
  input {
    File bed
    Int size
    String prefix
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 10,
                               disk_gb: ceil(10 + size(bed, "GB") * 2),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    mkdir out
    awk '0!~"#"' ~{bed} | split -a 5 -d -l ~{size} - out/~{prefix}
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

  output{
    Array[File] out = glob("out/~{prefix}*")
  }
}

task RunIGV {
  input {
    File bed
    File ref_fasta
    File ref_fai
    File bam_or_cram
    File bam_or_cram_index
    String prefix
    String igv_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 10,
                               disk_gb: ceil(100 + size([bed, ref_fasta, bam_or_cram], "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    python3 /src/make_igv_script.py \
      --bed ~{bed} \
      --reads ~{bam_or_cram} \
      --snapshot-dir "~{prefix}-screenshots" \
      --snapshot-prefix "~{prefix}" \
      --out "~{prefix}-screenshot-script.txt" \
      --reference ~{ref_fasta}
    xvfb-run --server-args="-screen 0, 1920x3000x24" bash /IGV_2.4.14/igv.sh -b "~{prefix}-screenshot-script.txt"
    tar -czf "~{prefix}-screenshots.tar.gz" "~{prefix}-screenshots"
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: igv_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  output{
    File out = "~{prefix}-screenshots.tar.gz"
    File script = "~{prefix}-screenshot-script.txt"
  }
}
