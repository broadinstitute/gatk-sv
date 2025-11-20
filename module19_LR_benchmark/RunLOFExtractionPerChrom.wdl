version 1.0

workflow RunLOFExtractionPerChrom {
    input {
        Array[File] vcfs
        Array[File] vcf_indexes   # matching 1-to-1 with vcfs
        File extract_script       # python extract_SV_gene.lof.py
        String prefix
        String sv_pipeline_base_docker
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_extract_lof
        RuntimeAttr? runtime_attr_concat_bed
    }

    scatter (i in range(length(vcfs))) {
        call ExtractLOF {
            input:
                vcf = vcfs[i],
                index = vcf_indexes[i],
                extract_script = extract_script,
                sv_pipeline_base_docker = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_extract_lof
        }
    }

    call ConcatBeds{
        input:
            shard_bed_files = ExtractLOF.bed_file,
            prefix  = prefix,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat_bed
    }

    output {
        File bed_file = ConcatBeds.merged_bed_file
    }
}

task ExtractLOF {
    input {
        File vcf
        File index
        File extract_script
        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 5,
        disk_gb: 10+ceil(size(vcf, "GiB") *2),
        boot_disk_gb: 10,
        preemptible_tries: 0,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(vcf, ".vcf.bgz")

    command <<<
        set -euo pipefail

        python3 ~{extract_script} ~{vcf} ~{prefix}.bed
        bgzip "~{prefix}.bed"
    >>>

    output {
        File bed_file = "~{prefix}.bed.gz"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_base_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}


task ConcatBeds {
  input {
    Array[File] shard_bed_files
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }
  
  Float input_size = size(shard_bed_files, "GB")

  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(10.0 + input_size * 3.0),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }   

  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eux

    # note head -n1 stops reading early and sends SIGPIPE to zcat,
    # so setting pipefail here would result in early termination
    zcat ~{shard_bed_files[0]} | head -n1 > header.txt

    # no more early stopping
    set -o pipefail

    while read SPLIT; do
      zcat $SPLIT | tail -n+2
    done < ~{write_lines(shard_bed_files)} \
      | cat header.txt - \
      | bgzip -c \
      > "~{prefix}.bed.gz"

   >>>

  output {
    File merged_bed_file = "~{prefix}.bed.gz"
  }
}



