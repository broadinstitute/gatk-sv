version 1.0

import "Structs.wdl"

workflow VCFToBEDWorkflow {
    input {
        Array[File] vcfs
        Array[File] idxs
        String output_prefix
        String sv_base_mini_docker
        String sv_base_pipeline_docker
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_concat_bed
    }

    scatter (i in range(length(vcfs))) {
        call VCFToBED {
            input:
                vcf = vcfs[i],
                index = idxs[i],
                sv_base_pipeline_docker = sv_base_pipeline_docker,
                runtime_attr_override = runtime_attr_vcf2bed
        }
    }

    call ConcatBeds as concat_beds{
        input:
            shard_bed_files = VCFToBED.bed,
            prefix = output_prefix,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat_bed
    }

    output {
        File output_bed = concat_beds.merged_bed_file
    }

}


task VCFToBED {
    input {
        File vcf
        File index
        String sv_base_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 10.0,
        disk_gb: ceil(10.0 + size(vcf, "GiB") * 2.0),
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10
    }   

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    String prefix = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail
        svtk vcf2bed -i ALL --include-filters ~{vcf} ~{prefix}.bed
    >>>

    output {
        File bed = "~{prefix}.bed"
    }

    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_base_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
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

  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }
}


