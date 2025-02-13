version 1.0

import "Structs.wdl"
import "TasksBenchmark.wdl" as mini_tasks


workflow AnnotateILFeatures{
    input{
        Array[File] query_file_list
        Array[File] ref_file_list
        Array[String] contig_list
        File compare_callset_sh
        File compare_callse_helper_R
        String prefix

        String rdpesr_benchmark_docker
        String sv_base_mini_docker
    }

    scatter(i in range(length(contig_list))){
        call BedComparisons{
            input:
                query_file = query_file_list[i],
                ref_file = ref_file_list[i],
                compare_callset_sh = compare_callset_sh,
                compare_callse_helper_R = compare_callse_helper_R,
                prefix = "~{prefix}.~{contig_list[i]}",
                rdpesr_benchmark_docker=rdpesr_benchmark_docker
        }
    }

    call ConcatBeds{
        input:
            shard_bed_files = BedComparisons.vs_manta,
            prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker
    }

    output{
        File comparison_results = ConcatBeds.merged_bed_file
    }
}



task BedComparisons{
    input{
        File query_file
        File ref_file
        String prefix
        File compare_callset_sh
        File compare_callse_helper_R
        String rdpesr_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File vs_manta = "~{prefix}.bed.gz"
    }

    command <<<
        bash ~{compare_callset_sh} -r ~{compare_callse_helper_R} -O ~{prefix}.bed -p ~{prefix} ~{query_file} ~{ref_file}
        bgzip ~{prefix}.bed
    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: rdpesr_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}


task ConcatBeds {
  input {
    Array[File] shard_bed_files
    String prefix
    Boolean? index_output
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Boolean call_tabix = select_first([index_output, false])
  String output_file="~{prefix}.bed.gz"

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(shard_bed_files, "GB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + compression_factor)),
    cpu_cores: 1,
    preemptible_tries: 3,
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
      zcat $SPLIT
    done < ~{write_lines(shard_bed_files)} \
      | (grep -Ev "^#" || printf "") \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | cat header.txt - \
      | bgzip -c \
      > ~{output_file}

    if ~{call_tabix}; then
      tabix -p bed ~{output_file}
    else
      touch ~{output_file}.tbi
    fi
  >>>

  output {
    File merged_bed_file = output_file
    File merged_bed_idx = output_file + ".tbi"
  }
  }


