version 1.0

import "Structs.wdl"

workflow MakeBincovMatrix {
  input {
    Array[String] samples
    Array[File] count_files
    Array[String]? bincov_matrix_samples
    File? bincov_matrix
    File reference_dict
    String batch
    Int? binsize
    Boolean skip_bin_size_filter = false
    String gatk_docker
    String sv_base_mini_docker
    String sv_base_docker
    RuntimeAttr? runtime_attr_override
  }

  Int output_binsize = select_first([binsize, 100])
  Array[String] all_samples = flatten([samples, select_first([bincov_matrix_samples, []])])
  Array[File] all_count_files = flatten([count_files, select_all([bincov_matrix])])

  scatter (i in range(length(all_count_files))) {
    call PrepareBincovEvidenceFile {
      input:
        count_file = all_count_files[i],
        sample_id = samples[i],
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_override
    }
  }

  call BuildBincovMatrix {
    input:
      evidence_files = PrepareBincovEvidenceFile.rd_file,
      evidence_file_indices = PrepareBincovEvidenceFile.rd_file_index,
      samples = all_samples,
      reference_dict = reference_dict,
      batch = batch,
      binsize = output_binsize,
      skip_bin_size_filter = skip_bin_size_filter,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File merged_bincov = BuildBincovMatrix.merged_bincov
    File merged_bincov_idx = BuildBincovMatrix.merged_bincov_idx
  }
}

task PrepareBincovEvidenceFile {
  input {
    File count_file
    String sample_id
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  Float file_size = size(count_file, "GiB")
  Int disk_gb = 50 + ceil(file_size * 3.0)
  Int java_heap_size_mb = 1024
  Float mem_gb = java_heap_size_mb / 1024.0 + 1.5
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: mem_gb,
    disk_gb: disk_gb,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File rd_file = "~{sample_id}.rd.txt.gz"
    File rd_file_index = "~{sample_id}.rd.txt.gz.tbi"
  }

  command <<<
    set -euo pipefail
    
    ln -s "~{count_file}" "sample.counts.tsv.gz"
    /gatk/gatk --java-options "-Xmx~{java_heap_size_mb}m" IndexFeatureFile -I "sample.counts.tsv.gz"

    /gatk/gatk --java-options "-Xmx~{java_heap_size_mb}m" ConvertCountsToDepthFile \
      --sample-name "~{sample_id}" \
      --counts-filename "sample.counts.tsv.gz" \
      --output "~{sample_id}.rd.txt.gz"
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }
}

task BuildBincovMatrix {
  input {
    Array[File] evidence_files
    Array[File] evidence_file_indices
    Array[String] samples
    File reference_dict
    String batch
    Int binsize
    Boolean skip_bin_size_filter = false
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  Float file_size = size(evidence_files, "GiB")
  Int disk_gb = 100 + ceil(file_size * 4.0)
  Int java_heap_size_mb = round(42.0 * length(evidence_files) + 1024.0)
  Float mem_gb = java_heap_size_mb / 1024.0 + 2.5
  RuntimeAttr default_attr = object {
    cpu_cores: 2,
    mem_gb: mem_gb,
    disk_gb: disk_gb,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String merged_bincov_pre_filter = "~{batch}_pre_filter.RD.txt.gz"
  String merged_bincov_file_name = "~{batch}.RD.txt.gz"

  output {
    File merged_bincov = merged_bincov_file_name
    File merged_bincov_idx = "~{merged_bincov_file_name}.tbi"
  }

  command <<<
    set -euo pipefail
    ulimit -n 100000

    cp ~{write_lines(evidence_files)} evidence.list
    cp ~{write_lines(samples)} samples.list

    /gatk/gatk --java-options "-Xmx~{java_heap_size_mb}m" PrintSVEvidence \
      -F evidence.list \
      --sample-names samples.list \
      --sequence-dictionary ~{reference_dict} \
      --output "~{merged_bincov_pre_filter}"

    tabix -f -p bed "~{merged_bincov_pre_filter}"

    if [[ "~{skip_bin_size_filter}" == "false" ]]; then
      zcat "~{merged_bincov_pre_filter}" \
        | awk -v bin_size="~{binsize}" 'NR==1 || ($3 - $2) == bin_size' \
        | bgzip -c > "~{merged_bincov_file_name}"

      tabix -f -p bed "~{merged_bincov_file_name}"
    else
      mv "~{merged_bincov_pre_filter}" "~{merged_bincov_file_name}"
      mv "~{merged_bincov_pre_filter}.tbi" "~{merged_bincov_file_name}.tbi"
    fi
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }
}
