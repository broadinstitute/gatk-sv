version 1.0

import "Structs.wdl"

workflow BatchEvidenceMerging {
  input {
    Array[String] samples
    Array[File?]? BAF_files
    Array[File] PE_files
    Array[File] SR_files
    Array[File]? SD_files
    File? sd_locs_vcf
    File reference_dict
    String batch
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  call MergeEvidence as MergeSREvidence {
    input:
      files = SR_files,
      batch = batch,
      evidence = "sr",
      samples = samples,
      reference_dict = reference_dict,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_override
  }

  call MergeEvidence as MergePEEvidence {
    input:
      files = PE_files,
      batch = batch,
      evidence = "pe",
      samples = samples,
      reference_dict = reference_dict,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_override
  }

  if (defined(BAF_files)) {
    call MergeEvidence as MergeBAFEvidence {
      input:
        files = select_all(select_first([BAF_files])),
        batch = batch,
        evidence = "baf",
        samples = samples,
        reference_dict = reference_dict,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_override
    }
  }
  if (!defined(BAF_files)) {
    call SDtoBAF {
      input:
        SD_files = select_first([SD_files]),
        sd_locs_vcf = select_first([sd_locs_vcf]),
        batch = batch,
        samples = samples,
        reference_dict = reference_dict,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_override
    }
  }

  output {
    File merged_BAF = select_first([MergeBAFEvidence.out, SDtoBAF.out])
    File merged_BAF_index = select_first([MergeBAFEvidence.out_index, SDtoBAF.out_index])
    File merged_SR = MergeSREvidence.out
    File merged_SR_index = MergeSREvidence.out_index
    File merged_PE = MergePEEvidence.out
    File merged_PE_index = MergePEEvidence.out_index
  }
}

task MergeEvidence {
  input {
    Array[File] files
    String batch
    String evidence
    Array[String] samples
    File reference_dict
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  Int disk_size = 10 + ceil(size(files, "GiB") * 2)
  Int java_heap_size_mb = round(42.0 * length(files) + 1024.0)
  Float mem_size_gb = java_heap_size_mb / 1024.0 + 2.5

  RuntimeAttr default_attr = object {
    cpu_cores: 2,
    mem_gb: mem_size_gb,
    disk_gb: disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{batch}.~{evidence}.txt.gz"
    File out_index = "~{batch}.~{evidence}.txt.gz.tbi"
  }
  command <<<

    set -euo pipefail
    ulimit -n 100000

    mv ~{write_lines(files)} evidence.list
    mv ~{write_lines(samples)} samples.list

    awk '/txt\.gz$/' evidence.list | while read fil; do
      tabix -f -s1 -b2 -e2 $fil
    done

    /gatk/gatk --java-options "-Xmx~{java_heap_size_mb}m" PrintSVEvidence -F evidence.list --sample-names samples.list --sequence-dictionary ~{reference_dict} -O "~{batch}.~{evidence}.txt.gz"

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task SDtoBAF {
  input {
    Array[File] SD_files
    File sd_locs_vcf
    String batch
    Array[String] samples
    File reference_dict
    Float min_het_probability = 0.05
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  Int disk_size = 10 + ceil(size(SD_files, "GiB") * 2)
  Int java_heap_size_mb = round(42.0 * length(SD_files) + 1024.0)
  Float mem_size_gb = java_heap_size_mb / 1024.0 + 2.5

  RuntimeAttr default_attr = object {
    cpu_cores: 2,
    mem_gb: mem_size_gb,
    disk_gb: disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{batch}.baf.txt.gz"
    File out_index = "~{batch}.baf.txt.gz.tbi"
  }
  command <<<

    set -euo pipefail
    ulimit -n 100000

    mv ~{write_lines(SD_files)} inputs.list
    mv ~{write_lines(samples)} samples.list

    awk '/txt\.gz$/' inputs.list | while read fil; do
      tabix -f -s1 -b2 -e2 $fil
    done

    /gatk/gatk --java-options "-Xmx~{java_heap_size_mb}m" SiteDepthtoBAF \
        -F inputs.list \
        --sample-names samples.list \
        --sequence-dictionary "~{reference_dict}" \
        --baf-sites-vcf "~{sd_locs_vcf}" \
        --min-het-probability "~{min_het_probability}" \
        -O "~{batch}.baf.txt.gz"

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
