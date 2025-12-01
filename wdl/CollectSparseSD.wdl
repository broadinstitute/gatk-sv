version 1.0

import "Structs.wdl"

# Workflow to run PE/SR collection on a single sample
workflow CollectSparseSD {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_id
    File reference_fasta
    File reference_index
    File reference_dict
    File sparse_sd_locs_vcf
    File? primary_contigs_list
    File? gatk_jar_override
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  call RunCollectSparseSDEvidence {
    input:
      bam_or_cram_file = bam_or_cram_file,
      bam_or_cram_index = bam_or_cram_index,
      sample_id = sample_id,
      reference_fasta = reference_fasta,
      reference_index = reference_index,
      reference_dict = reference_dict,
      sparse_sd_locs_vcf = sparse_sd_locs_vcf,
      sparse_sd_locs_vcf_index = sparse_sd_locs_vcf + ".tbi",
      primary_contigs_list = primary_contigs_list,
      gatk_jar_override = gatk_jar_override,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File? sparse_sd_out = RunCollectSparseSDEvidence.sparse_sd_out
    File? sparse_sd_out_index = RunCollectSparseSDEvidence.sparse_sd_out_index
  }
}

task RunCollectSparseSDEvidence {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_id
    File reference_fasta
    File reference_index
    File reference_dict
    File sparse_sd_locs_vcf
    File sparse_sd_locs_vcf_index
    Int site_depth_min_mapq = 6
    Int site_depth_min_baseq = 10
    File? primary_contigs_list
    File? gatk_jar_override
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + 50)

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int command_mem_mb = ceil(mem_gb * 1000 - 500)

  output {
    File? sparse_sd_out = "${sample_id}.sparse.sd.txt.gz"
    File? sparse_sd_out_index = "${sample_id}.sparse.sd.txt.gz.tbi"
  }
  command <<<

    set -euo pipefail

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_jar_override}

    # Collect sparse site-depth for ploidy estimation (smaller site list)
    /gatk/gatk --java-options "-Xmx~{command_mem_mb}m" CollectSVEvidence \
        -I "~{bam_or_cram_file}" \
        --sample-name "~{sample_id}" \
        -F "~{sparse_sd_locs_vcf}" \
        -SD "~{sample_id}.sparse.sd.txt.gz" \
        --site-depth-min-mapq "~{site_depth_min_mapq}" \
        --site-depth-min-baseq "~{site_depth_min_baseq}" \
        ~{"-R " + reference_fasta} \
        ~{"-L " + primary_contigs_list} \
        --read-filter NonZeroReferenceLengthAlignmentReadFilter

    tabix -f -0 -s1 -b2 -e2 "~{sample_id}.sparse.sd.txt.gz"
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

