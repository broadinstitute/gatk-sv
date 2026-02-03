version 1.0

import "Structs.wdl"

workflow VcfVariantSizeSummary {
  input {
    File vcf_gz
    File vcf_tbi
    String sv_pipeline_base_docker

    RuntimeAttr? runtime_attr_detect_contigs
    RuntimeAttr? runtime_attr_count_contig
    RuntimeAttr? runtime_attr_merge
  }

  call DetectContigs {
    input:
      vcf_gz = vcf_gz,
      vcf_tbi = vcf_tbi,
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_detect_contigs
  }

  scatter (c in DetectContigs.contigs) {
    call CountVariantsPerContig {
      input:
        vcf_gz = vcf_gz,
        vcf_tbi = vcf_tbi,
        contig = c,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_count_contig
    }
  }

  call MergeCounts {
    input:
      count_tables = CountVariantsPerContig.count_table,
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_merge
  }

  output {
    File summary_table = MergeCounts.summary_table
  }
}

# ======================
# Detect contigs
# ======================
task DetectContigs {
  input {
    File vcf_gz
    File vcf_tbi
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(10 + size(vcf_gz, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    bcftools index -n ~{vcf_gz} > contigs.txt
  >>>

  output {
    Array[String] contigs = read_lines("contigs.txt")
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# ======================
# Count variants per contig
# ======================
task CountVariantsPerContig {
  input {
    File vcf_gz
    File vcf_tbi
    String contig
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(10 + size(vcf_gz, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    bcftools view -r ~{contig} ~{vcf_gz} -Ov | \
    awk -v c="~{contig}" '
    BEGIN {
      snv=0
      del_1_49=0
      ins_1_49=0
      del_gt49=0
      ins_gt49=0
    }
    !/^#/ {
      ref=$4
      alt=$5
      split(alt, alts, ",")
      for (i in alts) {
        rlen=length(ref)
        alen=length(alts[i])

        if (rlen==1 && alen==1) {
          snv++
        } else if (rlen > alen) {
          d = rlen - alen
          if (d <= 49) del_1_49++
          else del_gt49++
        } else if (alen > rlen) {
          d = alen - rlen
          if (d <= 49) ins_1_49++
          else ins_gt49++
        }
      }
    }
    END {
      print c "\t" snv "\t" del_1_49 "\t" ins_1_49 "\t" del_gt49 "\t" ins_gt49
    }' > ~{contig}.counts.tsv
  >>>

  output {
    File count_table = "~{contig}.counts.tsv"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# ======================
# Merge per-contig tables
# ======================
task MergeCounts {
  input {
    Array[File] count_tables
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: 20,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    echo -e "CONTIG\tSNV\tDEL_1_49\tINS_1_49\tDEL_GT49\tINS_GT49" > summary.tsv
    cat ~{sep=' ' count_tables} >> summary.tsv
  >>>

  output {
    File summary_table = "summary.tsv"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}