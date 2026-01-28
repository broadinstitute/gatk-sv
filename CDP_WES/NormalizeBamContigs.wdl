version 1.0

import "Structs.wdl"

workflow NormalizeBamContigsScatter {
  input {
    File input_bam
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_detect_contigs
    RuntimeAttr? runtime_attr_convert_contig
    RuntimeAttr? runtime_attr_merge_bams
  }

  call DetectContigs {
    input:
      bam = input_bam,
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_detect_contigs
  }

  if (DetectContigs.needs_conversion) {
    scatter (c in DetectContigs.canonical_contigs) {
      call ConvertOneContig {
        input:
          bam = input_bam,
          contig = c,
          docker_image = sv_pipeline_base_docker,
          runtime_attr_override = runtime_attr_convert_contig

      }
    }

    call MergeBams {
      input:
        bams = ConvertOneContig.out_bam,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_merge_bams
    }
  }

  output {
    File bam = select_first([
      MergeBams.out_bam,
      input_bam
    ])
    File bai = select_first([
      MergeBams.out_bai,
      input_bam + ".bai"
    ])
  }
}


task DetectContigs {
  input {
    File bam
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: ceil(10 + size(bam, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    samtools view -H ~{bam} > header.txt

    if grep -q "^@SQ.*SN:[0-9]" header.txt; then
      echo true > needs_conversion.txt
    else
      echo false > needs_conversion.txt
    fi

    # Canonical contigs only
    printf "1\n2\n3\n4\n5\n6\n7\n8\n9\n10\n11\n12\n13\n14\n15\n16\n17\n18\n19\n20\n21\n22\nX\nY\nMT\n" > contigs.txt
  >>>

  output {
    Boolean needs_conversion = read_boolean("needs_conversion.txt")
    Array[String] canonical_contigs = read_lines("contigs.txt")
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

task ConvertOneContig {
  input {
    File bam
    String contig
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: ceil(10 + size(bam, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    if [ "~{contig}" = "MT" ]; then
      new_contig="chrM"
    else
      new_contig="chr~{contig}"
    fi

    samtools view -h ~{bam} "~{contig}" | \
    awk -v nc="$new_contig" '
      BEGIN { OFS="\t" }
      /^@/ {
        if ($1 == "@SQ" && $0 ~ "SN:~{contig}") {
          sub("SN:~{contig}", "SN:" nc)
        }
        print
        next
      }
      {
        $3 = nc
        if ($7 != "=" && $7 != "*") {
          $7 = nc
        }
        print
      }
    ' | samtools view -b -o ~{contig}.bam

    samtools index ~{contig}.bam
  >>>

  output {
    File out_bam = "~{contig}.bam"
    File out_bai = "~{contig}.bam.bai"
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

task NormalizeContigs {
  input {
    File bam
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: ceil(10 + size(bam, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    samtools view -H ~{bam} > header.txt

    # Check whether contigs are numeric (1,2,3,...) rather than chr*
    if grep -q "^@SQ.*SN:[0-9]" header.txt; then
      echo "Detected numeric contigs, converting to chr-prefixed contigs"

      samtools view -h ~{bam} | \
      awk '
      BEGIN { OFS="\t" }

      # Header
      /^@/ {
        if ($1 == "@SQ") {
          for (i=1; i<=NF; i++) {
            if ($i ~ /^SN:/) {
              split($i, a, ":")
              c = a[2]
              if (c == "MT") {
                $i = "SN:chrM"
              } else if (c ~ /^(1[0-9]|2[0-2]|[1-9]|X|Y)$/) {
                $i = "SN:chr" c
              } else {
                next
              }
            }
          }
        }
        print
        next
      }

      # Alignment records
      {
        # RNAME
        if ($3 == "MT") {
          $3 = "chrM"
        } else if ($3 ~ /^(1[0-9]|2[0-2]|[1-9]|X|Y)$/) {
          $3 = "chr" $3
        } else {
          next
        }

        # RNEXT
        if ($7 == "=") {
          # do nothing
        } else if ($7 == "MT") {
          $7 = "chrM"
        } else if ($7 ~ /^(1[0-9]|2[0-2]|[1-9]|X|Y)$/) {
          $7 = "chr" $7
        } else {
          $7 = "*"
        }

        print
      }
      ' | samtools view -b -o out.bam

    else
      echo "Contigs already chr-prefixed; passing BAM through"
      cp ~{bam} out.bam
    fi

    samtools index out.bam
  >>>

  output {
    File out_bam = "out.bam"
    File out_bai = "out.bam.bai"
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

task MergeBams {
  input {
    Array[File] bams
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: ceil(10 + size(s, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    samtools merge -f merged.bam ~{sep=' ' bams}
    samtools index merged.bam
  >>>

  output {
    File out_bam = "merged.bam"
    File out_bai = "merged.bam.bai"
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
