version 1.0

workflow NormalizeBamContigs {
  input {
    File input_bam
  }

  call NormalizeContigs {
    input:
      bam = input_bam
  }

  output {
    File bam = NormalizeContigs.out_bam
    File bai = NormalizeContigs.out_bai
  }
}

task NormalizeContigs {
  input {
    File bam
  }

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
    docker: "biocontainers/samtools:v1.17-1-deb_cv1"
    cpu: 2
    memory: "4G"
  }
}


