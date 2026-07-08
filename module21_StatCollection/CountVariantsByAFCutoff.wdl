version 1.0

workflow CountVariantsByAFCutoff {
  input {
    Array[File] vcfs
    Float af_cutoff
    String output_prefix = "af_cutoff_variant_counts"
    String bcftools_docker = "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
  }

  scatter (vcf in vcfs) {
    call CountPerContigVariantsByAF {
      input:
        vcf = vcf,
        af_cutoff = af_cutoff,
        bcftools_docker = bcftools_docker
    }
  }

  call SumCountsAcrossContigs {
    input:
      per_contig_count_tables = CountPerContigVariantsByAF.per_contig_counts_tsv,
      output_prefix = output_prefix,
      bcftools_docker = bcftools_docker
  }

  output {
    Array[File] per_contig_sample_counts = CountPerContigVariantsByAF.per_contig_counts_tsv
    File sample_counts_tsv = SumCountsAcrossContigs.sample_counts_tsv
  }
}

task CountPerContigVariantsByAF {
  input {
    File vcf
    Float af_cutoff
    String bcftools_docker
  }

  command <<<
    set -euo pipefail

    bcftools query -l "~{vcf}" > samples.list

    bcftools query -f '%INFO/AF[\t%GT]\n' "~{vcf}" \
      | awk -v cutoff=~{af_cutoff} '
          BEGIN { OFS="\t" }
          {
            n_af = split($1, afs, ",")
            min_af = ""
            for (j = 1; j <= n_af; j++) {
              if (afs[j] == "." || afs[j] == "NA" || afs[j] == "") continue
              af_val = afs[j] + 0
              if (min_af == "" || af_val < min_af) min_af = af_val
            }
            if (min_af == "") next

            for (i = 2; i <= NF; i++) {
              gt = $i
              if (gt == "." || gt == "./." || gt == ".|.") continue
              if (gt ~ /[1-9]/) {
                total[i-1]++
                if (min_af < cutoff) low[i-1]++
              }
            }
          }
          END {
            for (idx in total) {
              print idx, total[idx], low[idx] + 0
            }
          }
        ' | sort -k1,1n > counts_by_sample_index.tsv

    awk '
      BEGIN { OFS="\t" }
      NR==FNR {
        total[$1] = $2
        low[$1] = $3
        next
      }
      {
        idx = FNR
        t = (idx in total ? total[idx] : 0)
        l = (idx in low ? low[idx] : 0)
        print $1, t, l
      }
    ' counts_by_sample_index.tsv samples.list > per_contig_counts.tsv
  >>>

  output {
    File per_contig_counts_tsv = "per_contig_counts.tsv"
  }

  runtime {
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk 20 HDD"
    docker: bcftools_docker
  }
}

task SumCountsAcrossContigs {
  input {
    Array[File] per_contig_count_tables
    String output_prefix
    String bcftools_docker
  }

  command <<<
    set -euo pipefail

    cat ~{sep=" " per_contig_count_tables} \
      | awk '
          BEGIN { OFS="\t" }
          { total[$1] += $2; low[$1] += $3 }
          END {
            for (sample in total) {
              print sample, total[sample], low[sample]
            }
          }
        ' \
      | sort -k1,1 > summed_counts.tsv

    {
      echo -e "sample_id\ttotal_variant_count\taf_lt_cutoff_variant_count"
      cat summed_counts.tsv
    } > "~{output_prefix}.tsv"
  >>>

  output {
    File sample_counts_tsv = "~{output_prefix}.tsv"
  }

  runtime {
    cpu: 1
    memory: "2 GiB"
    disks: "local-disk 10 HDD"
    docker: bcftools_docker
  }
}
