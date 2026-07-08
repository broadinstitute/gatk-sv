version 1.0

workflow CountVariantsByAFCutoff {
  input {
    Array[File] vcfs
    Float af_cutoff
    String output_prefix = "af_cutoff_variant_counts"
    String bcftools_docker = "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
  }

  call CountPerSampleVariantsByAF {
    input:
      vcfs = vcfs,
      af_cutoff = af_cutoff,
      output_prefix = output_prefix,
      bcftools_docker = bcftools_docker
  }

  output {
    File sample_counts_tsv = CountPerSampleVariantsByAF.sample_counts_tsv
  }
}

task CountPerSampleVariantsByAF {
  input {
    Array[File] vcfs
    Float af_cutoff
    String output_prefix
    String bcftools_docker
  }

  command <<<
    set -euo pipefail

    VCF_LIST="~{write_lines(vcfs)}"

    # Use the first VCF to define sample IDs and sample order.
    first_vcf=$(head -n 1 "${VCF_LIST}")
    bcftools query -l "${first_vcf}" > samples.list

    # Aggregate per-contig contributions as: sample_index, total_nonref_calls, low_af_nonref_calls
    : > per_vcf_counts.tsv

    while read -r vcf; do
      [ -z "${vcf}" ] && continue

      bcftools query -f '%INFO/AF[\t%GT]\n' "${vcf}" \
        | awk -v cutoff=~{af_cutoff} '
            BEGIN { OFS="\t" }
            {
              # AF can be comma-separated for multiallelic sites; use the minimum AF.
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
          ' >> per_vcf_counts.tsv
    done < "${VCF_LIST}"

    # Collapse counts across all VCFs.
    awk '
      BEGIN { OFS="\t" }
      { total[$1] += $2; low[$1] += $3 }
      END {
        for (idx in total) {
          print idx, total[idx], low[idx]
        }
      }
    ' per_vcf_counts.tsv | sort -k1,1n > counts_by_sample_index.tsv

    # Final output: sample_id, total_count, low_af_count
    awk '
      BEGIN { OFS="\t"; print "sample_id", "total_variant_count", "af_lt_cutoff_variant_count" }
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
    ' counts_by_sample_index.tsv samples.list > "~{output_prefix}.tsv"
  >>>

  output {
    File sample_counts_tsv = "~{output_prefix}.tsv"
  }

  runtime {
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk 30 HDD"
    docker: bcftools_docker
  }
}
