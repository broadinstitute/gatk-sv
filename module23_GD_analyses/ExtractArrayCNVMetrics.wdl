version 1.0

## Given a multi-sample SNP-array VCF (+ index), a list of sample IDs, and a
## region of interest, scatter across samples -- subsetting each to its own
## single-sample VCF for that region -- and run extract_array_cnv_metrics.sh
## per sample to pull out copy-number-relevant metrics (GT, BAF, LRR by
## default). Aggregate the per-sample tables into one combined long-format
## table (one row per variant x sample).

workflow ExtractArrayCNVMetrics {
  input {
    File vcf
    File vcf_idx
    File sample_list
    String region
    File script
    Boolean extra_fields = false
    String output_basename = "array_cnv_metrics"
    String docker = "staphb/bcftools:1.19"
  }

  Array[String] samples = read_lines(sample_list)

  scatter (sample in samples) {
    call ExtractPerSample {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        sample = sample,
        region = region,
        script = script,
        extra_fields = extra_fields,
        docker = docker,
    }
  }

  call AggregateTables {
    input:
      tables = ExtractPerSample.metrics_tsv,
      output_basename = output_basename,
      docker = docker,
  }

  output {
    Array[File] per_sample_tables = ExtractPerSample.metrics_tsv
    File aggregated_table = AggregateTables.combined_tsv
  }
}

task ExtractPerSample {
  input {
    File vcf
    File vcf_idx
    String sample
    String region
    File script
    Boolean extra_fields
    String docker
  }

  Int disk_gb = ceil(size(vcf, "GB") * 2) + 20

  command <<<
    set -euo pipefail

    bcftools view -s ~{sample} -r ~{region} ~{vcf} -Oz -o ~{sample}.subset.vcf.gz
    bcftools index -t ~{sample}.subset.vcf.gz

    bash ~{script} ~{sample}.subset.vcf.gz ~{region} ~{sample}.array_cnv_metrics.tsv ~{true="--extra" false="" extra_fields}
  >>>

  output {
    File metrics_tsv = "~{sample}.array_cnv_metrics.tsv"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}

task AggregateTables {
  input {
    Array[File] tables
    String output_basename
    String docker
  }

  command <<<
    set -euo pipefail

    files=(~{sep=" " tables})
    head -n 1 "${files[0]}" > ~{output_basename}.tsv
    for f in "${files[@]}"; do
      tail -n +2 "$f" >> ~{output_basename}.tsv
    done
  >>>

  output {
    File combined_tsv = "~{output_basename}.tsv"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GiB"
    disks: "local-disk 20 HDD"
    preemptible: 2
  }
}
