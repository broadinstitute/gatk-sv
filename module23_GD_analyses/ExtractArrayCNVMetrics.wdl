version 1.0

## Given one single-sample SNP-array VCF per sample and a region of
## interest, scatter across samples and run extract_array_cnv_metrics.sh
## per sample to pull out copy-number-relevant metrics (GT, BAF, LRR by
## default). Aggregate the per-sample tables into one combined long-format
## table (one row per variant x sample).
##
## extract_array_cnv_metrics.sh handles compression and indexing itself, so
## each input VCF can be plain, gzipped, or bgzipped, indexed or not.

workflow ExtractArrayCNVMetrics {
  input {
    Array[File] vcfs
    String region
    File script
    Boolean extra_fields = false
    String output_basename = "array_cnv_metrics"
    String docker = "staphb/bcftools:1.19"
  }

  scatter (vcf in vcfs) {
    call ExtractPerSample {
      input:
        vcf = vcf,
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
    String region
    File script
    Boolean extra_fields
    String docker
  }

  Int disk_gb = ceil(size(vcf, "GB") * 2) + 20

  command <<<
    set -euo pipefail

    sample=$(bcftools query -l ~{vcf} | head -n1)
    bash ~{script} ~{vcf} ~{region} "${sample}.array_cnv_metrics.tsv" ~{true="--extra" false="" extra_fields}
  >>>

  output {
    File metrics_tsv = glob("*.array_cnv_metrics.tsv")[0]
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
