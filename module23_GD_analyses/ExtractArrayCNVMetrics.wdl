version 1.0

## Given one single-sample SNP-array VCF per sample and a region of
## interest, scatter across samples and run extract_array_cnv_metrics.sh
## per sample to pull out copy-number-relevant metrics (GT, BAF, LRR by
## default). Pivot the per-sample long-format tables into one wide matrix
## per metric: rows are sites, columns are samples, "." for a sample
## missing at a site.
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
    String python_docker = "python:3.11-slim"
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

  call PivotToWideMatrices {
    input:
      tables = ExtractPerSample.metrics_tsv,
      output_basename = output_basename,
      docker = python_docker,
  }

  output {
    Array[File] per_sample_tables = ExtractPerSample.metrics_tsv
    Array[File] wide_matrices = PivotToWideMatrices.wide_matrices
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

task PivotToWideMatrices {
  input {
    Array[File] tables
    String output_basename
    String docker
  }

  command <<<
    set -euo pipefail

    python3 <<CODE
import csv

table_files = "~{sep=',' tables}".split(",")
site_key_cols = ["CHROM", "POS", "ID", "REF", "ALT"]

sites = set()
samples = set()
metric_cols = None
data = {}

for path in table_files:
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if metric_cols is None:
            metric_cols = [c for c in reader.fieldnames if c not in site_key_cols and c != "SAMPLE"]
            for m in metric_cols:
                data[m] = {}
        for row in reader:
            key = tuple(row[c] for c in site_key_cols)
            sites.add(key)
            sample = row["SAMPLE"]
            samples.add(sample)
            for m in metric_cols:
                data[m].setdefault(key, {})[sample] = row[m]

sorted_sites = sorted(sites, key=lambda k: (k[0], int(k[1])))
sorted_samples = sorted(samples)

for m in metric_cols:
    out_path = f"~{output_basename}.{m}.tsv"
    with open(out_path, "w") as out:
        out.write("\t".join(site_key_cols + sorted_samples) + "\n")
        for key in sorted_sites:
            row_vals = list(key) + [data[m].get(key, {}).get(s, ".") for s in sorted_samples]
            out.write("\t".join(row_vals) + "\n")
    print(f"Wrote {out_path}: {len(sorted_sites)} sites x {len(sorted_samples)} samples")
CODE
  >>>

  output {
    Array[File] wide_matrices = glob("~{output_basename}.*.tsv")
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk 20 HDD"
    preemptible: 2
  }
}
