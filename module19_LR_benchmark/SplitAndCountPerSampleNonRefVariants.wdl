version 1.0

import "Structs.wdl"
import "SplitAndCountPerSampleNonRefVariantsPerContig.wdl" as PerContig

workflow SplitAndCountPerSampleNonRefVariants {
  input {
    Array[File] input_vcfs

    # If provided, only these samples are processed per contig.
    # If omitted, the sample list is derived from the first VCF header.
    Array[String]? sample_list

    String output_prefix = "per_sample_nonref"
    String bcftools_docker = "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    String python_docker = "python:3.11-slim"

    RuntimeAttr? runtime_attr_sample_list
    RuntimeAttr? runtime_attr_extract
    RuntimeAttr? runtime_attr_split_count
    RuntimeAttr? runtime_attr_merge
    RuntimeAttr? runtime_attr_sum
  }

  scatter (contig_vcf in input_vcfs) {
    call PerContig.SplitAndCountPerSampleNonRefVariantsPerContig as run_per_contig {
      input:
        input_vcf                 = contig_vcf,
        sample_list               = sample_list,
        output_prefix             = output_prefix,
        bcftools_docker           = bcftools_docker,
        python_docker             = python_docker,
        runtime_attr_sample_list  = runtime_attr_sample_list,
        runtime_attr_extract      = runtime_attr_extract,
        runtime_attr_split_count  = runtime_attr_split_count,
        runtime_attr_merge        = runtime_attr_merge
    }
  }

  call SumAcrossContigs {
    input:
      per_contig_count_tables = run_per_contig.merged_sample_count_table,
      per_contig_output       = output_prefix + ".sample_category_counts.per_contig.tsv",
      summed_output           = output_prefix + ".sample_category_counts.summed_across_contigs.tsv",
      python_docker           = python_docker,
      runtime_attr_override   = runtime_attr_sum
  }

  output {
    Array[File] per_contig_merged_count_tables = run_per_contig.merged_sample_count_table
    Array[Array[File]] per_contig_per_sample_count_tables = run_per_contig.per_sample_count_tables

    File all_contigs_per_sample_table = SumAcrossContigs.per_contig_table
    File summed_across_contigs_per_sample_table = SumAcrossContigs.summed_table
  }
}


task SumAcrossContigs {
  input {
    Array[File] per_contig_count_tables
    String per_contig_output
    String summed_output
    String python_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: 20,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python3 << 'PY'
import csv

inputs = "~{sep=' ' per_contig_count_tables}".split()
per_contig_path = "~{per_contig_output}"
summed_path = "~{summed_output}"

if not inputs:
    raise RuntimeError("No per-contig count tables provided")

header = None
numeric_cols = []
sum_by_sample = {}
rows_written = 0

with open(per_contig_path, "wt", newline="") as out_all:
    writer_all = None

    for path in inputs:
        with open(path, "rt", newline="") as fh:
            reader = csv.reader(fh, delimiter='\t')
            this_header = next(reader)

            if header is None:
                header = this_header
                if len(header) < 3 or header[0] != "contig" or header[1] != "sample":
                    raise RuntimeError("Expected per-contig table header: contig, sample, ...")
                numeric_cols = header[2:]
                writer_all = csv.writer(out_all, delimiter='\t')
                writer_all.writerow(header)
            elif this_header != header:
                raise RuntimeError(f"Header mismatch in {path}")

            for row in reader:
                writer_all.writerow(row)
                rows_written += 1

                sample = row[1]
                if sample not in sum_by_sample:
                    sum_by_sample[sample] = [0] * len(numeric_cols)

                for i in range(len(numeric_cols)):
                    raw = row[i + 2].strip()
                    val = int(float(raw)) if raw != "" else 0
                    sum_by_sample[sample][i] += val

with open(summed_path, "wt", newline="") as out_sum:
    writer = csv.writer(out_sum, delimiter='\t')
    writer.writerow(["sample"] + numeric_cols)
    for sample in sorted(sum_by_sample.keys()):
        writer.writerow([sample] + sum_by_sample[sample])

print({
    "inputs": len(inputs),
    "rows_per_contig_table": rows_written,
    "samples": len(sum_by_sample),
    "per_contig_output": per_contig_path,
    "summed_output": summed_path,
})
PY
  >>>

  output {
    File per_contig_table = per_contig_output
    File summed_table = summed_output
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: python_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
