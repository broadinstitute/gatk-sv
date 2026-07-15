version 1.0

import "Structs.wdl"
import "SplitAndCountPerSampleNonRefVariantsPerContig.wdl" as PerContig

workflow SplitAndCountPerSampleNonRefVariants {
  input {
    Array[File] input_vcfs

    # If provided, only these samples are processed.
    # If omitted, all samples are derived from the first VCF header and
    # consistently processed across all contigs.
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

  if (!defined(sample_list)) {
    call PerContig.GetVcfSampleList as GetWorkflowSampleList {
      input:
        input_vcf             = input_vcfs[0],
        bcftools_docker       = bcftools_docker,
        runtime_attr_override = runtime_attr_sample_list
    }
  }

  Array[String] effective_samples = select_first([sample_list, GetWorkflowSampleList.sample_names])

  scatter (contig_vcf in input_vcfs) {
    call PerContig.SplitAndCountPerSampleNonRefVariantsPerContig as run_per_contig {
      input:
        input_vcf                 = contig_vcf,
        sample_list               = effective_samples,
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

  Array[File] all_per_contig_sample_nonref_vcfs = flatten(run_per_contig.sample_nonref_vcfs)

  scatter (sample_name in effective_samples) {
    call MergePerSampleNonRefVcfsAcrossContigs {
      input:
        sample_name             = sample_name,
        all_sample_nonref_vcfs  = all_per_contig_sample_nonref_vcfs,
        expected_contigs        = length(input_vcfs),
        output_prefix           = output_prefix,
        bcftools_docker         = bcftools_docker,
        runtime_attr_override   = runtime_attr_merge
    }
  }

  output {
    File? vcf_derived_sample_list = GetWorkflowSampleList.sample_list
    Array[String] processed_samples = effective_samples
    Array[Array[File]] per_contig_sample_nonref_vcfs = run_per_contig.sample_nonref_vcfs
    Array[File] per_sample_merged_nonref_vcfs = MergePerSampleNonRefVcfsAcrossContigs.merged_sample_nonref_vcf
    Array[File] per_sample_merged_nonref_vcf_indexes = MergePerSampleNonRefVcfsAcrossContigs.merged_sample_nonref_vcf_tbi
    Array[File] per_contig_merged_count_tables = run_per_contig.merged_sample_count_table
    Array[Array[File]] per_contig_per_sample_count_tables = run_per_contig.per_sample_count_tables

    File all_contigs_per_sample_table = SumAcrossContigs.per_contig_table
    File summed_across_contigs_per_sample_table = SumAcrossContigs.summed_table
  }
}


task MergePerSampleNonRefVcfsAcrossContigs {
  input {
    String sample_name
    Array[File] all_sample_nonref_vcfs
    Int expected_contigs
    String output_prefix
    String bcftools_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: 60,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  String merged_sample_nonref_vcf_name = output_prefix + "." + sample_name + ".nonref.vcf.gz"

  command <<<
    set -euo pipefail

    python3 << 'PY'
import json

sample = "~{sample_name}"
expected_contigs = int("~{expected_contigs}")
paths = "~{sep=' ' all_sample_nonref_vcfs}".split()
suffix = "." + sample + ".nonref.vcf.gz"
selected = [p for p in paths if p.endswith(suffix)]

if not selected:
    raise RuntimeError(f"No per-contig sample non-ref VCFs found for sample: {sample}")

if len(selected) != expected_contigs:
    raise RuntimeError(
        f"Expected {expected_contigs} per-contig VCFs for sample {sample}, found {len(selected)}"
    )

with open("sample_vcfs.list", "wt") as out:
    for p in selected:
        out.write(p + "\n")

print(json.dumps({
    "sample": sample,
    "expected_contigs": expected_contigs,
    "matched_vcfs": len(selected)
}))
PY

    bcftools concat \
      -f sample_vcfs.list \
      -Oz \
      -o ~{merged_sample_nonref_vcf_name}

    bcftools index -t ~{merged_sample_nonref_vcf_name}
  >>>

  output {
    File merged_sample_nonref_vcf = merged_sample_nonref_vcf_name
    File merged_sample_nonref_vcf_tbi = merged_sample_nonref_vcf_name + ".tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: bcftools_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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
