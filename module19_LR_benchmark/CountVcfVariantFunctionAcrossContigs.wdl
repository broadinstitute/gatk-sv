version 1.0

import "Structs.wdl"

workflow CountVcfVariantFunctionAcrossContigs {
  input {
    Array[File] contig_vcfs
    File count_script

    String output_prefix = "vcf_variant_function"
    String vep_key = "vep"
    String python_docker = "python:3.11-slim"

    RuntimeAttr? runtime_attr_process_contig
    RuntimeAttr? runtime_attr_aggregate
  }

  scatter (i in range(length(contig_vcfs))) {
    String contig_label = basename(contig_vcfs[i], ".vcf.gz")

    call ProcessContigVcfFunctionSummary {
      input:
        contig_vcf            = contig_vcfs[i],
        contig_label          = contig_label,
        count_script          = count_script,
        output_prefix         = output_prefix,
        vep_key               = vep_key,
        python_docker         = python_docker,
        runtime_attr_override = runtime_attr_process_contig
    }
  }

  call AggregateContigFunctionTables {
    input:
      summary_tables         = ProcessContigVcfFunctionSummary.summary_table,
      output_prefix          = output_prefix,
      python_docker          = python_docker,
      runtime_attr_override  = runtime_attr_aggregate
  }

  output {
    Array[File] per_contig_tables = ProcessContigVcfFunctionSummary.summary_table
    File merged_table = AggregateContigFunctionTables.merged_table
    File concatenated_table = AggregateContigFunctionTables.concatenated_table
  }
}


task ProcessContigVcfFunctionSummary {
  input {
    File contig_vcf
    String contig_label
    File count_script
    String output_prefix
    String vep_key
    String python_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  String out_file = output_prefix + "." + contig_label + ".variant_function_summary.tsv"

  command <<<
    set -euo pipefail

    python3 ~{count_script} \
      --input ~{contig_vcf} \
      --output ~{out_file} \
      --vep-key ~{vep_key}
  >>>

  output {
    File summary_table = out_file
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: python_docker
  }
}


task AggregateContigFunctionTables {
  input {
    Array[File] summary_tables
    String output_prefix
    String python_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: 20,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  String merged_out = output_prefix + ".across_contigs.variant_function_summary.tsv"
  String concatenated_out = output_prefix + ".per_contig.variant_function_summary.tsv"

  command <<<
    set -euo pipefail

    python3 - <<'PY'
import csv
import os
from collections import OrderedDict

summary_files = """~{sep='\n' summary_tables}""".strip().splitlines()
merged_out = "~{merged_out}"
concatenated_out = "~{concatenated_out}"

if not summary_files or summary_files == ['']:
    raise RuntimeError("No summary tables provided.")

header = None
metric_cols = None
agg = OrderedDict()

with open(concatenated_out, "w", newline="") as concat_fh:
    concat_writer = None

    for fpath in summary_files:
        contig_label = os.path.basename(fpath)
        suffix = ".variant_function_summary.tsv"
        if contig_label.endswith(suffix):
            contig_label = contig_label[:-len(suffix)]

        with open(fpath, "r") as fh:
            data_lines = [line.rstrip("\n") for line in fh if line.strip() and not line.startswith("#")]

        if not data_lines:
            continue

        local_header = data_lines[0].split("\t")
        if header is None:
            header = local_header
            metric_cols = header[1:]
            concat_writer = csv.writer(concat_fh, delimiter="\t")
            concat_writer.writerow(["contig"] + header)
        elif local_header != header:
            raise RuntimeError(f"Header mismatch in {fpath}: {local_header} != {header}")

        for line in data_lines[1:]:
            fields = line.split("\t")
            vclass = fields[0]
            vals = [int(x) for x in fields[1:]]

            if vclass not in agg:
                agg[vclass] = [0] * len(metric_cols)
            for i, v in enumerate(vals):
                agg[vclass][i] += v

            concat_writer.writerow([contig_label] + fields)

if header is None:
    raise RuntimeError("No data rows found in summary tables.")

with open(merged_out, "w", newline="") as out_fh:
    writer = csv.writer(out_fh, delimiter="\t")
    writer.writerow(header)
    for vclass, vals in agg.items():
        writer.writerow([vclass] + vals)
PY
  >>>

  output {
    File merged_table = merged_out
    File concatenated_table = concatenated_out
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: python_docker
  }
}
