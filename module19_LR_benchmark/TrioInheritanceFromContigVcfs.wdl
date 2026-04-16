version 1.0

import "Structs.wdl"

workflow TrioInheritanceFromContigVcfs {
  input {
    Array[File] contig_vcfs

    File split_vcf_script
    File split_trv_by_motifs_script
    File count_trio_script
    File analyze_inheritance_script
    File trio_table_tsv

    Boolean pass_only = false
    String output_prefix = "trio_inheritance"
    String python_docker = "python:3.11-slim"

    RuntimeAttr? runtime_attr_process_contig
    RuntimeAttr? runtime_attr_aggregate_subset
  }

  scatter (i in range(length(contig_vcfs))) {
    String contig_label = basename(contig_vcfs[i], ".vcf.gz")

    call ProcessContigVcf {
      input:
        contig_vcf                = contig_vcfs[i],
        contig_label              = contig_label,
        split_vcf_script          = split_vcf_script,
        split_trv_by_motifs_script = split_trv_by_motifs_script,
        count_trio_script         = count_trio_script,
        analyze_inheritance_script = analyze_inheritance_script,
        trio_table_tsv            = trio_table_tsv,
        pass_only                 = pass_only,
        python_docker             = python_docker,
        runtime_attr_override     = runtime_attr_process_contig
    }
  }

  call AggregateSubsetSummary as aggregate_trv {
    input:
      summary_files          = ProcessContigVcf.summary_trv,
      subset_name            = "trv",
      output_prefix          = output_prefix,
      python_docker          = python_docker,
      runtime_attr_override  = runtime_attr_aggregate_subset
  }

  call AggregateSubsetSummary as aggregate_non_trv_snv {
    input:
      summary_files          = ProcessContigVcf.summary_non_trv_snv,
      subset_name            = "non_trv_snv",
      output_prefix          = output_prefix,
      python_docker          = python_docker,
      runtime_attr_override  = runtime_attr_aggregate_subset
  }

  call AggregateSubsetSummary as aggregate_non_trv_ins_lt50 {
    input:
      summary_files          = ProcessContigVcf.summary_non_trv_ins_lt50,
      subset_name            = "non_trv_ins_lt50",
      output_prefix          = output_prefix,
      python_docker          = python_docker,
      runtime_attr_override  = runtime_attr_aggregate_subset
  }

  call AggregateSubsetSummary as aggregate_non_trv_del_lt50 {
    input:
      summary_files          = ProcessContigVcf.summary_non_trv_del_lt50,
      subset_name            = "non_trv_del_lt50",
      output_prefix          = output_prefix,
      python_docker          = python_docker,
      runtime_attr_override  = runtime_attr_aggregate_subset
  }

  call AggregateSubsetSummary as aggregate_non_trv_ins_ge50 {
    input:
      summary_files          = ProcessContigVcf.summary_non_trv_ins_ge50,
      subset_name            = "non_trv_ins_ge50",
      output_prefix          = output_prefix,
      python_docker          = python_docker,
      runtime_attr_override  = runtime_attr_aggregate_subset
  }

  call AggregateSubsetSummary as aggregate_non_trv_del_ge50 {
    input:
      summary_files          = ProcessContigVcf.summary_non_trv_del_ge50,
      subset_name            = "non_trv_del_ge50",
      output_prefix          = output_prefix,
      python_docker          = python_docker,
      runtime_attr_override  = runtime_attr_aggregate_subset
  }

  call AggregateSubsetSummary as aggregate_trv_motif_1bp {
    input:
      summary_files          = ProcessContigVcf.summary_trv_motif_1bp,
      subset_name            = "trv_motif_1bp",
      output_prefix          = output_prefix,
      python_docker          = python_docker,
      runtime_attr_override  = runtime_attr_aggregate_subset
  }

  call AggregateSubsetSummary as aggregate_trv_motif_2bp {
    input:
      summary_files          = ProcessContigVcf.summary_trv_motif_2bp,
      subset_name            = "trv_motif_2bp",
      output_prefix          = output_prefix,
      python_docker          = python_docker,
      runtime_attr_override  = runtime_attr_aggregate_subset
  }

  call AggregateSubsetSummary as aggregate_trv_motif_rest {
    input:
      summary_files          = ProcessContigVcf.summary_trv_motif_rest,
      subset_name            = "trv_motif_rest",
      output_prefix          = output_prefix,
      python_docker          = python_docker,
      runtime_attr_override  = runtime_attr_aggregate_subset
  }

  output {
    Array[File] contig_summary_trv = ProcessContigVcf.summary_trv
    Array[File] contig_summary_non_trv_snv = ProcessContigVcf.summary_non_trv_snv
    Array[File] contig_summary_non_trv_ins_lt50 = ProcessContigVcf.summary_non_trv_ins_lt50
    Array[File] contig_summary_non_trv_del_lt50 = ProcessContigVcf.summary_non_trv_del_lt50
    Array[File] contig_summary_non_trv_ins_ge50 = ProcessContigVcf.summary_non_trv_ins_ge50
    Array[File] contig_summary_non_trv_del_ge50 = ProcessContigVcf.summary_non_trv_del_ge50
    Array[File] contig_summary_trv_motif_1bp = ProcessContigVcf.summary_trv_motif_1bp
    Array[File] contig_summary_trv_motif_2bp = ProcessContigVcf.summary_trv_motif_2bp
    Array[File] contig_summary_trv_motif_rest = ProcessContigVcf.summary_trv_motif_rest

    File merged_summary_trv = aggregate_trv.merged_summary
    File merged_summary_non_trv_snv = aggregate_non_trv_snv.merged_summary
    File merged_summary_non_trv_ins_lt50 = aggregate_non_trv_ins_lt50.merged_summary
    File merged_summary_non_trv_del_lt50 = aggregate_non_trv_del_lt50.merged_summary
    File merged_summary_non_trv_ins_ge50 = aggregate_non_trv_ins_ge50.merged_summary
    File merged_summary_non_trv_del_ge50 = aggregate_non_trv_del_ge50.merged_summary
    File merged_summary_trv_motif_1bp = aggregate_trv_motif_1bp.merged_summary
    File merged_summary_trv_motif_2bp = aggregate_trv_motif_2bp.merged_summary
    File merged_summary_trv_motif_rest = aggregate_trv_motif_rest.merged_summary
  }
}


task ProcessContigVcf {
  input {
    File contig_vcf
    String contig_label

    File split_vcf_script
    File split_trv_by_motifs_script
    File count_trio_script
    File analyze_inheritance_script
    File trio_table_tsv

    Boolean pass_only
    String python_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 2,
    mem_gb: 8,
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String subset_prefix = contig_label + ".split"

  command <<<
    set -euo pipefail

    python3 ~{split_vcf_script} \
      --input ~{contig_vcf} \
      --out-prefix ~{subset_prefix}

    python3 ~{split_trv_by_motifs_script} \
      --input ~{contig_vcf} \
      --out-prefix ~{subset_prefix}

    for subset in trv non_trv.snv non_trv.ins_lt50 non_trv.del_lt50 non_trv.ins_ge50 non_trv.del_ge50; do
      subset_vcf="~{subset_prefix}.${subset}.vcf.gz"
      pattern_out="~{subset_prefix}.${subset}.inheritance.tsv"
      summary_out="~{subset_prefix}.${subset}.summary.tsv"
      annotated_out="~{subset_prefix}.${subset}.annotated.tsv"

      python3 ~{count_trio_script} \
        --vcf "$subset_vcf" \
        --trios ~{trio_table_tsv} \
        --output "$pattern_out" \
        ~{if pass_only then "--pass-only" else ""}

      python3 ~{analyze_inheritance_script} \
        "$pattern_out" \
        "$summary_out" \
        "$annotated_out"
    done

    for subset in trv.motif_1bp trv.motif_2bp trv.motif_rest; do
      subset_vcf="~{subset_prefix}.${subset}.vcf.gz"
      pattern_out="~{subset_prefix}.${subset}.inheritance.tsv"
      summary_out="~{subset_prefix}.${subset}.summary.tsv"
      annotated_out="~{subset_prefix}.${subset}.annotated.tsv"

      python3 ~{count_trio_script} \
        --vcf "$subset_vcf" \
        --trios ~{trio_table_tsv} \
        --output "$pattern_out" \
        ~{if pass_only then "--pass-only" else ""}

      python3 ~{analyze_inheritance_script} \
        "$pattern_out" \
        "$summary_out" \
        "$annotated_out"
    done
  >>>

  output {
    File summary_trv = subset_prefix + ".trv.summary.tsv"
    File summary_non_trv_snv = subset_prefix + ".non_trv.snv.summary.tsv"
    File summary_non_trv_ins_lt50 = subset_prefix + ".non_trv.ins_lt50.summary.tsv"
    File summary_non_trv_del_lt50 = subset_prefix + ".non_trv.del_lt50.summary.tsv"
    File summary_non_trv_ins_ge50 = subset_prefix + ".non_trv.ins_ge50.summary.tsv"
    File summary_non_trv_del_ge50 = subset_prefix + ".non_trv.del_ge50.summary.tsv"
    File summary_trv_motif_1bp = subset_prefix + ".trv.motif_1bp.summary.tsv"
    File summary_trv_motif_2bp = subset_prefix + ".trv.motif_2bp.summary.tsv"
    File summary_trv_motif_rest = subset_prefix + ".trv.motif_rest.summary.tsv"
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


task AggregateSubsetSummary {
  input {
    Array[File] summary_files
    String subset_name
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
  String out_file = output_prefix + "." + subset_name + ".across_contigs.summary.tsv"

  command <<<
    set -euo pipefail

    python3 << 'PY'
import csv
from collections import defaultdict

files = "~{sep=' ' summary_files}".split()
outfile = "~{out_file}"

if not files:
    raise RuntimeError("No summary files were provided to aggregate.")

header = None
metric_cols = None
acc = defaultdict(lambda: defaultdict(int))

for path in files:
    with open(path, 'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if reader.fieldnames is None:
            continue

        if header is None:
            header = reader.fieldnames
            metric_cols = [c for c in header if c not in ('family_ID', 'summary_type')]

        for row in reader:
            key = (row['family_ID'], row['summary_type'])
            for col in metric_cols:
                val = row.get(col, '')
                if val == '' or val is None:
                    continue
                acc[key][col] += int(val)

with open(outfile, 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(header)
    for family_id, summary_type in sorted(acc.keys()):
        writer.writerow([
            family_id,
            summary_type,
            *[str(acc[(family_id, summary_type)].get(col, 0)) for col in metric_cols],
        ])
PY
  >>>

  output {
    File merged_summary = out_file
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
