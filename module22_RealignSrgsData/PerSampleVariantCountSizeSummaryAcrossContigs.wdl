version 1.0

import "Structs.wdl"
import "PerSampleVariantCountSizeSummary.wdl" as PerContig

## Scatters PerSampleVariantCountSizeSummary across one VCF per chromosome,
## then sums the per-sample variant COUNT and SIZE tables across all contigs
## into genome-wide totals (rows = sample, cols = SVTYPE / component type).

workflow PerSampleVariantCountSizeSummaryAcrossContigs {
  input {
    Array[File] contig_vcf_gzs
    Array[File] contig_vcf_tbis
    String output_prefix
    String docker = "python:3.11-slim"

    RuntimeAttr? runtime_attr_sample_ids
    RuntimeAttr? runtime_attr_extract
    RuntimeAttr? runtime_attr_summarize
    RuntimeAttr? runtime_attr_sum
  }

  scatter (vcf_pair in zip(contig_vcf_gzs, contig_vcf_tbis)) {
    call PerContig.PerSampleVariantCountSizeSummary as run_per_contig {
      input:
        vcf_gz = vcf_pair.left,
        vcf_tbi = vcf_pair.right,
        output_prefix = output_prefix,
        docker = docker,
        runtime_attr_sample_ids = runtime_attr_sample_ids,
        runtime_attr_extract = runtime_attr_extract,
        runtime_attr_summarize = runtime_attr_summarize
    }
  }

  call SumWideTablesAcrossContigs as SumCounts {
    input:
      per_contig_tables = run_per_contig.variant_counts_table,
      output_name = output_prefix + ".variant_counts_per_sample.genome_wide.tsv",
      docker = docker,
      runtime_attr_override = runtime_attr_sum
  }

  call SumWideTablesAcrossContigs as SumSizes {
    input:
      per_contig_tables = run_per_contig.variant_sizes_table,
      output_name = output_prefix + ".variant_sizes_per_sample.genome_wide.tsv",
      docker = docker,
      runtime_attr_override = runtime_attr_sum
  }

  output {
    Array[File] per_contig_nonref_variants_tables = run_per_contig.nonref_variants_table
    Array[File] per_contig_variant_counts_tables = run_per_contig.variant_counts_table
    Array[File] per_contig_variant_sizes_tables = run_per_contig.variant_sizes_table
    File genome_wide_variant_counts_table = SumCounts.summed_table
    File genome_wide_variant_sizes_table = SumSizes.summed_table
  }
}

task SumWideTablesAcrossContigs {
  input {
    Array[File] per_contig_tables
    String output_name
    String docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(10 + size(per_contig_tables, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 2,
    max_retries: 0
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python3 << 'EOF'
    out_path = "~{output_name}"
    paths = "~{sep=' ' per_contig_tables}".split()

    if not paths:
        raise RuntimeError("No per-contig tables provided")

    all_cols = []
    seen_cols = set()
    sums = {}

    for path in paths:
        with open(path) as f:
            header = f.readline().rstrip("\n").split("\t")
            cols = header[1:]
            for c in cols:
                if c not in seen_cols:
                    seen_cols.add(c)
                    all_cols.append(c)

            for line in f:
                fields = line.rstrip("\n").split("\t")
                sample_id = fields[0]
                values = fields[1:]
                row = sums.setdefault(sample_id, {})
                for col, val in zip(cols, values):
                    row[col] = row.get(col, 0) + int(val)

    with open(out_path, "w") as out:
        out.write("sample_id\t" + "\t".join(all_cols) + "\n")
        for sample_id in sorted(sums.keys()):
            row = sums[sample_id]
            out.write(sample_id + "\t" + "\t".join(str(row.get(c, 0)) for c in all_cols) + "\n")

    print(f"n_contigs={len(paths)} n_samples={len(sums)} n_cols={len(all_cols)}")
    EOF
  >>>

  output {
    File summed_table = output_name
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
