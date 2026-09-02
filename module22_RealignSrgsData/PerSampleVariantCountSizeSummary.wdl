version 1.0

import "Structs.wdl"

## Given a single (per-contig) SV VCF:
##   1. ExtractSampleIds        -- list all sample IDs in the VCF header.
##   2. ExtractNonrefPassVariants -- one pass over the VCF, emitting a long
##      table (sample_id, ID, SVTYPE, SVLEN, CPX_INTERVALS) for every
##      (sample, variant) pair where the sample's genotype is non-ref and
##      the record's FILTER is PASS.
##   3. SummarizeCountsAndSizes -- from that table, compute per-sample:
##        - variant COUNT per SVTYPE (DEL, INS, DUP, INV, BND, CPX, ... --
##          a CPX record counts once, toward "CPX", regardless of how many
##          sub-events it contains)
##        - altered SIZE per component type. For non-CPX records this is
##          abs(SVLEN) added to that record's own SVTYPE. For CPX records,
##          CPX_INTERVALS (e.g. "DEL_chr1:100-200,DUP_chr1:300-400") is
##          parsed and each sub-interval's span is added to ITS OWN type's
##          size total (so a CPX event's size never appears under "CPX" --
##          it's distributed across the DEL/DUP/INV/etc. size buckets).
##
## A single (non-scattered) design is used for tasks 2-3 because SV cohorts
## can have tens of thousands of samples; scattering per sample would mean
## re-reading the whole VCF once per sample.

workflow PerSampleVariantCountSizeSummary {
  input {
    File vcf_gz
    File vcf_tbi
    String output_prefix
    String bcftools_docker = "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    String python_docker = "python:3.11-slim"

    RuntimeAttr? runtime_attr_sample_ids
    RuntimeAttr? runtime_attr_extract
    RuntimeAttr? runtime_attr_summarize
  }

  call ExtractSampleIds {
    input:
      vcf_gz = vcf_gz,
      docker = bcftools_docker,
      runtime_attr_override = runtime_attr_sample_ids
  }

  call ExtractNonrefPassVariants {
    input:
      vcf_gz = vcf_gz,
      vcf_tbi = vcf_tbi,
      output_prefix = output_prefix,
      docker = bcftools_docker,
      runtime_attr_override = runtime_attr_extract
  }

  call SummarizeCountsAndSizes {
    input:
      nonref_variants_table = ExtractNonrefPassVariants.nonref_variants_table,
      sample_ids = ExtractSampleIds.sample_ids,
      output_prefix = output_prefix,
      docker = python_docker,
      runtime_attr_override = runtime_attr_summarize
  }

  output {
    File sample_id_list = ExtractSampleIds.sample_id_list
    File nonref_variants_table = ExtractNonrefPassVariants.nonref_variants_table
    File variant_counts_table = SummarizeCountsAndSizes.variant_counts_table
    File variant_sizes_table = SummarizeCountsAndSizes.variant_sizes_table
  }
}

task ExtractSampleIds {
  input {
    File vcf_gz
    String docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: ceil(10 + size(vcf_gz, "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 2,
    max_retries: 0
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    bcftools query -l ~{vcf_gz} > sample_ids.list
    wc -l < sample_ids.list
  >>>

  output {
    File sample_id_list = "sample_ids.list"
    Array[String] sample_ids = read_lines("sample_ids.list")
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

task ExtractNonrefPassVariants {
  input {
    File vcf_gz
    File vcf_tbi
    String output_prefix
    String docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(20 + size(vcf_gz, "GB") * 3),
    boot_disk_gb: 10,
    preemptible_tries: 2,
    max_retries: 0
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  String out_name = output_prefix + ".nonref_pass_variants.tsv.gz"

  # NOTE: GT="alt" (not GT!="ref") is required here -- GT!="ref" was tried
  # first and found to incorrectly retain missing (./.) genotypes alongside
  # true alt-carrying ones, silently inflating counts by ~80x. GT="alt"
  # matches only genotypes with at least one non-ref, non-missing allele,
  # which is exactly the non-ref/non-missing definition used throughout
  # this pipeline.
  command <<<
    set -euo pipefail

    { echo -e "sample_id\tID\tSVTYPE\tSVLEN\tCPX_INTERVALS"; \
      bcftools query \
        -i 'FILTER="PASS" && GT="alt"' \
        -f '[%SAMPLE\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/CPX_INTERVALS\n]' \
        ~{vcf_gz}; \
    } | gzip -c > ~{out_name}
  >>>

  output {
    File nonref_variants_table = out_name
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

task SummarizeCountsAndSizes {
  input {
    File nonref_variants_table
    Array[String] sample_ids
    String output_prefix
    String docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 8,
    disk_gb: ceil(20 + size(nonref_variants_table, "GB") * 3),
    boot_disk_gb: 10,
    preemptible_tries: 2,
    max_retries: 0
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  String counts_out = output_prefix + ".variant_counts_per_sample.tsv"
  String sizes_out = output_prefix + ".variant_sizes_per_sample.tsv"

  command <<<
    set -euo pipefail

    python3 << 'EOF'
    import gzip

    in_path = "~{nonref_variants_table}"
    counts_path = "~{counts_out}"
    sizes_path = "~{sizes_out}"

    def open_maybe_gz(path):
        if path.endswith(".gz"):
            return gzip.open(path, "rt")
        return open(path, "rt")

    counts = {}
    sizes = {}
    count_types = set()
    size_types = set()

    with open_maybe_gz(in_path) as f:
        header = f.readline()
        for line in f:
            sample_id, variant_id, svtype, svlen_str, cpx_intervals = line.rstrip("\n").split("\t")

            counts.setdefault(sample_id, {})
            counts[sample_id][svtype] = counts[sample_id].get(svtype, 0) + 1
            count_types.add(svtype)

            sizes.setdefault(sample_id, {})
            if svtype == "CPX" and cpx_intervals != ".":
                for interval in cpx_intervals.split(","):
                    comp_type, _, coord = interval.partition("_")
                    _, _, span = coord.rpartition(":")
                    start_str, _, end_str = span.partition("-")
                    if start_str.isdigit() and end_str.isdigit():
                        comp_size = abs(int(end_str) - int(start_str))
                        sizes[sample_id][comp_type] = sizes[sample_id].get(comp_type, 0) + comp_size
                        size_types.add(comp_type)
            elif svlen_str not in (".", ""):
                sizes[sample_id][svtype] = sizes[sample_id].get(svtype, 0) + abs(int(svlen_str))
                size_types.add(svtype)

    all_samples = "~{sep=',' sample_ids}".split(",")
    count_types = sorted(count_types)
    size_types = sorted(size_types)

    with open(counts_path, "w") as out:
        out.write("sample_id\t" + "\t".join(count_types) + "\n")
        for sample_id in all_samples:
            row = counts.get(sample_id, {})
            out.write(sample_id + "\t" + "\t".join(str(row.get(t, 0)) for t in count_types) + "\n")

    with open(sizes_path, "w") as out:
        out.write("sample_id\t" + "\t".join(size_types) + "\n")
        for sample_id in all_samples:
            row = sizes.get(sample_id, {})
            out.write(sample_id + "\t" + "\t".join(str(row.get(t, 0)) for t in size_types) + "\n")

    print(f"n_samples={len(all_samples)} count_types={count_types} size_types={size_types}")
    EOF
  >>>

  output {
    File variant_counts_table = counts_out
    File variant_sizes_table = sizes_out
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
