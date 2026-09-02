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
    String docker = "python:3.11-slim"

    RuntimeAttr? runtime_attr_sample_ids
    RuntimeAttr? runtime_attr_extract
    RuntimeAttr? runtime_attr_summarize
  }

  call ExtractSampleIds {
    input:
      vcf_gz = vcf_gz,
      docker = docker,
      runtime_attr_override = runtime_attr_sample_ids
  }

  call ExtractNonrefPassVariants {
    input:
      vcf_gz = vcf_gz,
      output_prefix = output_prefix,
      docker = docker,
      runtime_attr_override = runtime_attr_extract
  }

  call SummarizeCountsAndSizes {
    input:
      nonref_variants_table = ExtractNonrefPassVariants.nonref_variants_table,
      sample_ids = ExtractSampleIds.sample_ids,
      output_prefix = output_prefix,
      docker = docker,
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

    python3 << 'EOF'
    import gzip

    with gzip.open("~{vcf_gz}", "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                samples = fields[9:]
                break
        else:
            raise RuntimeError("No #CHROM header line found in VCF")

    with open("sample_ids.list", "w") as out:
        for s in samples:
            out.write(s + "\n")

    print(f"n_samples={len(samples)}")
    EOF
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
    String output_prefix
    String docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 8,
    disk_gb: ceil(20 + size(vcf_gz, "GB") * 3),
    boot_disk_gb: 10,
    preemptible_tries: 2,
    max_retries: 0
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  String out_name = output_prefix + ".nonref_pass_variants.tsv.gz"

  command <<<
    set -euo pipefail

    python3 << 'EOF'
    import gzip

    in_path = "~{vcf_gz}"
    out_path = "~{out_name}"

    def parse_info(info_str):
        info = {}
        for entry in info_str.split(";"):
            if "=" in entry:
                key, _, value = entry.partition("=")
                info[key] = value
        return info

    def is_nonref(gt_str):
        gt_str = gt_str.split(":")[0]
        alleles = gt_str.replace("|", "/").split("/")
        return any(a not in (".", "0") for a in alleles)

    samples = None
    gt_idx = 0
    n_out = 0

    with gzip.open(in_path, "rt") as fin, gzip.open(out_path, "wt") as fout:
        fout.write("sample_id\tID\tSVTYPE\tSVLEN\tCPX_INTERVALS\n")
        for line in fin:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                samples = fields[9:]
                continue

            fields = line.rstrip("\n").split("\t")
            variant_id = fields[2]
            filter_val = fields[6]
            info = parse_info(fields[7])
            format_keys = fields[8].split(":")
            gt_idx = format_keys.index("GT")

            if filter_val != "PASS":
                continue

            svtype = info.get("SVTYPE", ".")
            svlen = info.get("SVLEN", ".")
            cpx_intervals = info.get("CPX_INTERVALS", ".")

            for sample, sample_field in zip(samples, fields[9:]):
                sample_gt = sample_field.split(":")[gt_idx]
                if is_nonref(sample_gt):
                    fout.write(f"{sample}\t{variant_id}\t{svtype}\t{svlen}\t{cpx_intervals}\n")
                    n_out += 1

    print(f"n_nonref_pass_records_written={n_out}")
    EOF
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
