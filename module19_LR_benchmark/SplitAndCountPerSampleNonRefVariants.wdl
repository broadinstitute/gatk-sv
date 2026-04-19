version 1.0

import "Structs.wdl"

workflow SplitAndCountPerSampleNonRefVariants {
  input {
    File input_vcf

    String output_prefix = "per_sample_nonref"
    String bcftools_docker = "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    String python_docker = "python:3.11-slim"

    RuntimeAttr? runtime_attr_sample_list
    RuntimeAttr? runtime_attr_extract
    RuntimeAttr? runtime_attr_split_count
    RuntimeAttr? runtime_attr_merge
  }

  call GetVcfSampleList {
    input:
      input_vcf             = input_vcf,
      bcftools_docker       = bcftools_docker,
      runtime_attr_override = runtime_attr_sample_list
  }

  scatter (sample_name in GetVcfSampleList.sample_names) {
    call ExtractSampleNonRefVcf {
      input:
        input_vcf             = input_vcf,
        sample_name           = sample_name,
        bcftools_docker       = bcftools_docker,
        runtime_attr_override = runtime_attr_extract
    }

    call SplitAndCountSampleVariants {
      input:
        sample_name            = sample_name,
        sample_nonref_vcf      = ExtractSampleNonRefVcf.sample_nonref_vcf,
        python_docker          = python_docker,
        runtime_attr_override  = runtime_attr_split_count
    }
  }

  call MergeSampleCountTables {
    input:
      sample_count_tables    = SplitAndCountSampleVariants.count_tsv,
      merged_output          = output_prefix + ".sample_category_counts.tsv",
      python_docker          = python_docker,
      runtime_attr_override  = runtime_attr_merge
  }

  output {
    File sample_list = GetVcfSampleList.sample_list
    Array[String] sample_names = GetVcfSampleList.sample_names

    Array[File] sample_nonref_vcfs = ExtractSampleNonRefVcf.sample_nonref_vcf
    Array[File] sample_nonref_vcf_indexes = ExtractSampleNonRefVcf.sample_nonref_vcf_tbi

    Array[File] per_sample_count_tables = SplitAndCountSampleVariants.count_tsv
    Array[Array[File]] per_sample_category_vcfs = SplitAndCountSampleVariants.category_vcfs

    File merged_sample_count_table = MergeSampleCountTables.merged_count_table
  }
}


task GetVcfSampleList {
  input {
    File input_vcf
    String bcftools_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: 30,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    bcftools query -l ~{input_vcf} > sample_list.txt
  >>>

  output {
    File sample_list = "sample_list.txt"
    Array[String] sample_names = read_lines("sample_list.txt")
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


task ExtractSampleNonRefVcf {
  input {
    File input_vcf
    String sample_name
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

  String sample_nonref_vcf_name = sample_name + ".nonref.vcf.gz"

  command <<<
    set -euo pipefail

    # Keep only sites where selected sample carries any ALT allele.
    bcftools view \
      -s ~{sample_name} \
      -i 'GT~"[1-9]"' \
      -Oz \
      -o ~{sample_nonref_vcf_name} \
      ~{input_vcf}

    bcftools index -t ~{sample_nonref_vcf_name}
  >>>

  output {
    File sample_nonref_vcf = sample_nonref_vcf_name
    File sample_nonref_vcf_tbi = sample_nonref_vcf_name + ".tbi"
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


task SplitAndCountSampleVariants {
  input {
    String sample_name
    File sample_nonref_vcf
    String python_docker
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

  String count_tsv_name = sample_name + ".variant_category_counts.tsv"

  command <<<
    set -euo pipefail

    python3 << 'PY'
import csv
import gzip
import re
from collections import OrderedDict

sample_name = "~{sample_name}"
input_vcf = "~{sample_nonref_vcf}"
count_tsv = "~{count_tsv_name}"

# Categories requested by user.
categories = [
    "non_trv_snv",
    "non_trv_del_1_49bp",
    "non_trv_ins_1_49bp",
    "non_trv_del_gt49bp",
    "non_trv_ins_gt49bp",
    "non_trid_snv",
    "non_trid_del_1_49bp",
    "non_trid_ins_1_49bp",
    "non_trid_del_gt49bp",
    "non_trid_ins_gt49bp",
    "trv_repeat_unit_le6bp",
    "trv_repeat_unit_gt6bp",
]

extra_audit_categories = [
    "trv_repeat_unit_unknown",
    "uncategorized_non_trv",
    "uncategorized_non_trid",
]

all_categories = categories + extra_audit_categories
counts = OrderedDict((k, 0) for k in all_categories)

# One split VCF per category.
vcf_writers = {}
vcf_handles = {}

def info_to_dict(info_str):
    d = {}
    for token in info_str.split(';'):
        if not token:
            continue
        if '=' in token:
            k, v = token.split('=', 1)
            d[k] = v
        else:
            d[token] = True
    return d


def get_first_numeric(value):
    if value is None:
        return None
    token = str(value).split(',', 1)[0].strip()
    if token == '' or token == '.':
        return None
    try:
        return abs(int(float(token)))
    except ValueError:
        return None


def is_simple_allele(allele):
    return allele not in ('.', '*') and not allele.startswith('<') and '[' not in allele and ']' not in allele


def classify_core(ref, alts, info):
    if len(alts) == 0:
        return None

    if all(len(ref) == 1 and len(a) == 1 and is_simple_allele(a) for a in alts):
        return "snv"

    svtype = str(info.get("SVTYPE", "")).upper().split(',', 1)[0]
    svlen = get_first_numeric(info.get("SVLEN"))

    inferred_type = None
    if svtype in ("DEL", "INS"):
        inferred_type = svtype

    if inferred_type is None:
        if len(alts) == 1 and is_simple_allele(alts[0]) and is_simple_allele(ref):
            if len(alts[0]) > len(ref):
                inferred_type = "INS"
                if svlen is None:
                    svlen = len(alts[0]) - len(ref)
            elif len(ref) > len(alts[0]):
                inferred_type = "DEL"
                if svlen is None:
                    svlen = len(ref) - len(alts[0])

    if inferred_type not in ("DEL", "INS"):
        return None

    if svlen is None or svlen <= 0:
        return None

    size_label = "1_49bp" if svlen <= 49 else "gt49bp"
    return inferred_type.lower() + "_" + size_label


def repeat_unit_size(info):
    motifs = info.get("MOTIFS")
    if motifs:
        for m in str(motifs).split(','):
            m = m.strip()
            if m:
                return len(m)

    ru = info.get("RU")
    if ru:
        for m in str(ru).split(','):
            m = m.strip()
            if m:
                return len(m)

    period = get_first_numeric(info.get("PERIOD"))
    if period is not None:
        return period

    return None


def write_record(category, line):
    counts[category] += 1
    vcf_writers[category].write(line)


header_lines = []

with gzip.open(input_vcf, "rt") as fin:
    for line in fin:
        if line.startswith('#'):
            header_lines.append(line)
            continue

        if not vcf_writers:
            for cat in all_categories:
                path = f"{sample_name}.{cat}.vcf.gz"
                fh = gzip.open(path, "wt")
                vcf_handles[cat] = fh
                vcf_writers[cat] = fh
                for h in header_lines:
                    fh.write(h)

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            continue

        vid = fields[2]
        ref = fields[3]
        alts = fields[4].split(',') if fields[4] else []
        info = info_to_dict(fields[7])

        svid = str(info.get("SVID", vid))
        has_trv_in_svid = "TRV" in svid
        has_trid = "TRID" in info

        core = classify_core(ref, alts, info)

        # Group 1: variants without TRV in SVID.
        if not has_trv_in_svid:
            if core is None:
                write_record("uncategorized_non_trv", line)
            else:
                write_record("non_trv_" + core, line)
        else:
            # Group 3: variants with TRV in SVID split by repeat unit size.
            ru = repeat_unit_size(info)
            if ru is None:
                write_record("trv_repeat_unit_unknown", line)
            elif ru <= 6:
                write_record("trv_repeat_unit_le6bp", line)
            else:
                write_record("trv_repeat_unit_gt6bp", line)

        # Group 2: variants without TRID in INFO, same SNV/INDEL subcategories.
        if not has_trid:
            if core is None:
                write_record("uncategorized_non_trid", line)
            else:
                write_record("non_trid_" + core, line)

# Emit empty category VCFs if no records were seen.
if not vcf_writers:
    for cat in all_categories:
        path = f"{sample_name}.{cat}.vcf.gz"
        fh = gzip.open(path, "wt")
        vcf_handles[cat] = fh
        for h in header_lines:
            fh.write(h)
        fh.close()
else:
    for fh in vcf_handles.values():
        fh.close()

with open(count_tsv, "wt", newline="") as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(["sample"] + all_categories)
    writer.writerow([sample_name] + [counts[k] for k in all_categories])

print({"sample": sample_name, "records_by_category": counts})
PY
  >>>

  output {
    File count_tsv = count_tsv_name
    Array[File] category_vcfs = glob("*.non_tr*.vcf.gz") + glob("*.trv_*.vcf.gz") + glob("*.uncategorized_*.vcf.gz")
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


task MergeSampleCountTables {
  input {
    Array[File] sample_count_tables
    String merged_output
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

inputs = "~{sep=' ' sample_count_tables}".split()
out_path = "~{merged_output}"

if not inputs:
    raise RuntimeError("No sample count tables provided")

header = None
rows_written = 0

with open(out_path, "wt", newline="") as out:
    writer = None
    for path in inputs:
        with open(path, "rt", newline="") as fh:
            reader = csv.reader(fh, delimiter='\t')
            this_header = next(reader)
            if header is None:
                header = this_header
                writer = csv.writer(out, delimiter='\t')
                writer.writerow(header)
            elif this_header != header:
                raise RuntimeError(f"Header mismatch in {path}")

            for row in reader:
                writer.writerow(row)
                rows_written += 1

print({"inputs": len(inputs), "rows": rows_written, "output": out_path})
PY
  >>>

  output {
    File merged_count_table = merged_output
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
