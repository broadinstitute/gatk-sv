version 1.0

import "Structs.wdl"

workflow VcfToBedBySvlenShards {
  input {
    File input_vcf  # Accepts VCF.gz, BCF, or uncompressed VCF. Format is auto-detected from file extension.
    File vcf_to_bed_script

    Int variants_per_shard = 100000
    String output_prefix = "vcf_to_bed"
    String python_docker = "python:3.11-slim"

    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_convert_split
    RuntimeAttr? runtime_attr_concat
  }

  call SplitVcfByVariantCount {
    input:
      input_vcf             = input_vcf,
      output_prefix         = output_prefix,
      variants_per_shard    = variants_per_shard,
      python_docker         = python_docker,
      runtime_attr_override = runtime_attr_split_vcf
  }

  scatter (shard_vcf in SplitVcfByVariantCount.shard_vcfs) {
    call ConvertAndSplitBedBySvlen {
      input:
        shard_vcf              = shard_vcf,
        vcf_to_bed_script      = vcf_to_bed_script,
        python_docker          = python_docker,
        runtime_attr_override  = runtime_attr_convert_split
    }
  }

  call ConcatBedGz as concat_lt50 {
    input:
      bed_files              = ConvertAndSplitBedBySvlen.bed_lt50,
      output_file            = output_prefix + ".svlen_lt50.bed.gz",
      python_docker          = python_docker,
      runtime_attr_override  = runtime_attr_concat
  }

  call ConcatBedGz as concat_ge50 {
    input:
      bed_files              = ConvertAndSplitBedBySvlen.bed_ge50,
      output_file            = output_prefix + ".svlen_ge50.bed.gz",
      python_docker          = python_docker,
      runtime_attr_override  = runtime_attr_concat
  }

  output {
    Array[File] shard_vcfs = SplitVcfByVariantCount.shard_vcfs

    Array[File] shard_all_beds = ConvertAndSplitBedBySvlen.bed_all
    Array[File] shard_lt50_beds = ConvertAndSplitBedBySvlen.bed_lt50
    Array[File] shard_ge50_beds = ConvertAndSplitBedBySvlen.bed_ge50

    File merged_lt50_bed = concat_lt50.merged_bed
    File merged_ge50_bed = concat_ge50.merged_bed
  }
}


task SplitVcfByVariantCount {
  input {
    File input_vcf
    String output_prefix
    Int variants_per_shard
    String python_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 2,
    mem_gb: 8,
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    # Install bcftools if needed (for BCF format support)
    which bcftools > /dev/null 2>&1 || (apt-get update && apt-get install -y bcftools)
    
    set -euo pipefail

    python3 << 'PY'
import gzip
import subprocess
import os

input_vcf = "~{input_vcf}"
prefix = "~{output_prefix}"
chunk_size = int("~{variants_per_shard}")

# Detect format based on file extension
input_format = "unknown"
if input_vcf.endswith('.bcf'):
    input_format = "bcf"
elif input_vcf.endswith('.vcf.gz'):
    input_format = "vcf.gz"
elif input_vcf.endswith('.vcf'):
    input_format = "vcf"

def open_text(path, mode, file_format):
    """Open file in text mode, handling gzip and bcf formats"""
    if file_format == "bcf":
        # Use bcftools to convert BCF to uncompressed VCF stream
        proc = subprocess.Popen(['bcftools', 'view', path], stdout=subprocess.PIPE, text=True)
        return proc.stdout
    elif file_format == "vcf.gz":
        return gzip.open(path, mode)
    else:  # vcf or unknown
        return open(path, mode)

header_lines = []
shard_idx = 0
records_in_shard = 0
records_total = 0
out = None

with open_text(input_vcf, 'rt', input_format) as fh:
    for line in fh:
        if line.startswith('#'):
            header_lines.append(line)
            continue

        if out is None or records_in_shard >= chunk_size:
            if out is not None:
                out.close()
            shard_path = f"{prefix}.shard_{shard_idx:05d}.vcf.gz"
            out = gzip.open(shard_path, 'wt')
            for h in header_lines:
                out.write(h)
            shard_idx += 1
            records_in_shard = 0

        out.write(line)
        records_in_shard += 1
        records_total += 1

if out is not None:
    out.close()
elif records_total == 0:
    # Empty VCF body: still emit one header-only shard.
    shard_path = f"{prefix}.shard_{0:05d}.vcf.gz"
    with gzip.open(shard_path, 'wt') as out_empty:
        for h in header_lines:
            out_empty.write(h)

print({"records_total": records_total, "shards": max(1, shard_idx), "input_format": input_format})
PY
  >>>

  output {
    Array[File] shard_vcfs = glob("*.shard_*.vcf.gz")
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


task ConvertAndSplitBedBySvlen {
  input {
    File shard_vcf
    File vcf_to_bed_script
    String python_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 2,
    mem_gb: 8,
    disk_gb: 60,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String shard_label = basename(shard_vcf, ".vcf.gz")
  String bed_all_out = shard_label + ".all_info.bed.gz"
  String bed_lt50_out = shard_label + ".svlen_lt50.bed.gz"
  String bed_ge50_out = shard_label + ".svlen_ge50.bed.gz"

  command <<<
    set -euo pipefail

    python3 ~{vcf_to_bed_script} \
      --input ~{shard_vcf} \
      --output ~{bed_all_out}

    python3 << 'PY'
import gzip

in_bed = "~{bed_all_out}"
lt50_out = "~{bed_lt50_out}"
ge50_out = "~{bed_ge50_out}"

lt50 = 0
ge50 = 0
unknown = 0

def parse_svlen(raw):
    if raw is None or raw == '.' or raw == '':
        return None
    first = raw.split(',', 1)[0].strip()
    if first == '':
        return None
    try:
        return abs(int(float(first)))
    except ValueError:
        return None

with gzip.open(in_bed, 'rt') as fin, \
     gzip.open(lt50_out, 'wt') as flt, \
     gzip.open(ge50_out, 'wt') as fge:

    header = fin.readline()
    if not header:
        raise RuntimeError('Input BED has no header')

    cols = header.rstrip('\n').split('\t')
    if 'SVLEN' not in cols:
        raise RuntimeError('SVLEN column not found in BED header')

    idx = cols.index('SVLEN')
    flt.write(header)
    fge.write(header)

    for line in fin:
        if not line.strip():
            continue
        fields = line.rstrip('\n').split('\t')
        svlen = parse_svlen(fields[idx] if idx < len(fields) else None)
        if svlen is None:
            unknown += 1
            continue
        if svlen < 50:
            flt.write(line)
            lt50 += 1
        else:
            fge.write(line)
            ge50 += 1

print({"lt50": lt50, "ge50": ge50, "unknown": unknown})
PY
  >>>

  output {
    File bed_all = bed_all_out
    File bed_lt50 = bed_lt50_out
    File bed_ge50 = bed_ge50_out
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


task ConcatBedGz {
  input {
    Array[File] bed_files
    String output_file
    String python_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: 30,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python3 << 'PY'
import gzip

inputs = "~{sep=' ' bed_files}".split()
out_path = "~{output_file}"

if not inputs:
    raise RuntimeError('No BED files to concatenate')

header_written = False
rows = 0

with gzip.open(out_path, 'wt') as out:
    for path in inputs:
        with gzip.open(path, 'rt') as fin:
            header = fin.readline()
            if not header:
                continue
            if not header_written:
                out.write(header)
                header_written = True
            for line in fin:
                if line.strip():
                    out.write(line)
                    rows += 1

print({"inputs": len(inputs), "rows": rows})
PY
  >>>

  output {
    File merged_bed = output_file
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
