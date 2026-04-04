version 1.0

import "Structs.wdl"

# Workflow: for a list of GTFs, compute per-gene proportion covered by
# segmental duplications (SD) and simple repeats (SR), and produce
# integrated gene x sample tables.
workflow GeneTranscriptRepeatOverlap {
  input {
    Array[File]   gtf_files           # one GTF per sample/annotation
    Array[String] column_labels       # one label per GTF (becomes column name)
    File          sd_file             # segmental duplications BED
    File          sr_file             # simple repeats BED
    String        docker              # docker with bedtools + python3
    RuntimeAttr?  runtime_attr_gtf_to_bed
    RuntimeAttr?  runtime_attr_calc_overlap
    RuntimeAttr?  runtime_attr_integrate
  }

  scatter (i in range(length(gtf_files))) {

    # Task 1: GTF -> transcript-level BED
    call GTFTranscriptsToBed {
      input:
        gtf                   = gtf_files[i],
        label                 = column_labels[i],
        docker                = docker,
        runtime_attr_override = runtime_attr_gtf_to_bed
    }

    # Task 2a: proportion of each gene covered by SD
    call CalcRepeatOverlap as CalcSD {
      input:
        transcript_bed        = GTFTranscriptsToBed.transcript_bed,
        repeat_file           = sd_file,
        label                 = column_labels[i],
        repeat_type           = "SD",
        docker                = docker,
        runtime_attr_override = runtime_attr_calc_overlap
    }

    # Task 2b: proportion of each gene covered by SR
    call CalcRepeatOverlap as CalcSR {
      input:
        transcript_bed        = GTFTranscriptsToBed.transcript_bed,
        repeat_file           = sr_file,
        label                 = column_labels[i],
        repeat_type           = "SR",
        docker                = docker,
        runtime_attr_override = runtime_attr_calc_overlap
    }
  }

  # Task 3a: merge SD proportions into one gene x sample table
  call IntegrateOverlapTable as IntegrateSD {
    input:
      overlap_tsvs          = CalcSD.overlap_tsv,
      column_labels         = column_labels,
      output_prefix         = "gene_SD_overlap",
      docker                = docker,
      runtime_attr_override = runtime_attr_integrate
  }

  # Task 3b: merge SR proportions into one gene x sample table
  call IntegrateOverlapTable as IntegrateSR {
    input:
      overlap_tsvs          = CalcSR.overlap_tsv,
      column_labels         = column_labels,
      output_prefix         = "gene_SR_overlap",
      docker                = docker,
      runtime_attr_override = runtime_attr_integrate
  }

  output {
    Array[File] transcript_beds     = GTFTranscriptsToBed.transcript_bed
    Array[File] sd_overlap_tsvs     = CalcSD.overlap_tsv
    Array[File] sr_overlap_tsvs     = CalcSR.overlap_tsv
    File        sd_integrated_table = IntegrateSD.integrated_table
    File        sr_integrated_table = IntegrateSR.integrated_table
  }
}


# ---------------------------------------------------------------------------
# Task 1: Parse GTF → sorted BED6 of transcript boundaries
#   col4 = gene_name, col5 = transcript_id
#   Falls back to exon bounding-box if no 'transcript' feature rows exist.
# ---------------------------------------------------------------------------
task GTFTranscriptsToBed {
  input {
    File         gtf
    String       label
    String       docker
    RuntimeAttr? runtime_attr_override
  }

  String out_bed = label + ".transcripts.bed"

  command <<<
    set -euo pipefail

    python3 <<'PYEOF'
import gzip
import sys
from collections import defaultdict

gtf_file = "~{gtf}"
out_file = "~{out_bed}"

opener = gzip.open if gtf_file.endswith('.gz') else open

records = []

# First pass: collect 'transcript' feature rows
with opener(gtf_file, 'rt') as fin:
    for line in fin:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 9 or fields[2] != 'transcript':
            continue
        chrom     = fields[0]
        start     = int(fields[3]) - 1   # GTF is 1-based closed; BED is 0-based half-open
        end       = int(fields[4])
        strand    = fields[6]
        attr_str  = fields[8]

        gene_name = ''
        tid       = ''
        for attr in attr_str.split(';'):
            attr = attr.strip()
            if attr.startswith('gene_name'):
                parts = attr.split('"')
                if len(parts) >= 2:
                    gene_name = parts[1]
            elif attr.startswith('transcript_id'):
                parts = attr.split('"')
                if len(parts) >= 2:
                    tid = parts[1]
        if gene_name and tid:
            records.append((chrom, start, end, gene_name, tid, strand))

# Fallback: derive transcript bounding boxes from exon entries
if not records:
    sys.stderr.write("Warning: no 'transcript' features found in GTF; inferring from exons.\n")
    bounds = {}  # tid -> [chrom, min_start, max_end, gene_name, strand]
    with opener(gtf_file, 'rt') as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'exon':
                continue
            chrom    = fields[0]
            start    = int(fields[3]) - 1
            end      = int(fields[4])
            strand   = fields[6]
            attr_str = fields[8]
            gene_name = ''
            tid       = ''
            for attr in attr_str.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_name'):
                    parts = attr.split('"')
                    if len(parts) >= 2:
                        gene_name = parts[1]
                elif attr.startswith('transcript_id'):
                    parts = attr.split('"')
                    if len(parts) >= 2:
                        tid = parts[1]
            if gene_name and tid:
                if tid not in bounds:
                    bounds[tid] = [chrom, start, end, gene_name, strand]
                else:
                    bounds[tid][1] = min(bounds[tid][1], start)
                    bounds[tid][2] = max(bounds[tid][2], end)
    for tid, (chrom, s, e, gene, strand) in bounds.items():
        records.append((chrom, s, e, gene, tid, strand))

# Sort by chrom then start
records.sort(key=lambda r: (r[0], r[1]))

with open(out_file, 'w') as fout:
    for chrom, start, end, gene_name, tid, strand in records:
        fout.write(f"{chrom}\t{start}\t{end}\t{gene_name}\t{tid}\t{strand}\n")
PYEOF
  >>>

  output {
    File transcript_bed = out_bed
  }

  RuntimeAttr default_attr = object {
    cpu_cores:         2,
    mem_gb:            4,
    disk_gb:           20 + ceil(size(gtf, "GiB") * 4),
    boot_disk_gb:      10,
    preemptible_tries: 3,
    max_retries:       1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu:            select_first([runtime_attr.cpu_cores,          default_attr.cpu_cores])
    memory:         select_first([runtime_attr.mem_gb,             default_attr.mem_gb]) + " GiB"
    disks:          "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    docker:         docker
    preemptible:    select_first([runtime_attr.preemptible_tries,  default_attr.preemptible_tries])
    maxRetries:     select_first([runtime_attr.max_retries,        default_attr.max_retries])
  }
}


# ---------------------------------------------------------------------------
# Task 2: For each transcript BED, compute proportion of each gene covered
#         by a repeat BED (SD or SR).  Repeat intervals are merged first to
#         avoid double-counting overlapping entries.
#
# Output TSV: gene <tab> <label>
# ---------------------------------------------------------------------------
task CalcRepeatOverlap {
  input {
    File         transcript_bed
    File         repeat_file
    String       label
    String       repeat_type    # "SD" or "SR" — used only in output filename
    String       docker
    RuntimeAttr? runtime_attr_override
  }

  String out_tsv = label + "." + repeat_type + ".overlap.tsv"

  command <<<
    set -euo pipefail

    # Merge overlapping repeat intervals to avoid double-counting
    bedtools sort -i "~{repeat_file}" | bedtools merge -i - > repeat_merged.bed

    # bedtools coverage output (BED6 input → 10 columns):
    #   0:chrom 1:start 2:end 3:gene_name 4:tid 5:strand
    #   6:#features  7:covered_bases  8:transcript_len  9:fraction
    bedtools coverage \
      -a "~{transcript_bed}" \
      -b repeat_merged.bed \
      > coverage_raw.tsv

    python3 <<'PYEOF'
from collections import defaultdict

gene_covered = defaultdict(int)
gene_total   = defaultdict(int)

with open('coverage_raw.tsv') as f:
    for line in f:
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 10:
            continue
        gene_name     = fields[3]
        covered_bases = int(fields[7])
        total_length  = int(fields[8])
        gene_covered[gene_name] += covered_bases
        gene_total[gene_name]   += total_length

label   = "~{label}"
out_tsv = "~{out_tsv}"

with open(out_tsv, 'w') as out:
    out.write(f"gene\t{label}\n")
    for gene in sorted(gene_covered):
        covered = gene_covered[gene]
        total   = gene_total[gene]
        prop    = covered / total if total > 0 else 0.0
        out.write(f"{gene}\t{prop:.6f}\n")
PYEOF
  >>>

  output {
    File overlap_tsv = out_tsv
  }

  RuntimeAttr default_attr = object {
    cpu_cores:         2,
    mem_gb:            8,
    disk_gb:           20 + ceil((size(transcript_bed, "GiB") + size(repeat_file, "GiB")) * 4),
    boot_disk_gb:      10,
    preemptible_tries: 3,
    max_retries:       1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu:            select_first([runtime_attr.cpu_cores,          default_attr.cpu_cores])
    memory:         select_first([runtime_attr.mem_gb,             default_attr.mem_gb]) + " GiB"
    disks:          "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    docker:         docker
    preemptible:    select_first([runtime_attr.preemptible_tries,  default_attr.preemptible_tries])
    maxRetries:     select_first([runtime_attr.max_retries,        default_attr.max_retries])
  }
}


# ---------------------------------------------------------------------------
# Task 3: Merge per-sample overlap TSVs into one wide table.
#   Rows = genes (union across all inputs)
#   Cols = gene, <column_labels[0]>, <column_labels[1]>, ...
#   Missing values filled with NA.
# ---------------------------------------------------------------------------
task IntegrateOverlapTable {
  input {
    Array[File]   overlap_tsvs
    Array[String] column_labels
    String        output_prefix   # e.g. "gene_SD_overlap" → gene_SD_overlap.tsv
    String        docker
    RuntimeAttr?  runtime_attr_override
  }

  String out_tsv = output_prefix + ".tsv"

  command <<<
    set -euo pipefail

    python3 <<'PYEOF'
tsv_files  = "~{sep=',' overlap_tsvs}".split(',')
col_labels = "~{sep=',' column_labels}".split(',')
out_file   = "~{out_tsv}"

all_genes = set()
tables = []

for fp in tsv_files:
    d = {}
    with open(fp) as fh:
        fh.readline()  # skip header line
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                d[parts[0]] = parts[1]
                all_genes.add(parts[0])
    tables.append(d)

with open(out_file, 'w') as out:
    out.write('gene\t' + '\t'.join(col_labels) + '\n')
    for gene in sorted(all_genes):
        row = [gene] + [d.get(gene, 'NA') for d in tables]
        out.write('\t'.join(row) + '\n')
PYEOF
  >>>

  output {
    File integrated_table = out_tsv
  }

  RuntimeAttr default_attr = object {
    cpu_cores:         2,
    mem_gb:            8,
    disk_gb:           20,
    boot_disk_gb:      10,
    preemptible_tries: 3,
    max_retries:       1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
    cpu:            select_first([runtime_attr.cpu_cores,          default_attr.cpu_cores])
    memory:         select_first([runtime_attr.mem_gb,             default_attr.mem_gb]) + " GiB"
    disks:          "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    docker:         docker
    preemptible:    select_first([runtime_attr.preemptible_tries,  default_attr.preemptible_tries])
    maxRetries:     select_first([runtime_attr.max_retries,        default_attr.max_retries])
  }
}
