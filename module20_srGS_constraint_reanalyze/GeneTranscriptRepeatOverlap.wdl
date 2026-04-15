version 1.0

import "Structs.wdl"
import "RDataToTsv.wdl" as RDataToTsv

# Workflow: for a list of GTFs, compute per-gene proportion covered by
# segmental duplications (SD) and simple repeats (SR), and produce
# integrated gene x sample tables.
workflow GeneTranscriptRepeatOverlap {
  input {
    Array[File]   gtf_files           # one GTF per sample/annotation (same order as rdata_files)
    Array[File]   rdata_files         # one .rData per sample, same order as gtf_files
    String?        rdata_obj_name = "gene.data.reanno.permu"  # R object name inside each .rData file
    Array[String] column_labels       # one label per GTF/rData (becomes column name)
    File          sd_file             # segmental duplications BED
    File          sr_file             # simple repeats BED
    String        docker              # docker with bedtools + python3
    String        r_docker            # docker with R installed (for RData conversion)
    RuntimeAttr?  runtime_attr_gtf_to_bed
    RuntimeAttr?  runtime_attr_calc_overlap
    RuntimeAttr?  runtime_attr_integrate
    RuntimeAttr?  runtime_attr_rdata_to_tsv
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

    # Task 1b: convert rData to TSV.gz (same index as GTF)
    call RDataToTsv.ConvertRDataToTsv {
      input:
        rdata_file            = rdata_files[i],
        rdata_obj_name        = rdata_obj_name,
        r_docker              = r_docker,
        runtime_attr_override = runtime_attr_rdata_to_tsv
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
    Array[File] rdata_tsv_gz        = ConvertRDataToTsv.tsv_gz
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

  command <<<
    set -euo pipefail


    zcat  ~{gtf} \
     | awk '{if ($3=="transcript") print}' \
     | cut -f1,4,5,9 \
     | awk '{print $1,$2,$3,$5, $9}' \
     | sed -e 's/"//g' | sed -e 's/;//g' | sed -e 's/ /\t/g' \
     > "~{label}.transcripts.bed"
  >>>

  output {
    File transcript_bed = "~{label}.transcripts.bed"
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

    awk '{print $5,$NF}' coverage_raw.tsv | sed -e 's/ /\t/g' >  "~{label}.~{repeat_type}.overlap.tsv"
  >>>

  output {
    File overlap_tsv = "~{label}.~{repeat_type}.overlap.tsv"
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
# Task 3: Merge per-sample overlap TSVs into one wide bgzipped table.
#   Rows = genes (union across all inputs)
#   Cols = gene, <column_labels[0]>, <column_labels[1]>, ...
#   Missing values filled with NA.
# ---------------------------------------------------------------------------
task IntegrateOverlapTable {
  input {
    Array[File]   overlap_tsvs
    Array[String] column_labels
    String        output_prefix   # e.g. "gene_SD_overlap" → gene_SD_overlap.tsv.gz
    String        docker
    RuntimeAttr?  runtime_attr_override
  }

  String out_tsv = output_prefix + ".tsv"
  String out_tsv_gz = output_prefix + ".tsv.gz"

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

    bgzip -c "~{out_tsv}" > "~{out_tsv_gz}"
  >>>

  output {
    File integrated_table = out_tsv_gz
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
