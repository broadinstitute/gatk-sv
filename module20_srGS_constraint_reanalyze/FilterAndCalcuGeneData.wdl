version 1.0

import "Structs.wdl"
import "PermutateSVAnnotation.wdl" as Permu

# Workflow:
#   For each integrated_file in the input list:
#     1. FilterIntegratedBySvid  — keep only rows whose SVID (col 1) appears in
#                                  the 4th column of a supplied bed.gz
#     2. Task8_CalcuGeneData     — run the per-gene SV calculation R script on
#                                  the filtered file
#
# The bed.gz (sv_bed) is shared across all scattered items.
workflow FilterAndCalcuGeneData {
  input {
    Array[File]  integrated_files   # one per permutation seed
    File         sv_bed             # BED.gz; col 4 = SVID whitelist
    File         sv_info            # passed to Task8
    File         gene_info          # passed to Task8
    File         calcu_r_script     # the calcu.gene.data.reanno.permu.R script
    String       seed_suffix        # e.g. "permu_seed1"  (one per file — see note)

    String       docker             # docker with bedtools / awk
    String       r_docker           # docker with R

    RuntimeAttr? runtime_attr_filter
    RuntimeAttr? runtime_attr_calcu
  }

  scatter (integrated_file in integrated_files) {

    # Step 1: filter integrated_file to SVIDs in sv_bed col 4
    call FilterIntegratedBySvid {
      input:
        integrated_file       = integrated_file,
        sv_bed                = sv_bed,
        docker                = docker,
        runtime_attr_override = runtime_attr_filter
    }

    # Step 2: run Task8 on the filtered file
    call Permu.Task8_CalcuGeneData {
      input:
        r_script         = calcu_r_script,
        seed_suffix      = seed_suffix,
        integrated_file  = FilterIntegratedBySvid.filtered_integrated,
        sv_info          = sv_info,
        gene_info        = gene_info,
        docker           = r_docker,
        runtime_attr_override = runtime_attr_calcu
    }
  }

  output {
    Array[File] filtered_integrated = FilterIntegratedBySvid.filtered_integrated
    Array[File] result_rdata        = Task8_CalcuGeneData.result
    Array[File] result_tsv_gz       = Task8_CalcuGeneData.result_tsv_gz
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Task: filter integrated_file rows by SVID whitelist from a bed.gz col 4
# ─────────────────────────────────────────────────────────────────────────────
task FilterIntegratedBySvid {
  input {
    File    integrated_file   # TSV (possibly .gz); col 1 = SVID
    File    sv_bed            # BED.gz; col 4 = SVID whitelist
    String  docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    mem_gb:            4,
    cpu_cores:         1,
    disk_gb:           20,
    boot_disk_gb:      10,
    preemptible_tries: 3,
    max_retries:       1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  # Derive output name: strip .gz if present, append .filtered
  String base = if sub(basename(integrated_file), "\\.gz$", "") != basename(integrated_file)
                then sub(basename(integrated_file), "\\.gz$", ".filtered.gz")
                else basename(integrated_file) + ".filtered"

  command <<<
    set -euo pipefail

    # Extract SVID whitelist from col 4 of the bed.gz
    zcat ~{sv_bed} | awk -F'\t' '{print $4}' | sort -u > svid_whitelist.txt

    # Determine if integrated_file is gzipped
    if file ~{integrated_file} | grep -q gzip; then
      READER="zcat"
    else
      READER="cat"
    fi

    # Filter: keep header line(s) starting with # or the first line,
    # then keep data rows whose col 1 is in the whitelist
    $READER ~{integrated_file} | awk -F'\t' '
      NR == FNR { whitelist[$1] = 1; next }
      FNR == 1  { print; next }          # always keep header
      $1 in whitelist { print }
    ' svid_whitelist.txt - \
    | gzip -c > ~{base}
  >>>

  output {
    File filtered_integrated = base
  }

  runtime {
    docker:         docker
    memory:         select_first([runtime_attr.mem_gb, 4]) + " GiB"
    cpu:            select_first([runtime_attr.cpu_cores, 1])
    disks:          "local-disk " + select_first([runtime_attr.disk_gb, 20]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, 10])
    preemptible:    select_first([runtime_attr.preemptible_tries, 3])
    maxRetries:     select_first([runtime_attr.max_retries, 1])
  }
}
