version 1.0

# Reads two variant files (VCF or BED format), finds sample IDs common to both,
# then for each common sample:
#   1. Extracts non-ref variants in two schemas (query / ref) from each file.
#   2. Compares the two callsets in both directions with compare_callsets_V2.sh.
#   3. Annotates each comparison BED with genomic context (repeat, seg-dup,
#      simple-repeat) using annotate_genomic_context.sh.
#   4. Appends a sample_id column to every annotated BED.
#   5. Concatenates all per-sample, per-direction results into one merged BED.
#
# Supported variant-file formats (plain or gzip-compressed):
#
#   VCF / BCF  (.vcf, .vcf.gz, .bcf)
#     - Sample IDs: genotype columns in the VCF header.
#     - Non-ref: GT not ref-hom and not missing.
#     - Coordinates: start = POS-1 (0-based), end = INFO/END.
#     - SVTYPE, SVLEN, AF from INFO field; SVLEN stored as absolute value.
#
#   BED  (.bed, .bed.gz)  — layout matches example.bed in this directory:
#       #chrom  start  end  name  svtype  samples  SVTYPE  SVLEN  AF
#     - Sample IDs: comma-separated values in the `bed_sample_col` column.
#     - Non-ref: rows where the sample appears in the samples column.
#
# Per-sample extraction output schemas:
#   query  : #chrom  start  end  name  SVTYPE  SVLEN
#   ref    : #chrom  start  end  VID   svtype  length  AF  samples
#
# compare_callsets_V2.sh output schema (input to annotation step):
#   #chr  start  end  VID  svtype  length  AF  ovr1a  ovr1b  ovr2a  ovr2b  ovr3
#
# Final merged output adds columns: GC  sample_id

workflow ExtractCommonSamples {
  input {
    File variant_file_1
    File variant_file_2

    # Column name in the BED header holding comma-separated sample IDs.
    String bed_sample_col = "samples"

    # Short labels used in output filenames (e.g. caller names or cohort IDs).
    String file_label_1 = "file1"
    String file_label_2 = "file2"

    String output_prefix = "common_samples"

    # Scripts — both compare_callsets_V2.sh and its companion R script must be
    # provided; they are co-localised in the task working directory at runtime.
    File compare_callsets_script          # compare_callsets_V2.sh
    File companion_r_script              # R script called by compare_callsets_V2.sh
    File annotate_genomic_context_script  # annotate_genomic_context.sh

    # hg38 reference annotation BEDs consumed by annotate_genomic_context.sh.
    File hg38_repeat_mask_bed
    File hg38_seg_dup_bed
    File hg38_simp_rep_bed

    # Docker for extraction steps (needs bcftools + bgzip).
    String bcftools_docker = "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    # Docker for compare / annotate steps (needs bash, R, bgzip, awk, bedtools).
    String scripts_docker
  }

  call ExtractSamples as ExtractSamples1 {
    input:
      variant_file    = variant_file_1,
      bed_sample_col  = bed_sample_col,
      output_prefix   = "file1",
      bcftools_docker = bcftools_docker
  }

  call ExtractSamples as ExtractSamples2 {
    input:
      variant_file    = variant_file_2,
      bed_sample_col  = bed_sample_col,
      output_prefix   = "file2",
      bcftools_docker = bcftools_docker
  }

  call FindCommonSamples {
    input:
      samples_file_1  = ExtractSamples1.samples_file,
      samples_file_2  = ExtractSamples2.samples_file,
      output_prefix   = output_prefix,
      bcftools_docker = bcftools_docker
  }

  Array[String] common_sample_ids = read_lines(FindCommonSamples.common_samples_file)

  scatter (sample_id in common_sample_ids) {

    # Extract non-ref variants for this sample from each file.
    call ExtractVariantsForSample as ExtractFromFile1 {
      input:
        variant_file    = variant_file_1,
        sample_id       = sample_id,
        bed_sample_col  = bed_sample_col,
        file_label      = file_label_1,
        bcftools_docker = bcftools_docker
    }
    call ExtractVariantsForSample as ExtractFromFile2 {
      input:
        variant_file    = variant_file_2,
        sample_id       = sample_id,
        bed_sample_col  = bed_sample_col,
        file_label      = file_label_2,
        bcftools_docker = bcftools_docker
    }

    # Compare both directions, annotate genomic context, add sample_id column.
    call CompareAndAnnotateSample {
      input:
        file1_query                      = ExtractFromFile1.query_file,
        file1_ref                        = ExtractFromFile1.ref_file,
        file2_query                      = ExtractFromFile2.query_file,
        file2_ref                        = ExtractFromFile2.ref_file,
        sample_id                        = sample_id,
        file_label_1                     = file_label_1,
        file_label_2                     = file_label_2,
        compare_callsets_script          = compare_callsets_script,
        companion_r_script               = companion_r_script,
        annotate_genomic_context_script  = annotate_genomic_context_script,
        hg38_repeat_mask_bed             = hg38_repeat_mask_bed,
        hg38_seg_dup_bed                 = hg38_seg_dup_bed,
        hg38_simp_rep_bed                = hg38_simp_rep_bed,
        scripts_docker                   = scripts_docker
    }

    # Compute per-sample summary statistics from both annotated comparison BEDs.
    call ComputeSampleStats {
      input:
        annotated_bed_1v2 = CompareAndAnnotateSample.annotated_bed_1v2,
        annotated_bed_2v1 = CompareAndAnnotateSample.annotated_bed_2v1,
        sample_id         = sample_id,
        bcftools_docker   = bcftools_docker
    }
  }

  # Concatenate all per-sample annotated BEDs (both directions) into one file.
  call ConcatResults {
    input:
      annotated_beds  = flatten([
                          CompareAndAnnotateSample.annotated_bed_1v2,
                          CompareAndAnnotateSample.annotated_bed_2v1
                        ]),
      output_prefix   = output_prefix,
      bcftools_docker = bcftools_docker
  }

  # Assemble per-sample stat rows into two summary tables (one per direction).
  call ConcatStats {
    input:
      stat_rows_1v2   = ComputeSampleStats.stats_row_1v2,
      stat_rows_2v1   = ComputeSampleStats.stats_row_2v1,
      file_label_1    = file_label_1,
      file_label_2    = file_label_2,
      output_prefix   = output_prefix,
      bcftools_docker = bcftools_docker
  }

  output {
    # Sample lists and intersection.
    File         samples_file_1       = ExtractSamples1.samples_file
    File         samples_file_2       = ExtractSamples2.samples_file
    File         common_samples       = FindCommonSamples.common_samples_file
    Int          n_samples_file1      = ExtractSamples1.n_samples
    Int          n_samples_file2      = ExtractSamples2.n_samples
    Int          n_common             = FindCommonSamples.n_common

    # Per-sample query/ref extractions.
    # query schema : #chrom start end name   SVTYPE SVLEN
    # ref   schema : #chrom start end VID    svtype length AF samples
    Array[File]  file1_query_files    = ExtractFromFile1.query_file
    Array[File]  file1_ref_files      = ExtractFromFile1.ref_file
    Array[File]  file2_query_files    = ExtractFromFile2.query_file
    Array[File]  file2_ref_files      = ExtractFromFile2.ref_file

    # Per-sample annotated comparison BEDs (one per direction per sample).
    # Schema: #chr start end VID svtype length AF ovr1a ovr1b ovr2a ovr2b ovr3 GC sample_id
    Array[File]  annotated_beds_1v2   = CompareAndAnnotateSample.annotated_bed_1v2
    Array[File]  annotated_beds_2v1   = CompareAndAnnotateSample.annotated_bed_2v1

    # All samples + both directions merged into one BED.
    File         merged_bed           = ConcatResults.merged_bed

    # Per-direction summary stat tables (one row per sample).
    # Columns: sample_id  total_SVs  ovr_SVs  total_US_RM  ovr_US_RM  total_SD_SR  ovr_SD_SR
    File         stats_table_1v2      = ConcatStats.stats_table_1v2
    File         stats_table_2v1      = ConcatStats.stats_table_2v1
  }
}


# ---------------------------------------------------------------------------
# Task: ExtractSamples
#   Extracts the full sorted unique sample list from a VCF or BED file.
# ---------------------------------------------------------------------------
task ExtractSamples {
  input {
    File   variant_file
    String bed_sample_col
    String output_prefix
    String bcftools_docker
  }

  String samples_out = output_prefix + ".samples.txt"

  command <<<
    set -euo pipefail

    filename="~{variant_file}"
    base="${filename%.gz}"

    if [[ "$base" == *.vcf || "$base" == *.bcf ]]; then
      bcftools query -l "~{variant_file}" | sort -u > "~{samples_out}"

    elif [[ "$base" == *.bed ]]; then
      col_idx=$(zcat -f "~{variant_file}" \
                  | grep '^#' \
                  | grep -v '^##' \
                  | head -1 \
                  | awk -v col="~{bed_sample_col}" '
                      {
                        sub(/^#/, "", $1)
                        for (i = 1; i <= NF; i++) {
                          if ($i == col) { print i; exit }
                        }
                      }')

      if [[ -z "$col_idx" ]]; then
        echo "ERROR: column '~{bed_sample_col}' not found in BED header" >&2
        exit 1
      fi

      zcat -f "~{variant_file}" \
        | grep -v '^#' \
        | awk -v col="${col_idx}" '{ print $col }' \
        | tr ',' '\n' \
        | awk 'NF' \
        | sort -u > "~{samples_out}"
    else
      echo "ERROR: unrecognised extension for '~{variant_file}'." \
           "Supported: .vcf, .vcf.gz, .bcf, .bed, .bed.gz" >&2
      exit 1
    fi

    wc -l < "~{samples_out}" > n_samples.txt
    echo "Extracted $(cat n_samples.txt) unique sample(s) from ~{variant_file}"
  >>>

  output {
    File samples_file = samples_out
    Int  n_samples    = read_int("n_samples.txt")
  }

  runtime {
    cpu:    1
    memory: "4 GiB"
    disks:  "local-disk 50 HDD"
    docker: bcftools_docker
  }
}


# ---------------------------------------------------------------------------
# Task: FindCommonSamples
#   Intersects two sample lists and writes the common sample IDs.
# ---------------------------------------------------------------------------
task FindCommonSamples {
  input {
    File   samples_file_1
    File   samples_file_2
    String output_prefix
    String bcftools_docker
  }

  String common_out = output_prefix + ".common_samples.txt"

  command <<<
    set -euo pipefail

    comm -12 \
      <(sort "~{samples_file_1}") \
      <(sort "~{samples_file_2}") \
      > "~{common_out}"

    n1=$(wc -l < "~{samples_file_1}")
    n2=$(wc -l < "~{samples_file_2}")
    nc=$(wc -l < "~{common_out}")

    echo "Samples in file 1 : ${n1}"
    echo "Samples in file 2 : ${n2}"
    echo "Common samples     : ${nc}"

    echo "${nc}" > n_common.txt
  >>>

  output {
    File common_samples_file = common_out
    Int  n_common            = read_int("n_common.txt")
  }

  runtime {
    cpu:    1
    memory: "2 GiB"
    disks:  "local-disk 10 HDD"
    docker: bcftools_docker
  }
}


# ---------------------------------------------------------------------------
# Task: ExtractVariantsForSample
#   For one sample, extracts its non-ref variants from a single VCF or BED
#   file and writes two output files:
#
#   query.gz  columns: #chrom  start  end  name   SVTYPE  SVLEN
#   ref.gz    columns: #chrom  start  end  VID    svtype  length  AF  samples
#
#   Output filenames: {sample_id}.{file_label}.query.gz
#                     {sample_id}.{file_label}.ref.gz
#
#   VCF notes:
#     - start  = POS - 1  (convert to 0-based BED start)
#     - end    = INFO/END  (already equivalent to 0-based exclusive BED end)
#     - SVTYPE, SVLEN, AF pulled from the INFO field
#     - SVLEN stored as absolute value in both output files
#
#   BED notes (layout from example.bed):
#     #chrom  start  end  name  svtype  samples  SVTYPE  SVLEN  AF
#     - query SVTYPE  ← SVTYPE column (uppercase)
#     - ref   svtype  ← svtype column (lowercase)
#     - length        ← SVLEN column (absolute value)
#     - samples field in ref output contains only the current sample_id
# ---------------------------------------------------------------------------
task ExtractVariantsForSample {
  input {
    File   variant_file
    String sample_id
    String bed_sample_col
    String file_label
    String bcftools_docker
  }

  String query_out = sample_id + "." + file_label + ".query.gz"
  String ref_out   = sample_id + "." + file_label + ".ref.gz"

  command <<<
    set -euo pipefail

    filename="~{variant_file}"
    base="${filename%.gz}"

    if [[ "$base" == *.vcf || "$base" == *.bcf ]]; then
      # ------------------------------------------------------------------ VCF
      # Subset to this sample and keep only non-ref, non-missing genotypes.
      bcftools view \
        -s "~{sample_id}" \
        -i 'GT[0]!="RR" && GT[0]!="mis"' \
        -O z -o filtered.vcf.gz \
        "~{variant_file}"

      # query format: #chrom start end name SVTYPE SVLEN
      {
        printf "#chrom\tstart\tend\tname\tSVTYPE\tSVLEN\n"
        bcftools query \
          -f '%CHROM\t%POS\t%END\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN\n' \
          filtered.vcf.gz \
        | awk 'BEGIN{OFS="\t"} {
            $2 = $2 - 1                    # POS (1-based) → start (0-based)
            $6 = ($6 < 0 ? -$6 : $6)      # abs(SVLEN)
            print
          }'
      } | bgzip -c > "~{query_out}"

      # ref format: #chrom start end VID svtype length AF samples
      {
        printf "#chrom\tstart\tend\tVID\tsvtype\tlength\tAF\tsamples\n"
        bcftools query \
          -f '%CHROM\t%POS\t%END\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/AF\n' \
          filtered.vcf.gz \
        | awk -v samp="~{sample_id}" 'BEGIN{OFS="\t"} {
            $2 = $2 - 1
            len = ($6 < 0 ? -$6 : $6)
            print $1, $2, $3, $4, $5, len, $7, samp
          }'
      } | bgzip -c > "~{ref_out}"

    elif [[ "$base" == *.bed ]]; then
      # ------------------------------------------------------------------ BED
      # Single awk pass: parse header to locate column indices, filter rows
      # to those carrying this sample, write query and ref to separate files.
      zcat -f "~{variant_file}" | awk \
        -v sample="~{sample_id}" \
        -v samp_col="~{bed_sample_col}" \
        -v qout="query.tsv" \
        -v rout="ref.tsv" \
        'BEGIN { OFS="\t" }
        NR == 1 {
          sub(/^#/, "", $1)
          for (i = 1; i <= NF; i++) col[$i] = i
          chrom_c        = col["chrom"]
          start_c        = col["start"]
          end_c          = col["end"]
          name_c         = col["name"]
          svtype_upper_c = (col["SVTYPE"] ? col["SVTYPE"] : col["svtype"])
          svtype_lower_c = (col["svtype"] ? col["svtype"] : col["SVTYPE"])
          svlen_c        = (col["SVLEN"]  ? col["SVLEN"]  : col["svlen"])
          af_c           = col["AF"]
          samp_c         = col[samp_col]
          print "#chrom", "start", "end", "name",  "SVTYPE", "SVLEN"                  > qout
          print "#chrom", "start", "end", "VID",   "svtype", "length", "AF", "samples" > rout
          next
        }
        {
          n = split($samp_c, sa, ",")
          found = 0
          for (i = 1; i <= n; i++) if (sa[i] == sample) { found = 1; break }
          if (!found) next
          svlen = $svlen_c + 0
          len   = (svlen < 0 ? -svlen : svlen)
          print $chrom_c, $start_c, $end_c, $name_c, $svtype_upper_c, $svlen_c             >> qout
          print $chrom_c, $start_c, $end_c, $name_c, $svtype_lower_c, len, $af_c, sample   >> rout
        }'

      bgzip -c query.tsv > "~{query_out}"
      bgzip -c ref.tsv   > "~{ref_out}"

    else
      echo "ERROR: unrecognised extension for '~{variant_file}'." \
           "Supported: .vcf, .vcf.gz, .bcf, .bed, .bed.gz" >&2
      exit 1
    fi
  >>>

  output {
    File query_file = query_out
    File ref_file   = ref_out
  }

  runtime {
    cpu:    1
    memory: "4 GiB"
    disks:  "local-disk 50 HDD"
    docker: bcftools_docker
  }
}


# ---------------------------------------------------------------------------
# Task: CompareAndAnnotateSample
#   For one common sample, runs four steps:
#
#   1. compare_callsets_V2.sh in both directions:
#        {label1}vs{label2}.{sample}.bed  ← file1.query vs file2.ref
#        {label2}vs{label1}.{sample}.bed  ← file2.query vs file1.ref
#      compare output schema:
#        #chr start end VID svtype length AF ovr1a ovr1b ovr2a ovr2b ovr3
#
#   2. annotate_genomic_context.sh on each comparison BED:
#        output is a 2-column file: VID  GC
#
#   3. Join GC into the comparison BEDs by matching on VID (col 4).
#
#   4. Append a sample_id column and bgzip.
#
#   Both scripts and the companion R script are copied to the working
#   directory so relative-path lookups inside the scripts resolve correctly.
# ---------------------------------------------------------------------------
task CompareAndAnnotateSample {
  input {
    File   file1_query                    # {sample}.{label1}.query.gz
    File   file1_ref                      # {sample}.{label1}.ref.gz
    File   file2_query                    # {sample}.{label2}.query.gz
    File   file2_ref                      # {sample}.{label2}.ref.gz
    String sample_id
    String file_label_1
    String file_label_2
    File   compare_callsets_script        # compare_callsets_V2.sh
    File   companion_r_script            # R script required by compare_callsets_V2.sh
    File   annotate_genomic_context_script
    File   hg38_repeat_mask_bed
    File   hg38_seg_dup_bed
    File   hg38_simp_rep_bed
    String scripts_docker
  }

  String label_1v2 = file_label_1 + "vs" + file_label_2
  String label_2v1 = file_label_2 + "vs" + file_label_1
  String bed_1v2   = label_1v2 + "." + sample_id + ".bed"
  String bed_2v1   = label_2v1 + "." + sample_id + ".bed"
  String gc_1v2    = label_1v2 + "." + sample_id + ".SVID_GC"
  String gc_2v1    = label_2v1 + "." + sample_id + ".SVID_GC"
  String out_1v2   = label_1v2 + "." + sample_id + ".annotated.bed.gz"
  String out_2v1   = label_2v1 + "." + sample_id + ".annotated.bed.gz"

  command <<<
    set -euo pipefail

    # Co-localise scripts so any relative-path sourcing inside them works.
    cp "~{compare_callsets_script}"         compare_callsets_V2.sh
    cp "~{companion_r_script}"              "$(basename ~{companion_r_script})"
    cp "~{annotate_genomic_context_script}" annotate_genomic_context.sh
    chmod +x compare_callsets_V2.sh annotate_genomic_context.sh

    # ---- Step 1: Compare callsets (both directions) ------------------------
    bash compare_callsets_V2.sh \
      -O "~{bed_1v2}" \
      -p "~{label_1v2}" \
      ~{file1_query} ~{file2_ref}

    bash compare_callsets_V2.sh \
      -O "~{bed_2v1}" \
      -p "~{label_2v1}" \
      ~{file2_query} ~{file1_ref}

    # ---- Step 2: Annotate genomic context (both directions) ----------------
    bash annotate_genomic_context.sh \
      "~{bed_1v2}" "~{gc_1v2}" \
      --rm ~{hg38_repeat_mask_bed} \
      --sd ~{hg38_seg_dup_bed} \
      --sr ~{hg38_simp_rep_bed}

    bash annotate_genomic_context.sh \
      "~{bed_2v1}" "~{gc_2v1}" \
      --rm ~{hg38_repeat_mask_bed} \
      --sd ~{hg38_seg_dup_bed} \
      --sr ~{hg38_simp_rep_bed}

    # ---- Step 3: Join GC column into BED by matching VID (col 4) ----------
    # GC file format: VID  GC  (optional single-line header where col1=="VID")
    join_gc() {
      local bed="$1"
      local gc_file="$2"
      local tmp="${bed}.tmp"
      awk 'BEGIN { OFS="\t" }
           NR == FNR {
             if (NR == 1 && $1 == "VID") next   # skip optional header
             gc[$1] = $2
             next
           }
           FNR == 1 { print $0 "\tGC"; next }
           { print $0 "\t" (($4 in gc) ? gc[$4] : "NA") }' \
        "$gc_file" "$bed" > "$tmp"
      mv "$tmp" "$bed"
    }

    join_gc "~{bed_1v2}" "~{gc_1v2}"
    join_gc "~{bed_2v1}" "~{gc_2v1}"

    # ---- Step 4: Append sample_id column and bgzip -------------------------
    add_sample() {
      local bed="$1"
      local outgz="$2"
      awk -v samp="~{sample_id}" 'BEGIN { OFS="\t" }
        FNR == 1 { print $0 "\tsample_id"; next }
        { print $0 "\t" samp }' \
        "$bed" | bgzip -c > "$outgz"
    }

    add_sample "~{bed_1v2}" "~{out_1v2}"
    add_sample "~{bed_2v1}" "~{out_2v1}"
  >>>

  output {
    # Schema: #chr start end VID svtype length AF ovr1a ovr1b ovr2a ovr2b ovr3 GC sample_id
    File annotated_bed_1v2 = out_1v2
    File annotated_bed_2v1 = out_2v1
  }

  runtime {
    cpu:    2
    memory: "8 GiB"
    disks:  "local-disk 100 HDD"
    docker: scripts_docker
  }
}


# ---------------------------------------------------------------------------
# Task: ConcatResults
#   Concatenates all per-sample annotated BED.gz files (both comparison
#   directions, all samples) into a single merged BED.gz.
#   Header is taken from the first file; subsequent files have their header
#   line skipped.
# ---------------------------------------------------------------------------
task ConcatResults {
  input {
    Array[File] annotated_beds   # output of CompareAndAnnotateSample, all samples + directions
    String      output_prefix
    String      bcftools_docker
  }

  String merged_out = output_prefix + ".merged.bed.gz"

  command <<<
    set -euo pipefail

    files=(~{sep=" " annotated_beds})
    n=${#files[@]}

    if [[ $n -eq 0 ]]; then
      echo "ERROR: no annotated BEDs to concatenate" >&2
      exit 1
    fi

    # Header from the first file; data rows (skip header) from every file.
    {
      zcat "${files[0]}" | head -1
      for f in "${files[@]}"; do
        zcat "$f" | tail -n +2
      done
    } | bgzip -c > "~{merged_out}"

    total=$(zcat "~{merged_out}" | tail -n +2 | wc -l)
    echo "Merged ${n} files → ${total} data rows in ~{merged_out}"
  >>>

  output {
    File merged_bed = merged_out
  }

  runtime {
    cpu:    2
    memory: "8 GiB"
    disks:  "local-disk 100 HDD"
    docker: bcftools_docker
  }
}


# ---------------------------------------------------------------------------
# Task: ComputeSampleStats
#   For one common sample, computes six counts from each of the two annotated
#   comparison BEDs and writes a single header-free TSV row per direction.
#
#   Columns (written without header; header added by ConcatStats):
#     sample_id
#     total_SVs         — all data rows
#     ovr_SVs           — rows where ovr1a != "NO_OVR"
#     total_US_RM       — rows where GC is "US" or "RM"
#     ovr_US_RM         — rows where ovr1a != "NO_OVR" AND GC is "US" or "RM"
#     total_SD_SR       — rows where GC is "SD" or "SR"
#     ovr_SD_SR         — rows where ovr1a != "NO_OVR" AND GC is "SD" or "SR"
# ---------------------------------------------------------------------------
task ComputeSampleStats {
  input {
    File   annotated_bed_1v2   # {label1}vs{label2}.{sample}.annotated.bed.gz
    File   annotated_bed_2v1   # {label2}vs{label1}.{sample}.annotated.bed.gz
    String sample_id
    String bcftools_docker
  }

  String stats_out_1v2 = sample_id + ".stats_1v2.tsv"
  String stats_out_2v1 = sample_id + ".stats_2v1.tsv"

  command <<<
    set -euo pipefail

    # Shared awk program: dynamically locates ovr1a and GC columns from the
    # header, then accumulates the six counts and prints one data row.
    count_stats() {
      local bed_gz="$1"
      local sample="$2"
      local outfile="$3"

      zcat "$bed_gz" | awk \
        -v samp="$sample" \
        'BEGIN { OFS="\t" }
        NR == 1 {
          for (i = 1; i <= NF; i++) col[$i] = i
          ovr_c = col["ovr1a"]
          gc_c  = col["GC"]
          if (!ovr_c) { print "ERROR: ovr1a column not found" > "/dev/stderr"; exit 1 }
          if (!gc_c)  { print "ERROR: GC column not found"    > "/dev/stderr"; exit 1 }
          next
        }
        {
          total++
          ovr    = ($ovr_c != "NO_OVR")
          gc     = $gc_c
          us_rm  = (gc == "US" || gc == "RM")
          sd_sr  = (gc == "SD" || gc == "SR")

          if (ovr)           ovr_total++
          if (us_rm)         total_us_rm++
          if (ovr && us_rm)  ovr_us_rm++
          if (sd_sr)         total_sd_sr++
          if (ovr && sd_sr)  ovr_sd_sr++
        }
        END {
          print samp,
                total       + 0,
                ovr_total   + 0,
                total_us_rm + 0,
                ovr_us_rm   + 0,
                total_sd_sr + 0,
                ovr_sd_sr   + 0
        }' > "$outfile"
    }

    count_stats "~{annotated_bed_1v2}" "~{sample_id}" "~{stats_out_1v2}"
    count_stats "~{annotated_bed_2v1}" "~{sample_id}" "~{stats_out_2v1}"
  >>>

  output {
    File stats_row_1v2 = stats_out_1v2
    File stats_row_2v1 = stats_out_2v1
  }

  runtime {
    cpu:    1
    memory: "2 GiB"
    disks:  "local-disk 20 HDD"
    docker: bcftools_docker
  }
}


# ---------------------------------------------------------------------------
# Task: ConcatStats
#   Combines per-sample stat rows (from ComputeSampleStats) into two
#   tab-separated summary tables, one per comparison direction.
#
#   Output header:
#     sample_id  total_SVs  ovr_SVs  total_US_RM  ovr_US_RM  total_SD_SR  ovr_SD_SR
# ---------------------------------------------------------------------------
task ConcatStats {
  input {
    Array[File] stat_rows_1v2   # one row per sample for direction 1→2
    Array[File] stat_rows_2v1   # one row per sample for direction 2→1
    String      file_label_1
    String      file_label_2
    String      output_prefix
    String      bcftools_docker
  }

  String table_out_1v2 = output_prefix + "." + file_label_1 + "vs" + file_label_2 + ".stats.tsv"
  String table_out_2v1 = output_prefix + "." + file_label_2 + "vs" + file_label_1 + ".stats.tsv"

  command <<<
    set -euo pipefail

    HEADER="sample_id\ttotal_SVs\tovr_SVs\ttotal_US_RM\tovr_US_RM\ttotal_SD_SR\tovr_SD_SR"

    build_table() {
      local outfile="$1"
      shift
      local rows=("$@")
      printf "${HEADER}\n" > "$outfile"
      for f in "${rows[@]}"; do
        cat "$f" >> "$outfile"
      done
    }

    build_table "~{table_out_1v2}" ~{sep=" " stat_rows_1v2}
    build_table "~{table_out_2v1}" ~{sep=" " stat_rows_2v1}

    echo "Stats table 1v2: $(wc -l < ~{table_out_1v2}) lines (including header)"
    echo "Stats table 2v1: $(wc -l < ~{table_out_2v1}) lines (including header)"
  >>>

  output {
    File stats_table_1v2 = table_out_1v2
    File stats_table_2v1 = table_out_2v1
  }

  runtime {
    cpu:    1
    memory: "2 GiB"
    disks:  "local-disk 10 HDD"
    docker: bcftools_docker
  }
}
