version 1.0

## VcfToBedVepParsed
## Converts a list of annotated VCF(.gz) files to BED format, splitting the VEP
## INFO field into three columns (vep_Consequence, vep_IMPACT, vep_SYMBOL).
## Scatters over all input VCFs in parallel, then concatenates into one BED,
## sorts it, and bgzips + tabix-indexes the final merged BED.
##
## Output columns (40, plus any extra_info_fields appended at the end):
##   #CHROM START END ID REF ALT QUAL FILTER
##   allele_type allele_length SOURCE REGION TRID dbGaP_ID
##   gnomAD_V4_match_type gnomAD_V4_match_ID gnomAD_V4_match_source AF AC AN
##   PREDICTED_BREAKEND_EXONIC PREDICTED_COPY_GAIN PREDICTED_DUP_PARTIAL
##   PREDICTED_INTERGENIC PREDICTED_INTRAGENIC_EXON_DUP PREDICTED_INTRONIC
##   PREDICTED_INV_SPAN PREDICTED_LOF PREDICTED_MSV_EXON_OVERLAP
##   PREDICTED_NEAREST_TSS PREDICTED_NONCODING_BREAKPOINT PREDICTED_NONCODING_SPAN
##   PREDICTED_PARTIAL_DISPERSED_DUP PREDICTED_PARTIAL_EXON_DUP PREDICTED_PROMOTER
##   PREDICTED_TSS_DUP PREDICTED_UTR
##   vep_Consequence vep_IMPACT vep_SYMBOL
##   [extra_info_fields, e.g. nhomref nhet nhomalt]
##
## extra_info_fields (optional): comma-separated list of additional site-level
## INFO field names to pull straight through into extra BED columns, e.g.
##   extra_info_fields = "nhomref,nhet,nhomalt"
## for a VCF annotated with:
##   ##INFO=<ID=nhomref,Number=1,Type=Integer,Description="Number of samples with homozygous reference genotypes (biallelic sites only).">
##   ##INFO=<ID=nhet,Number=1,Type=Integer,Description="Number of samples with heterozygous genotypes (biallelic sites only).">
##   ##INFO=<ID=nhomalt,Number=1,Type=Integer,Description="Number of samples with homozygous alternate genotypes (biallelic sites only).">

workflow VcfToBedVepParsed {

    input {
        Array[File]  input_vcfs          # list of annotated VCF(.gz) files
        File         script              # vcf_to_bed_vep_parsed.py
        String       output_basename     # prefix for the merged output BED
        String       extra_info_fields = ""  # comma-separated extra INFO field names, e.g. "nhomref,nhet,nhomalt"
        String       docker = "python:3.11-slim"
        String       htslib_docker = "staphb/htslib:1.19"  # needs bgzip + tabix
        Int          mem_gb        = 8
        Int          cpu           = 2
        Int          disk_gb       = 100
        Int          preemptible   = 1
    }

    # ── scatter: one task per VCF ────────────────────────────────────────────
    scatter (vcf in input_vcfs) {
        call ConvertVcfToBed {
            input:
                vcf               = vcf,
                script            = script,
                extra_info_fields = extra_info_fields,
                docker            = docker,
                mem_gb            = mem_gb,
                cpu               = cpu,
                disk_gb           = disk_gb,
                preemptible       = preemptible
        }

        call ExtractIdFilter {
            input:
                bed          = ConvertVcfToBed.bed,
                docker       = docker,
                mem_gb       = mem_gb,
                cpu          = cpu,
                disk_gb      = disk_gb,
                preemptible  = preemptible
        }
    }

    # ── gather: merge all BED shards ─────────────────────────────────────────
    call ConcatBeds {
        input:
            bed_files       = ConvertVcfToBed.bed,
            output_basename = output_basename,
            docker          = docker,
            mem_gb          = mem_gb,
            cpu             = 2,
            disk_gb         = disk_gb,
            preemptible     = preemptible
    }

    # ── gather: merge all ID/FILTER shards ───────────────────────────────────
    call ConcatIdFilter {
        input:
            id_filter_files = ExtractIdFilter.id_filter,
            output_basename = output_basename,
            docker          = docker,
            mem_gb          = mem_gb,
            cpu             = 2,
            disk_gb         = disk_gb,
            preemptible     = preemptible
    }

    # ── sort + bgzip + tabix-index the merged BED ────────────────────────────
    call BgzipTabixBed {
        input:
            bed             = ConcatBeds.merged_bed,
            output_basename = output_basename,
            docker          = htslib_docker,
            mem_gb          = mem_gb,
            cpu             = cpu,
            disk_gb         = disk_gb,
            preemptible     = preemptible
    }

    output {
        File merged_bed_gz     = BgzipTabixBed.bed_gz
        File merged_bed_gz_tbi = BgzipTabixBed.bed_gz_tbi
        File merged_id_filter  = ConcatIdFilter.merged_id_filter
    }

    meta {
        author: "gnomAD LR analysis"
        description: "Convert annotated VCFs to BED with parsed VEP columns and merge."
    }
}

# ── Task: convert one VCF to BED ─────────────────────────────────────────────
task ConvertVcfToBed {

    input {
        File    vcf
        File    script
        String  extra_info_fields
        String  docker
        Int     mem_gb
        Int     cpu
        Int     disk_gb
        Int     preemptible
    }

    # derive output name from the VCF filename
    String vcf_basename = basename(vcf, ".vcf.gz")
    String out_bed      = vcf_basename + ".vep_parsed.bed"

    command <<<
        set -euo pipefail
        python3 ~{script} ~{vcf} ~{out_bed} \
            ~{if extra_info_fields != "" then "--extra-info-fields " + extra_info_fields else ""}
    >>>

    output {
        File bed = out_bed
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         cpu
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

# ── Task: extract ID and FILTER columns from one BED shard ───────────────────
task ExtractIdFilter {

    input {
        File    bed
        String  docker
        Int     mem_gb
        Int     cpu
        Int     disk_gb
        Int     preemptible
    }

    # BED columns are #CHROM START END ID REF ALT QUAL FILTER ... -- ID is
    # column 4, FILTER is column 8.
    String bed_basename = basename(bed, ".bed")
    String out_name      = bed_basename + ".id_filter.tsv"

    command <<<
        set -euo pipefail
        cut -f4,8 ~{bed} > ~{out_name}
    >>>

    output {
        File id_filter = out_name
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         cpu
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

# ── Task: concatenate shard BEDs into one file ────────────────────────────────
task ConcatBeds {

    input {
        Array[File] bed_files
        String      output_basename
        String      docker
        Int         mem_gb
        Int         cpu
        Int         disk_gb
        Int         preemptible
    }

    String merged = output_basename + ".vep_parsed.bed"

    command <<<
        set -euo pipefail

        BED_FILES=(~{sep=" " bed_files})
        OUT="~{merged}"

        # Write header from the first shard
        head -1 "${BED_FILES[0]}" > "$OUT"

        # Append data rows (skip header line 1) from every shard
        for f in "${BED_FILES[@]}"; do
            tail -n +2 "$f" >> "$OUT"
        done

        echo "Merged $(( $(wc -l < "$OUT") - 1 )) variants → $OUT"
    >>>

    output {
        File merged_bed = merged
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         cpu
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

# ── Task: concatenate shard ID/FILTER files into one file ────────────────────
task ConcatIdFilter {

    input {
        Array[File] id_filter_files
        String      output_basename
        String      docker
        Int         mem_gb
        Int         cpu
        Int         disk_gb
        Int         preemptible
    }

    String merged = output_basename + ".id_filter.tsv"

    command <<<
        set -euo pipefail

        ID_FILTER_FILES=(~{sep=" " id_filter_files})
        OUT="~{merged}"

        # Write header from the first shard
        head -1 "${ID_FILTER_FILES[0]}" > "$OUT"

        # Append data rows (skip header line 1) from every shard
        for f in "${ID_FILTER_FILES[@]}"; do
            tail -n +2 "$f" >> "$OUT"
        done

        echo "Merged $(( $(wc -l < "$OUT") - 1 )) variants → $OUT"
    >>>

    output {
        File merged_id_filter = merged
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         cpu
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

# ── Task: sort, bgzip, and tabix-index the merged BED ────────────────────────
task BgzipTabixBed {

    input {
        File    bed
        String  output_basename
        String  docker
        Int     mem_gb
        Int     cpu
        Int     disk_gb
        Int     preemptible
    }

    String out_bed_gz = output_basename + ".vep_parsed.bed.gz"

    command <<<
        set -euo pipefail

        # Keep the "#CHROM..." header as a leading comment line (tabix's
        # default meta-char is '#', so it's skipped on indexing), sort the
        # data rows by chrom/start/end, then bgzip + tabix as BED.
        head -1 ~{bed} > sorted.bed
        tail -n +2 ~{bed} | sort -k1,1 -k2,2n -k3,3n >> sorted.bed

        bgzip -c sorted.bed > ~{out_bed_gz}
        tabix -p bed ~{out_bed_gz}
    >>>

    output {
        File bed_gz     = out_bed_gz
        File bed_gz_tbi = out_bed_gz + ".tbi"
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         cpu
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}
