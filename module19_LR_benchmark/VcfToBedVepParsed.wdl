version 1.0

## VcfToBedVepParsed
## Converts a list of annotated VCF(.gz) files to BED format, splitting the VEP
## INFO field into three columns (vep_Consequence, vep_IMPACT, vep_SYMBOL).
## Scatters over all input VCFs in parallel, then concatenates into one BED.
##
## Output columns (40):
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

workflow VcfToBedVepParsed {

    input {
        Array[File]  input_vcfs          # list of annotated VCF(.gz) files
        File         script              # vcf_to_bed_vep_parsed.py
        String       output_basename     # prefix for the merged output BED
        String       docker = "python:3.11-slim"
        Int          mem_gb        = 8
        Int          cpu           = 2
        Int          disk_gb       = 100
        Int          preemptible   = 1
    }

    # ── scatter: one task per VCF ────────────────────────────────────────────
    scatter (vcf in input_vcfs) {
        call ConvertVcfToBed {
            input:
                vcf          = vcf,
                script       = script,
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

    output {
        File merged_bed = ConcatBeds.merged_bed
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
        python3 ~{script} ~{vcf} ~{out_bed}
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
