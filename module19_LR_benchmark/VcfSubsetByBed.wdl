version 1.0

## vcf_subset_by_bed.wdl
##
## For each input VCF (on GCS):
##   1. Expand each BED region by `padding_bp` on both sides
##   2. Stream the VCF through bcftools view --targets-file to extract variants
##      falling within the padded regions
##   3. Concatenate all per-VCF subsets into a single VCF
##   4. bgzip + tabix index the final output
##
## Notes:
##   - Uses bcftools --targets (not --regions) so no local index is needed;
##     works with gsutil-streamed input.
##   - The .tbi files are declared as inputs so Cromwell localises them
##     alongside the VCFs, satisfying any co-location requirements.
##   - Requires bcftools >= 1.10 (for --targets-file BED support).

workflow VcfSubsetByBed {

    input {
        Array[File] vcf_gz_files          # gs:// paths to VCF.gz files
        Array[File] vcf_tbi_files         # gs:// paths to matching .tbi files (same order)
        File        bed_file              # BED file defining regions of interest
        Int         padding_bp    = 100   # bp to pad each BED region on both sides
        String      output_prefix = "merged_subset"

        # Runtime
        String      bcftools_docker = "quay.io/biocontainers/bcftools:1.19--h8b25389_1"
        String      htslib_docker   = "quay.io/biocontainers/htslib:1.19--h81da01d_0"
        Int         subset_mem_gb   = 8
        Int         subset_cpu      = 2
        Int         concat_mem_gb   = 16
        Int         concat_cpu      = 4
    }

    # ──────────────────────────────────────────────────────────────────────────
    # Step 1 — Pad the BED file once (shared across all VCF tasks)
    # ──────────────────────────────────────────────────────────────────────────

    call PadBed {
        input:
            bed_file   = bed_file,
            padding_bp = padding_bp,
            docker     = bcftools_docker,
            mem_gb     = 2
    }

    # ──────────────────────────────────────────────────────────────────────────
    # Step 2 — Extract variants per VCF (scattered)
    # ──────────────────────────────────────────────────────────────────────────

    Array[Pair[File, File]] vcf_pairs = zip(vcf_gz_files, vcf_tbi_files)

    scatter (pair in vcf_pairs) {
        call SubsetVcf {
            input:
                vcf_gz      = pair.left,
                vcf_tbi     = pair.right,
                padded_bed  = PadBed.padded_bed,
                docker      = bcftools_docker,
                mem_gb      = subset_mem_gb,
                cpu         = subset_cpu
        }
    }

    # ──────────────────────────────────────────────────────────────────────────
    # Step 3 — Concatenate all subset VCFs
    # ──────────────────────────────────────────────────────────────────────────

    call ConcatVcfs {
        input:
            subset_vcfs   = SubsetVcf.subset_vcf,
            output_prefix = output_prefix,
            docker        = bcftools_docker,
            mem_gb        = concat_mem_gb,
            cpu           = concat_cpu
    }

    # ──────────────────────────────────────────────────────────────────────────
    # Step 4 — bgzip + tabix
    # ──────────────────────────────────────────────────────────────────────────

    call BgzipAndTabix {
        input:
            vcf_file      = ConcatVcfs.concat_vcf,
            output_prefix = output_prefix,
            docker        = htslib_docker,
            mem_gb        = 4
    }

    # ──────────────────────────────────────────────────────────────────────────
    # Outputs
    # ──────────────────────────────────────────────────────────────────────────

    output {
        File final_vcf_gz  = BgzipAndTabix.vcf_gz
        File final_vcf_tbi = BgzipAndTabix.vcf_tbi
        File padded_bed    = PadBed.padded_bed
    }
}


# ════════════════════════════════════════════════════════════════════════════
# TASK DEFINITIONS
# ════════════════════════════════════════════════════════════════════════════

task PadBed {
    ## Expand every BED interval by padding_bp on both sides.
    ## Negative coordinates are clamped to 0.
    ## Column 4+ (name, score, strand …) are preserved if present.

    input {
        File   bed_file
        Int    padding_bp
        String docker
        Int    mem_gb
    }

    command <<<
        set -euo pipefail

        python3 - <<'PYEOF'
import sys

pad = int("~{padding_bp}")
out_lines = []

with open("~{bed_file}") as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
            out_lines.append(line)
            continue
        cols = line.split("\t")
        chrom = cols[0]
        start = max(0, int(cols[1]) - pad)
        end   = int(cols[2]) + pad
        rest  = cols[3:] if len(cols) > 3 else []
        out_lines.append("\t".join([chrom, str(start), str(end)] + rest))

with open("padded_regions.bed", "w") as fh:
    fh.write("\n".join(out_lines) + "\n")

print(f"Padded {sum(1 for l in out_lines if l and not l.startswith('#'))} regions "
      f"by ±{pad} bp", file=sys.stderr)
PYEOF
    >>>

    output {
        File padded_bed = "padded_regions.bed"
    }

    runtime {
        docker: docker
        memory: "~{mem_gb} GB"
        cpu:    1
        disks:  "local-disk 10 HDD"
    }
}


task SubsetVcf {
    ## Stream a GCS VCF through bcftools view --targets-file.
    ## --targets works on stdin so no local index is required.
    ## The .tbi is declared as input so Cromwell stages it; even if not
    ## directly used by the stream command, some backends require co-location.

    input {
        File   vcf_gz
        File   vcf_tbi
        File   padded_bed
        String docker
        Int    mem_gb
        Int    cpu
    }

    String sample_name = basename(vcf_gz, ".vcf.gz")

    command <<<
        set -euo pipefail

        echo "Processing ~{sample_name}" >&2

        # Try direct bcftools access first (works if htslib has GCS support).
        # Fall back to gsutil cat streaming if it fails.
        if bcftools view --targets-file "~{padded_bed}" \
               "~{vcf_gz}" \
               -O v \
               -o "~{sample_name}.subset.vcf" 2>/dev/null; then
            echo "Used direct htslib GCS access" >&2
        else
            echo "Falling back to gsutil cat stream" >&2
            gsutil cat "~{vcf_gz}" \
                | bcftools view \
                    --targets-file "~{padded_bed}" \
                    -O v \
                    -o "~{sample_name}.subset.vcf"
        fi

        # Report variant count (excluding header)
        n=$(grep -vc "^#" "~{sample_name}.subset.vcf" || true)
        echo "Extracted ${n} variants from ~{sample_name}" >&2
    >>>

    output {
        File subset_vcf = "~{sample_name}.subset.vcf"
    }

    runtime {
        docker: docker
        memory: "~{mem_gb} GB"
        cpu:    cpu
        disks:  "local-disk 100 HDD"
    }
}


task ConcatVcfs {
    ## Concatenate subset VCFs from all samples.
    ## Uses bcftools concat --allow-overlaps to handle any duplicate positions
    ## that may appear across samples.

    input {
        Array[File] subset_vcfs
        String      output_prefix
        String      docker
        Int         mem_gb
        Int         cpu
    }

    command <<<
        set -euo pipefail

        # Write file list for bcftools concat
        VCF_LIST="vcf_list.txt"
        printf '%s\n' ~{sep=' ' subset_vcfs} > "${VCF_LIST}"

        echo "Concatenating $(wc -l < ${VCF_LIST}) VCF files..." >&2

        bcftools concat \
            --file-list "${VCF_LIST}" \
            --allow-overlaps \
            --remove-duplicates \
            -O v \
            -o "~{output_prefix}.concat.vcf"

        n=$(grep -vc "^#" "~{output_prefix}.concat.vcf" || true)
        echo "Total variants after concat: ${n}" >&2
    >>>

    output {
        File concat_vcf = "~{output_prefix}.concat.vcf"
    }

    runtime {
        docker: docker
        memory: "~{mem_gb} GB"
        cpu:    cpu
        disks:  "local-disk 200 HDD"
    }
}


task BgzipAndTabix {
    ## bgzip-compress and tabix-index the final concatenated VCF.

    input {
        File   vcf_file
        String output_prefix
        String docker
        Int    mem_gb
    }

    command <<<
        set -euo pipefail

        # Sort by coordinate before bgzip (concat may not be sorted if
        # input VCFs cover different chromosomes in different orders)
        bcftools sort \
            "~{vcf_file}" \
            -O z \
            -o "~{output_prefix}.vcf.gz"

        tabix -p vcf "~{output_prefix}.vcf.gz"

        echo "Final output: ~{output_prefix}.vcf.gz" >&2
        echo "Variants: $(bcftools view -H ~{output_prefix}.vcf.gz | wc -l)" >&2
    >>>

    output {
        File vcf_gz  = "~{output_prefix}.vcf.gz"
        File vcf_tbi = "~{output_prefix}.vcf.gz.tbi"
    }

    runtime {
        docker: docker
        memory: "~{mem_gb} GB"
        cpu:    2
        disks:  "local-disk 100 HDD"
    }
}
