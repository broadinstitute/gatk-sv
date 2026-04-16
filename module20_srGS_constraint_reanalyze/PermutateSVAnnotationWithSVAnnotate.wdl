version 1.0

import "Structs.wdl"

# ============================================================
# PermutateSVAnnotation.wdl
#
# Full pipeline to permute gene annotations and intersect SVs
# against each genomic context category.
#
# Inputs:
#   gtf_file          - Input GTF (can be .gtf or .gtf.gz)
#   permu_number      - Permutation seed (0 = original location)
#   sv_bed            - SV sites BED (gnomAD_SV_v3.sites.PASS.bed.gz)
#   tel_cen_bed       - Telomere/centromere exclusion BED
#   gtf_label         - Short label for the GTF (e.g. "r3.gencode.v39.ensembl.105")
#   sv_label          - Short label for the SV set (e.g. "gnomAD_SV_v3")
#   output_file_name  - Name for the final output file (step 8)
#
# Script inputs (localize from GCS or local paths):
#   permute_gtf_script, split_gtf_script,
#   categorize_r_script, reorganize_r_script,
#   integrate_r_script, calcu_r_script
# ============================================================

workflow PermutateSVAnnotationWithSVAnnotate {

    input {
        File   gtf_file
        File   vcf
        File   vcf_idx
        File   tel_cen_bed
        File   blacklist_bed
        File   gene_info

        Int    permu_number
        Array[String] contigs

        # Python scripts
        File   permute_gtf_script       # permute_gtf.py
        File   split_gtf_script         # split_gtf_annotations.py

        # R scripts
        File   categorize_r_script      # categorize_intact_vs_partial_exon_overlap.R
        File   reorganize_r_script      # reorganize_SVID_vs_gene.R
        File   integrate_r_script       # integrate_SVID_vs_genes.across_different_overlaps.R
        File   calcu_r_script           # calcu.gene.data.reanno.permu.R

        String gtf_label                # e.g. "r3.gencode.v39.ensembl.105"
        String sv_label                 # e.g. "gnomAD_SV_v3"

        # Docker images (override as needed)
        String python_docker   = "python:3.11-slim"
        String bedtools_docker = "quay.io/biocontainers/bedtools:2.31.1--h63f29b4_0"
        String r_docker        = "r-base:4.3.0"
        String gatk_docker =  "us.gcr.io/broad-dsde-methods/gatk-sv/gatk:2025-05-20-4.6.2.0-4-g1facd911e-NIGHTLY-SNAPSHOT"
        String utils_docker = "us-central1-docker.pkg.dev/talkowski-training/kj-development/utils:kj_V13"

        # Runtime attribute overrides per task
        RuntimeAttr? runtime_attr_permute_gtf
        RuntimeAttr? runtime_attr_split_gtf
        RuntimeAttr? runtime_attr_bedtools_intersect
        RuntimeAttr? runtime_attr_extract_overlaps
        RuntimeAttr? runtime_attr_categorize_exon
        RuntimeAttr? runtime_attr_reorganize_svid
        RuntimeAttr? runtime_attr_integrate_overlaps
        RuntimeAttr? runtime_attr_calcu_gene_data
    }

    String seed_suffix   = "permuted_seed" + permu_number
    String sv_gtf_prefix = sv_label + ".vs." + gtf_label

    # ── Task 1: Permute GTF ────────────────────────────────────────────────
    call Task1_PermuteGTF {
        input:
            gtf_file              = gtf_file,
            permu_number          = permu_number,
            tel_cen_bed           = tel_cen_bed,
            blacklist_bed         = blacklist_bed,
            permute_script        = permute_gtf_script,
            gtf_label             = gtf_label,
            seed_suffix           = seed_suffix,
            docker                = python_docker,
            runtime_attr_override = runtime_attr_permute_gtf
    }


    # ── Task 2: Annotate functional consequences of SVs with permutated gtf ─────────────
    
    call AnnotateFunctionalConsequences {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx, 
            prefix = "permu_~{permu_number}",
            coding_gtf = Task1_PermuteGTF.permuted_gtf,
            docker = gatk_docker
    }

}

task AnnotateFunctionalConsequences {
    input {
        File vcf
        File vcf_idx
        File coding_gtf
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Int java_mem_mb = 1000 * ceil(0.8 * select_first([runtime_attr.mem_gb, default_attr.mem_gb]))

    command <<<
        set -euo pipefail

        gatk --java-options "-Xmx~{java_mem_mb}m" SVAnnotate \
            -V ~{vcf} \
            --protein-coding-gtf ~{coding_gtf} \
            -O ~{prefix}.vcf.gz
    >>>

    output {
        File anno_vcf = "~{prefix}.vcf.gz"
        File anno_vcf_idx = "~{prefix}.vcf.gz.tbi"
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


# ======================================================================
# TASK 1 — Permute GTF
# ======================================================================
task Task1_PermuteGTF {
    input {
        File   gtf_file
        Int    permu_number
        File   tel_cen_bed
        File   blacklist_bed
        File   permute_script
        String gtf_label
        String seed_suffix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String out_gtf          = gtf_label + "." + seed_suffix + ".gtf.gz"
    String out_too_large     = gtf_label + "." + seed_suffix + ".too_large_genes.txt"

    command <<<
        python3 ~{permute_script} \
            ~{permu_number} \
            ~{gtf_file} \
            ~{tel_cen_bed} \
            ~{blacklist_bed} \
            ~{out_gtf} \
            ~{out_too_large}
    >>>

    output {
        File permuted_gtf         = out_gtf
        File too_large_genes_list = out_too_large
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         4,
        mem_gb:            8,
        disk_gb:           50,
        boot_disk_gb:      10,
        preemptible_tries: 1,
        max_retries:       1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu:            select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:         select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks:          "local-disk " + select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb,     default_attr.boot_disk_gb])
        docker:         docker
        preemptible:    select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:     select_first([runtime_attr.max_retries,       default_attr.max_retries])
    }
}

# ======================================================================
# TASK 2 — Split permuted GTF into 7 sorted, gzipped annotation BEDs
# ======================================================================
task Task2_SplitGTF {
    input {
        File   permuted_gtf
        File   split_script
        String gtf_label
        String seed_suffix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String prefix = gtf_label + "." + seed_suffix

    command <<<
        python3 ~{split_script} ~{permuted_gtf} ~{prefix}

        # Gzip all output BED files
        for f in ~{prefix}.*.bed; do
            gzip -f "$f"
        done
    >>>

    output {
        File transcript_bed        = prefix + ".transcript.bed.gz"
        File cds_bed               = prefix + ".CDS.bed.gz"
        File intron_bed            = prefix + ".intron.bed.gz"
        File utr3_bed              = prefix + ".utr_3.bed.gz"
        File utr5_bed              = prefix + ".utr_5.bed.gz"
        File promoter_bed          = prefix + ".promoter.bed.gz"
        File coding_transcript_bed = prefix + ".coding_transcript.bed.gz"
        Array[File] all_beds = [
            transcript_bed, cds_bed, intron_bed,
            utr3_bed, utr5_bed, promoter_bed, coding_transcript_bed
        ]
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         2,
        mem_gb:            8,
        disk_gb:           50,
        boot_disk_gb:      10,
        preemptible_tries: 1,
        max_retries:       1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu:            select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:         select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks:          "local-disk " + select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb,     default_attr.boot_disk_gb])
        docker:         docker
        preemptible:    select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:     select_first([runtime_attr.max_retries,       default_attr.max_retries])
    }
}

# ======================================================================
# TASK 3 — Bedtools intersect: SV vs each annotation BED
# ======================================================================
task Task3_BedtoolsIntersect {
    input {
        File   sv_bed
        File   transcript_bed
        File   cds_bed
        File   intron_bed
        File   utr3_bed
        File   utr5_bed
        File   promoter_bed
        File   coding_transcript_bed
        String sv_gtf_prefix
        String seed_suffix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        SV_COL5="<(zcat ~{sv_bed} | cut -f1-5)"

        bedtools intersect -wo \
            -a <(zcat ~{sv_bed} | cut -f1-5) \
            -b ~{transcript_bed} \
            | bgzip > ~{sv_gtf_prefix}.transcript.~{seed_suffix}.bed.gz

        bedtools intersect -wo \
            -a <(zcat ~{sv_bed} | cut -f1-5) \
            -b ~{cds_bed} \
            | bgzip > ~{sv_gtf_prefix}.CDS.~{seed_suffix}.bed.gz

        bedtools intersect -wo \
            -a <(zcat ~{sv_bed} | cut -f1-5) \
            -b ~{intron_bed} \
            | bgzip > ~{sv_gtf_prefix}.intron.~{seed_suffix}.bed.gz

        bedtools intersect -wo \
            -a <(zcat ~{sv_bed} | cut -f1-5) \
            -b ~{utr3_bed} \
            | bgzip > ~{sv_gtf_prefix}.utr_3.~{seed_suffix}.bed.gz

        bedtools intersect -wo \
            -a <(zcat ~{sv_bed} | cut -f1-5) \
            -b ~{utr5_bed} \
            | bgzip > ~{sv_gtf_prefix}.utr_5.~{seed_suffix}.bed.gz

        bedtools intersect -wo \
            -a <(zcat ~{sv_bed} | cut -f1-5) \
            -b ~{promoter_bed} \
            | bgzip > ~{sv_gtf_prefix}.promoter.~{seed_suffix}.bed.gz

        bedtools intersect -wo \
            -a <(zcat ~{sv_bed} | cut -f1-5) \
            -b ~{coding_transcript_bed} \
            | bgzip > ~{sv_gtf_prefix}.coding_transcript.~{seed_suffix}.bed.gz
    >>>

    output {
        File transcript_isec        = sv_gtf_prefix + ".transcript."        + seed_suffix + ".bed.gz"
        File cds_isec               = sv_gtf_prefix + ".CDS."               + seed_suffix + ".bed.gz"
        File intron_isec            = sv_gtf_prefix + ".intron."            + seed_suffix + ".bed.gz"
        File utr3_isec              = sv_gtf_prefix + ".utr_3."             + seed_suffix + ".bed.gz"
        File utr5_isec              = sv_gtf_prefix + ".utr_5."             + seed_suffix + ".bed.gz"
        File promoter_isec          = sv_gtf_prefix + ".promoter."          + seed_suffix + ".bed.gz"
        File coding_transcript_isec = sv_gtf_prefix + ".coding_transcript." + seed_suffix + ".bed.gz"
        Array[File] all_isec = [
            transcript_isec, cds_isec, intron_isec,
            utr3_isec, utr5_isec, promoter_isec, coding_transcript_isec
        ]
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         4,
        mem_gb:            16,
        disk_gb:           100,
        boot_disk_gb:      10,
        preemptible_tries: 1,
        max_retries:       1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu:            select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:         select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks:          "local-disk " + select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb,     default_attr.boot_disk_gb])
        docker:         docker
        preemptible:    select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:     select_first([runtime_attr.max_retries,       default_attr.max_retries])
    }
}

# ======================================================================
# TASK 4 — Extract overlap categories (awk/cut filters)
# Column layout after bedtools intersect -wo (SV cols 1-5, annot cols 6-12):
#   $1  chr_sv   $2  sv_start  $3  sv_end  $4  SVID  $5  SVTYPE
#   $6  chr_ann  $7  ann_start $8  ann_end $9  strand
#   $10 gene_id  $11 gene_type $12 gene_name  $13 overlap_bp
# ======================================================================
task Task4_ExtractOverlaps {
    input {
        File   transcript_isec
        File   cds_isec
        File   intron_isec
        File   utr3_isec
        File   utr5_isec
        File   promoter_isec
        File   coding_transcript_isec
        String sv_gtf_prefix
        String seed_suffix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String p = sv_gtf_prefix  # shorthand

    command <<<
        set -euo pipefail

        # ── whole transcript overlap: SV fully spans transcript ──────────
        zcat ~{transcript_isec} \
            | awk '$2 < $7 && $3 > $8' \
            | cut -f4,5,10,11 \
            > ~{p}.whole_transcript_overlap.~{seed_suffix}

        # ── SVs completely inside transcript body ────────────────────────
        zcat ~{transcript_isec} \
            | awk '$2 > $7 && $3 < $8' \
            | cut -f4,5,10,11 \
            > ~{p}.SVs_inside_transcripts.~{seed_suffix}

        # ── TSS overlap: SV reaches into the TSS end of transcript ───────
        zcat ~{transcript_isec} \
            | awk '$9=="+" && $2 < $7+1 && $3 < $8+1' \
            | cut -f4,5,10,11 \
            > ~{p}.tss_transcripts_overlap.~{seed_suffix}
        zcat ~{transcript_isec} \
            | awk '$9=="-" && $2 > $7-1 && $3 > $8-1' \
            | cut -f4,5,10,11 \
            >> ~{p}.tss_transcripts_overlap.~{seed_suffix}

        # ── partial transcript overlap: SV reaches into the distal end ───
        zcat ~{transcript_isec} \
            | awk '$9=="-" && $2 < $7+1 && $3 < $8+1' \
            | cut -f4,5,10,11 \
            > ~{p}.partial_transcripts_overlap.~{seed_suffix}
        zcat ~{transcript_isec} \
            | awk '$9=="+" && $2 > $7-1 && $3 > $8-1' \
            | cut -f4,5,10,11 \
            >> ~{p}.partial_transcripts_overlap.~{seed_suffix}

        # ── 5' UTR overlap (strand-aware) ────────────────────────────────
        zcat ~{utr5_isec} \
            | awk '$9=="+" && $3 < $8+1 && $3 > $7-1' \
            | cut -f4,5,10,11 \
            > ~{p}.5_prime_utr.~{seed_suffix}
        zcat ~{utr5_isec} \
            | awk '$9=="-" && $2 < $8+1 && $2 > $7-1' \
            | cut -f4,5,10,11 \
            >> ~{p}.5_prime_utr.~{seed_suffix}

        # ── 3' UTR overlap (strand-aware) ────────────────────────────────
        zcat ~{utr3_isec} \
            | awk '$9=="+" && $2 < $8+1 && $2 > $7-1' \
            | cut -f4,5,10,11 \
            > ~{p}.3_prime_utr.~{seed_suffix}
        zcat ~{utr3_isec} \
            | awk '$9=="-" && $3 < $8+1 && $3 > $7-1' \
            | cut -f4,5,10,11 \
            >> ~{p}.3_prime_utr.~{seed_suffix}

        # ── SV inside exon (CDS) ─────────────────────────────────────────
        zcat ~{cds_isec} \
            | awk '$2 > $7-1 && $3 < $8+1' \
            | cut -f4,5,10,11 \
            > ~{p}.inside_exons.~{seed_suffix}

        # ── SV inside intron ─────────────────────────────────────────────
        zcat ~{intron_isec} \
            | awk '$2 > $7-1 && $3 < $8+1' \
            | cut -f4,5,10,11 \
            > ~{p}.inside_introns.~{seed_suffix}

        # ── promoter overlap (all overlapping rows) ───────────────────────
        zcat ~{promoter_isec} \
            | cut -f4,5,10,11 \
            > ~{p}.promoter.~{seed_suffix}

        # ── TSS overlap for coding transcripts ────────────────────────────
        zcat ~{coding_transcript_isec} \
            | awk '$9=="+" && $2 < $7+1 && $3 < $8+1' \
            | cut -f4,5,10,11 \
            > ~{p}.tss_coding_transcripts_overlap.~{seed_suffix}
        zcat ~{coding_transcript_isec} \
            | awk '$9=="-" && $2 > $7-1 && $3 > $8-1' \
            | cut -f4,5,10,11 \
            >> ~{p}.tss_coding_transcripts_overlap.~{seed_suffix}

        # ── partial overlap for coding transcripts ────────────────────────
        zcat ~{coding_transcript_isec} \
            | awk '$9=="-" && $2 < $7+1 && $3 < $8+1' \
            | cut -f4,5,10,11 \
            > ~{p}.partial_coding_transcripts_overlap.~{seed_suffix}
        zcat ~{coding_transcript_isec} \
            | awk '$9=="+" && $2 > $7-1 && $3 > $8-1' \
            | cut -f4,5,10,11 \
            >> ~{p}.partial_coding_transcripts_overlap.~{seed_suffix}
    >>>

    output {
        File whole_transcript_overlap          = p + ".whole_transcript_overlap."          + seed_suffix
        File svs_inside_transcripts            = p + ".SVs_inside_transcripts."            + seed_suffix
        File tss_transcripts_overlap           = p + ".tss_transcripts_overlap."           + seed_suffix
        File partial_transcripts_overlap       = p + ".partial_transcripts_overlap."       + seed_suffix
        File utr5_overlap                      = p + ".5_prime_utr."                       + seed_suffix
        File utr3_overlap                      = p + ".3_prime_utr."                       + seed_suffix
        File inside_exons                      = p + ".inside_exons."                      + seed_suffix
        File inside_introns                    = p + ".inside_introns."                    + seed_suffix
        File promoter                  = p + ".promoter."                  + seed_suffix
        File tss_coding_transcripts_overlap    = p + ".tss_coding_transcripts_overlap."   + seed_suffix
        File partial_coding_transcripts_overlap = p + ".partial_coding_transcripts_overlap." + seed_suffix
        Array[File] all_overlaps = [
            whole_transcript_overlap, svs_inside_transcripts,
            tss_transcripts_overlap, partial_transcripts_overlap,
            utr5_overlap, utr3_overlap, inside_exons, inside_introns,
            promoter, tss_coding_transcripts_overlap,
            partial_coding_transcripts_overlap
        ]
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         2,
        mem_gb:            8,
        disk_gb:           50,
        boot_disk_gb:      10,
        preemptible_tries: 1,
        max_retries:       1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu:            select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:         select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks:          "local-disk " + select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb,     default_attr.boot_disk_gb])
        docker:         docker
        preemptible:    select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:     select_first([runtime_attr.max_retries,       default_attr.max_retries])
    }
}

# ======================================================================
# TASK 5 — Categorize intact vs partial exon overlap (R)
# ======================================================================
task Task5_CategorizeExonOverlap {
    input {
        File   cds_isec
        File   svs_inside_transcripts
        File   r_script
        String sv_gtf_prefix
        String seed_suffix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        Rscript ~{r_script} \
            -c ~{cds_isec} \
            -g ~{svs_inside_transcripts} \
            -p ~{sv_gtf_prefix}.~{seed_suffix}
    >>>

    output {
        File intact_exon_overlap  = "~{sv_gtf_prefix}.~{seed_suffix}.intact_exon_overlap"
        File partial_exon_overlap = "~{sv_gtf_prefix}.~{seed_suffix}.partial_exon_overlap"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         2,
        mem_gb:            16,
        disk_gb:           50,
        boot_disk_gb:      10,
        preemptible_tries: 1,
        max_retries:       1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu:            select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:         select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks:          "local-disk " + select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb,     default_attr.boot_disk_gb])
        docker:         docker
        preemptible:    select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:     select_first([runtime_attr.max_retries,       default_attr.max_retries])
    }
}

# ======================================================================
# TASK 6 — Reorganize SVID vs gene for a single overlap file (R)
# ======================================================================
task Task6_ReorganizeSVIDGene {
    input {
        File   transcript_bed
        File   overlap_file
        String output_name
        File   r_script
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        Rscript ~{r_script} \
            -g ~{transcript_bed} \
            -i <(cut -f1-4 ~{overlap_file}) \
            -o ~{output_name}
    >>>

    output {
        File reorganized = output_name
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         2,
        mem_gb:            16,
        disk_gb:           50,
        boot_disk_gb:      10,
        preemptible_tries: 1,
        max_retries:       1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu:            select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:         select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks:          "local-disk " + select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb,     default_attr.boot_disk_gb])
        docker:         docker
        preemptible:    select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:     select_first([runtime_attr.max_retries,       default_attr.max_retries])
    }
}

# ======================================================================
# TASK 7 — Integrate all reorganized overlaps (R)
# The R script reads all *.reorganized files by prefix+annotation.
# We symlink them into the working directory first.
# ======================================================================
task Task7_IntegrateOverlaps {
    input {
        Array[File] all_reorganized
        File        r_script
        String      sv_gtf_prefix
        String      seed_suffix
        String      docker
        RuntimeAttr? runtime_attr_override
    }

    String out_file = sv_gtf_prefix + "." + seed_suffix + ".integrated"

    command <<<
        set -euo pipefail

        # Symlink all reorganized files into working dir
        for f in ~{sep=' ' all_reorganized}; do
            ln -sf "$f" "$(basename $f)"
        done

        Rscript ~{r_script} \
            -p ~{sv_gtf_prefix} \
            -a ~{seed_suffix} \
            -o ~{out_file}
    >>>

    output {
        File integrated_file = out_file
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         2,
        mem_gb:            16,
        disk_gb:           50,
        boot_disk_gb:      10,
        preemptible_tries: 1,
        max_retries:       1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu:            select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:         select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks:          "local-disk " + select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb,     default_attr.boot_disk_gb])
        docker:         docker
        preemptible:    select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:     select_first([runtime_attr.max_retries,       default_attr.max_retries])
    }
}

# ======================================================================
# TASK 8 — Calculate per-gene SV data (R) — final output
# ======================================================================
task Task8_CalcuGeneData {
    input {
        File r_script
        String seed_suffix
        File integrated_file
        File sv_info
        File gene_info
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(integrated_file, ".integrated")

    command <<<
        set -euo pipefail

        Rscript ~{r_script} \
        -p ~{seed_suffix} \
        -s ~{sv_info} \
        -g ~{gene_info} \
        -r ~{integrated_file} \
        -o ~{prefix}.rData
    >>>

    output {
        File result        = "~{prefix}.rData"
        File result_tsv_gz = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         2,
        mem_gb:            16,
        disk_gb:           50,
        boot_disk_gb:      10,
        preemptible_tries: 1,
        max_retries:       1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu:            select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:         select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks:          "local-disk " + select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb,     default_attr.boot_disk_gb])
        docker:         docker
        preemptible:    select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:     select_first([runtime_attr.max_retries,       default_attr.max_retries])
    }
}
