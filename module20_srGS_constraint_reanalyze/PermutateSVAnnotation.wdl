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

workflow PermutateSVAnnotation {

    input {
        File   gtf_file
        Int    permu_number
        File   sv_bed
        File   tel_cen_bed
        File   gene_info

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
            permute_script        = permute_gtf_script,
            gtf_label             = gtf_label,
            seed_suffix           = seed_suffix,
            docker                = python_docker,
            runtime_attr_override = runtime_attr_permute_gtf
    }

    # ── Task 2: Split permuted GTF into 7 annotation BED files ───────────
    call Task2_SplitGTF {
        input:
            permuted_gtf          = Task1_PermuteGTF.permuted_gtf,
            split_script          = split_gtf_script,
            gtf_label             = gtf_label,
            seed_suffix           = seed_suffix,
            docker                = python_docker,
            runtime_attr_override = runtime_attr_split_gtf
    }

    # ── Task 3: Bedtools intersect SV vs each annotation BED ─────────────
    call Task3_BedtoolsIntersect {
        input:
            sv_bed                = sv_bed,
            transcript_bed        = Task2_SplitGTF.transcript_bed,
            cds_bed               = Task2_SplitGTF.cds_bed,
            intron_bed            = Task2_SplitGTF.intron_bed,
            utr3_bed              = Task2_SplitGTF.utr3_bed,
            utr5_bed              = Task2_SplitGTF.utr5_bed,
            promoter_bed          = Task2_SplitGTF.promoter_bed,
            coding_transcript_bed = Task2_SplitGTF.coding_transcript_bed,
            sv_gtf_prefix         = sv_gtf_prefix,
            seed_suffix           = seed_suffix,
            docker                = bedtools_docker,
            runtime_attr_override = runtime_attr_bedtools_intersect
    }

    # ── Task 4: Extract overlap categories via awk/cut ────────────────────
    call Task4_ExtractOverlaps {
        input:
            transcript_isec        = Task3_BedtoolsIntersect.transcript_isec,
            cds_isec               = Task3_BedtoolsIntersect.cds_isec,
            intron_isec            = Task3_BedtoolsIntersect.intron_isec,
            utr3_isec              = Task3_BedtoolsIntersect.utr3_isec,
            utr5_isec              = Task3_BedtoolsIntersect.utr5_isec,
            promoter_isec          = Task3_BedtoolsIntersect.promoter_isec,
            coding_transcript_isec = Task3_BedtoolsIntersect.coding_transcript_isec,
            sv_gtf_prefix          = sv_gtf_prefix,
            seed_suffix            = seed_suffix,
            docker                 = bedtools_docker,
            runtime_attr_override  = runtime_attr_extract_overlaps
    }

    # ── Task 5: Categorize intact vs partial exon overlap (R) ────────────
    call Task5_CategorizeExonOverlap {
        input:
            cds_isec               = Task3_BedtoolsIntersect.cds_isec,
            svs_inside_transcripts = Task4_ExtractOverlaps.svs_inside_transcripts,
            r_script               = categorize_r_script,
            sv_gtf_prefix          = sv_gtf_prefix,
            seed_suffix            = seed_suffix,
            docker                 = r_docker,
            runtime_attr_override  = runtime_attr_categorize_exon
    }

    # ── Task 6: Reorganize SVID vs gene — one call per overlap type (R) ────
    call Task6_ReorganizeSVIDGene as Task6_Reorg_WholeTranscript {
        input:
            transcript_bed        = Task2_SplitGTF.transcript_bed,
            overlap_file          = Task4_ExtractOverlaps.whole_transcript_overlap,
            output_name           = sv_gtf_prefix + ".whole_transcript_overlap." + seed_suffix + ".reorganized",
            r_script              = reorganize_r_script,
            docker                = r_docker,
            runtime_attr_override = runtime_attr_reorganize_svid
    }
    call Task6_ReorganizeSVIDGene as Task6_Reorg_Utr3 {
        input:
            transcript_bed        = Task2_SplitGTF.transcript_bed,
            overlap_file          = Task4_ExtractOverlaps.utr3_overlap,
            output_name           = sv_gtf_prefix + ".3_prime_utr." + seed_suffix + ".reorganized",
            r_script              = reorganize_r_script,
            docker                = r_docker,
            runtime_attr_override = runtime_attr_reorganize_svid
    }
    call Task6_ReorganizeSVIDGene as Task6_Reorg_Utr5 {
        input:
            transcript_bed        = Task2_SplitGTF.transcript_bed,
            overlap_file          = Task4_ExtractOverlaps.utr5_overlap,
            output_name           = sv_gtf_prefix + ".5_prime_utr." + seed_suffix + ".reorganized",
            r_script              = reorganize_r_script,
            docker                = r_docker,
            runtime_attr_override = runtime_attr_reorganize_svid
    }
    call Task6_ReorganizeSVIDGene as Task6_Reorg_IntactExon {
        input:
            transcript_bed        = Task2_SplitGTF.transcript_bed,
            overlap_file          = Task5_CategorizeExonOverlap.intact_exon_overlap,
            output_name           = sv_gtf_prefix + ".intact_exon_overlap." + seed_suffix + ".reorganized",
            r_script              = reorganize_r_script,
            docker                = r_docker,
            runtime_attr_override = runtime_attr_reorganize_svid
    }
    call Task6_ReorganizeSVIDGene as Task6_Reorg_PartialExon {
        input:
            transcript_bed        = Task2_SplitGTF.transcript_bed,
            overlap_file          = Task5_CategorizeExonOverlap.partial_exon_overlap,
            output_name           = sv_gtf_prefix + ".partial_exon_overlap." + seed_suffix + ".reorganized",
            r_script              = reorganize_r_script,
            docker                = r_docker,
            runtime_attr_override = runtime_attr_reorganize_svid
    }
    call Task6_ReorganizeSVIDGene as Task6_Reorg_TssTranscripts {
        input:
            transcript_bed        = Task2_SplitGTF.transcript_bed,
            overlap_file          = Task4_ExtractOverlaps.tss_transcripts_overlap,
            output_name           = sv_gtf_prefix + ".tss_transcripts_overlap." + seed_suffix + ".reorganized",
            r_script              = reorganize_r_script,
            docker                = r_docker,
            runtime_attr_override = runtime_attr_reorganize_svid
    }
    call Task6_ReorganizeSVIDGene as Task6_Reorg_PartialTranscripts {
        input:
            transcript_bed        = Task2_SplitGTF.transcript_bed,
            overlap_file          = Task4_ExtractOverlaps.partial_transcripts_overlap,
            output_name           = sv_gtf_prefix + ".partial_transcripts_overlap." + seed_suffix + ".reorganized",
            r_script              = reorganize_r_script,
            docker                = r_docker,
            runtime_attr_override = runtime_attr_reorganize_svid
    }
    call Task6_ReorganizeSVIDGene as Task6_Reorg_InsideExons {
        input:
            transcript_bed        = Task2_SplitGTF.transcript_bed,
            overlap_file          = Task4_ExtractOverlaps.inside_exons,
            output_name           = sv_gtf_prefix + ".inside_exons." + seed_suffix + ".reorganized",
            r_script              = reorganize_r_script,
            docker                = r_docker,
            runtime_attr_override = runtime_attr_reorganize_svid
    }
    call Task6_ReorganizeSVIDGene as Task6_Reorg_InsideIntrons {
        input:
            transcript_bed        = Task2_SplitGTF.transcript_bed,
            overlap_file          = Task4_ExtractOverlaps.inside_introns,
            output_name           = sv_gtf_prefix + ".inside_introns." + seed_suffix + ".reorganized",
            r_script              = reorganize_r_script,
            docker                = r_docker,
            runtime_attr_override = runtime_attr_reorganize_svid
    }
    call Task6_ReorganizeSVIDGene as Task6_Reorg_Promoter {
        input:
            transcript_bed        = Task2_SplitGTF.transcript_bed,
            overlap_file          = Task4_ExtractOverlaps.promoter,
            output_name           = sv_gtf_prefix + ".promoter." + seed_suffix + ".reorganized",
            r_script              = reorganize_r_script,
            docker                = r_docker,
            runtime_attr_override = runtime_attr_reorganize_svid
    }

    # ── Task 7: Integrate all reorganized overlaps (R) ────────────────────
    call Task7_IntegrateOverlaps {
        input:
            all_reorganized = [
                Task6_Reorg_WholeTranscript.reorganized,
                Task6_Reorg_Utr3.reorganized,
                Task6_Reorg_Utr5.reorganized,
                Task6_Reorg_IntactExon.reorganized,
                Task6_Reorg_PartialExon.reorganized,
                Task6_Reorg_TssTranscripts.reorganized,
                Task6_Reorg_PartialTranscripts.reorganized,
                Task6_Reorg_InsideExons.reorganized,
                Task6_Reorg_InsideIntrons.reorganized,
                Task6_Reorg_Promoter.reorganized
            ],
            r_script        = integrate_r_script,
            sv_gtf_prefix   = sv_gtf_prefix,
            seed_suffix     = seed_suffix,
            docker          = r_docker,
            runtime_attr_override = runtime_attr_integrate_overlaps
    }

    # ── Task 8: Calculate per-gene SV data (R) ────────────────────────────
    call Task8_CalcuGeneData {
        input:
            r_script         = calcu_r_script,
            seed_suffix      = seed_suffix,
            integrated_file  = Task7_IntegrateOverlaps.integrated_file,
            sv_info = sv_bed,
            gene_info  = gene_info, 
            docker           = r_docker,
            runtime_attr_override = runtime_attr_calcu_gene_data
    }

    output {
        File   permuted_gtf             = Task1_PermuteGTF.permuted_gtf
        Array[File] split_beds          = Task2_SplitGTF.all_beds
        Array[File] intersection_beds   = Task3_BedtoolsIntersect.all_isec
        Array[File] overlap_tables      = Task4_ExtractOverlaps.all_overlaps
        Array[File] reorganized         = [
            Task6_Reorg_WholeTranscript.reorganized,
            Task6_Reorg_Utr3.reorganized,
            Task6_Reorg_Utr5.reorganized,
            Task6_Reorg_IntactExon.reorganized,
            Task6_Reorg_PartialExon.reorganized,
            Task6_Reorg_TssTranscripts.reorganized,
            Task6_Reorg_PartialTranscripts.reorganized,
            Task6_Reorg_InsideExons.reorganized,
            Task6_Reorg_InsideIntrons.reorganized,
            Task6_Reorg_Promoter.reorganized
        ]
        File        integrated          = Task7_IntegrateOverlaps.integrated_file
        File        result              = Task8_CalcuGeneData.result
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
        File   permute_script
        String gtf_label
        String seed_suffix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String out_gtf = gtf_label + "." + seed_suffix + ".gtf.gz"

    command <<<
        python3 ~{permute_script} \
            ~{permu_number} \
            ~{gtf_file} \
            ~{tel_cen_bed} \
            ~{out_gtf}
    >>>

    output {
        File permuted_gtf = out_gtf
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         4,
        mem_gb:            8,
        disk_gb:           50,
        boot_disk_gb:      10,
        preemptible_tries: 3,
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
        preemptible_tries: 3,
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
        preemptible_tries: 3,
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
        preemptible_tries: 3,
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
        preemptible_tries: 3,
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
        preemptible_tries: 3,
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
        preemptible_tries: 3,
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
        File result = "~{prefix}.rData"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:         2,
        mem_gb:            16,
        disk_gb:           50,
        boot_disk_gb:      10,
        preemptible_tries: 3,
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
