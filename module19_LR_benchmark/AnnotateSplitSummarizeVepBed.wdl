version 1.0

## AnnotateSplitSummarizeVepBed
##
## Takes a gnomAD_LR.{cohort}.vep_parsed.bed.gz (produced by VcfToBedVepParsed,
## optionally already carrying extra INFO columns such as nhomref/nhet/nhomalt
## via its extra_info_fields option) and:
##
##   1. AnnotateOverlaps: adds 6 boolean columns flagging overlap with
##      telomere / centromere / N-masked / srGS-low-coverage / lrGS-blacklist
##      regions, plus a centromere-specific subset of the lrGS blacklist
##      (in_telomere, in_centromere, in_nmasked, in_srGS_lowcov,
##      in_lrGS_lowcov, in_centromere_blacklist).
##
##   2. AnalyzeVepBed: classifies every variant into at most one of
##      {pLoF, missense, synonymous, intronic, intergenic} using vep_IMPACT
##      (already the single most-severe tier per variant) as the primary,
##      inherently mutually-exclusive key -- see analyze_vep_bed.py for the
##      exact rule per category. Writes each category to its own bed.gz,
##      tallies per-gene variant counts (with SNV/indel/SV splits for pLoF
##      and intronic), tallies per-gene PREDICTED_LOF SV counts (DEL/INS,
##      with/without TRID) from SVAnnotate's independent PREDICTED_LOF INFO
##      field, and -- if nhomref/nhet/nhomalt are present -- converts all of
##      the above into average variants-per-genome (genome-wide and per gene)
##      using carrier_count = nhet + nhomalt per site and N_samples = the
##      median of (nhomref+nhet+nhomalt) across all sites.

workflow AnnotateSplitSummarizeVepBed {

    input {
        File    vep_bed              # gnomAD_LR.{cohort}.vep_parsed.bed.gz
        String  cohort                # e.g. "hprc_hgsvc" -- used to name outputs
        File    telomere_bed          # hg38.telomere.bed.gz
        File    centromere_bed        # hg38.centromere.bed.gz (has a header row)
        File    nmasked_bed           # hg38.n_masked_regions.bed.gz
        File    srgs_lowcov_bed       # gnomAD_srGS.coverage_bins.low4.bed.gz
        File    lrgs_blacklist_bed    # hgsvc_hprc.blacklist.p90.annotated.bed.gz (has a header row)
        File    annotate_overlaps_script  # annotate_overlaps.sh
        File    analyze_script            # analyze_vep_bed.py
        String  bedtools_docker = "staphb/bedtools:2.31.1"
        String  python_docker   = "python:3.11-slim"
        Int     mem_gb      = 16
        Int     disk_gb      = 500
        Int     preemptible  = 1
    }

    call AnnotateOverlaps {
        input:
            vep_bed                  = vep_bed,
            telomere_bed             = telomere_bed,
            centromere_bed           = centromere_bed,
            nmasked_bed               = nmasked_bed,
            srgs_lowcov_bed          = srgs_lowcov_bed,
            lrgs_blacklist_bed       = lrgs_blacklist_bed,
            annotate_overlaps_script = annotate_overlaps_script,
            cohort                   = cohort,
            docker                   = bedtools_docker,
            mem_gb                   = mem_gb,
            disk_gb                  = disk_gb,
            preemptible              = preemptible
    }

    call AnalyzeVepBed {
        input:
            annotated_bed  = AnnotateOverlaps.out_bed,
            analyze_script = analyze_script,
            cohort         = cohort,
            docker         = python_docker,
            mem_gb         = mem_gb,
            disk_gb        = disk_gb,
            preemptible    = preemptible
    }

    output {
        File annotated_bed                = AnnotateOverlaps.out_bed
        File plof_bed                      = AnalyzeVepBed.plof_bed
        File missense_bed                  = AnalyzeVepBed.missense_bed
        File synonymous_bed                = AnalyzeVepBed.synonymous_bed
        File intron_bed                    = AnalyzeVepBed.intron_bed
        File intergenic_bed                = AnalyzeVepBed.intergenic_bed
        File gene_variant_counts           = AnalyzeVepBed.gene_variant_counts
        File predicted_lof_sv_gene_counts  = AnalyzeVepBed.predicted_lof_sv_gene_counts
        File genome_wide_summary           = AnalyzeVepBed.genome_wide_summary
    }

    meta {
        author: "gnomAD LR analysis"
        description: "Annotate region overlaps, split by VEP category, and summarize gene/genome burden."
    }
}

task AnnotateOverlaps {
    input {
        File    vep_bed
        File    telomere_bed
        File    centromere_bed
        File    nmasked_bed
        File    srgs_lowcov_bed
        File    lrgs_blacklist_bed
        File    annotate_overlaps_script
        String  cohort
        String  docker
        Int     mem_gb
        Int     disk_gb
        Int     preemptible
    }

    String out_bed_name = "gnomAD_LR." + cohort + ".vep_parsed.annotated.bed.gz"

    command <<<
        set -euo pipefail
        bash ~{annotate_overlaps_script} \
            ~{vep_bed} ~{telomere_bed} ~{centromere_bed} ~{nmasked_bed} \
            ~{srgs_lowcov_bed} ~{lrgs_blacklist_bed} ~{out_bed_name}
    >>>

    output {
        File out_bed = out_bed_name
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         2
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

task AnalyzeVepBed {
    input {
        File    annotated_bed
        File    analyze_script
        String  cohort
        String  docker
        Int     mem_gb
        Int     disk_gb
        Int     preemptible
    }

    String out_prefix = "gnomAD_LR." + cohort + ".vep_parsed"

    command <<<
        set -euo pipefail
        python3 ~{analyze_script} ~{annotated_bed} ~{out_prefix}
    >>>

    output {
        File plof_bed                     = out_prefix + ".HIGH_impact.bed.gz"
        File missense_bed                 = out_prefix + ".missense.bed.gz"
        File synonymous_bed               = out_prefix + ".synonymous.bed.gz"
        File intron_bed                    = out_prefix + ".intron.bed.gz"
        File intergenic_bed                = out_prefix + ".intergenic.bed.gz"
        File gene_variant_counts          = out_prefix + ".gene_variant_counts.tsv"
        File predicted_lof_sv_gene_counts = out_prefix + ".predicted_lof_sv_gene_counts.tsv"
        File genome_wide_summary          = out_prefix + ".genome_wide_summary.tsv"
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         2
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}
