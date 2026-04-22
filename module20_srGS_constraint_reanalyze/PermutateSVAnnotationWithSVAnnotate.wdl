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
        File contigs_list

        # Python scripts
        File   permute_gtf_script       # permute_gtf.py
        File   split_gtf_script         # split_gtf_annotations.py
        File calculate_gene_data_script

        String gtf_label                # e.g. "r3.gencode.v39.ensembl.105"
        String sv_label                 # e.g. "gnomAD_SV_v3"

        # Docker images (override as needed)
        String python_docker   = "python:3.11-slim"
        String bedtools_docker = "quay.io/biocontainers/bedtools:2.31.1--h63f29b4_0"
        String r_docker        = "r-base:4.3.0"
        String gatk_docker =  "us.gcr.io/broad-dsde-methods/gatk-sv/gatk:2025-05-20-4.6.2.0-4-g1facd911e-NIGHTLY-SNAPSHOT"
        String utils_docker = "us-central1-docker.pkg.dev/talkowski-training/kj-development/utils:kj_V13"
        String sv_pipeline_docker = "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2026-02-06-v1.1-797b7604"

        # Runtime attribute overrides per task
        RuntimeAttr? runtime_attr_permute_gtf
        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_annotate_sv
        RuntimeAttr? runtime_attr_merge_vcf
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
    Array[String] contigs = read_lines(contigs_list)

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


    # ── Task 2: Split VCF by contig ─────────────────────────────────────────
    scatter (contig in contigs) {
        call SplitVcfByContig {
            input:
                vcf                   = vcf,
                vcf_idx               = vcf_idx,
                contig                = contig,
                prefix                = "permu_~{permu_number}",
                docker                = gatk_docker,
                runtime_attr_override = runtime_attr_split_vcf
        }

        # ── Task 3: Annotate functional consequences per contig ────────────
        call AnnotateFunctionalConsequences {
            input:
                vcf                   = SplitVcfByContig.contig_vcf,
                vcf_idx               = SplitVcfByContig.contig_vcf_idx,
                prefix                = "permu_~{permu_number}.~{contig}",
                coding_gtf            = Task1_PermuteGTF.permuted_gtf,
                docker                = gatk_docker,
                runtime_attr_override = runtime_attr_annotate_sv
        }
    }

    # ── Task 4: Merge annotated VCFs across contigs ─────────────────────────
    call MergeAnnotatedVcfs {
        input:
            vcfs                  = AnnotateFunctionalConsequences.anno_vcf,
            prefix                = "permu_~{permu_number}",
            docker                = gatk_docker,
            runtime_attr_override = runtime_attr_merge_vcf
    }

    call CalculateGeneData {
        input:
            vcf = MergeAnnotatedVcfs.merged_vcf,
            vcf_idx = MergeAnnotatedVcfs.merged_vcf_idx,
            script = calculate_gene_data_script,
            gene_info = gene_info,
            prefix = "permu_~{permu_number}",
            docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_calcu_gene_data
    }

    output {
        File permuted_gtf          = Task1_PermuteGTF.permuted_gtf
        File too_large_genes_list  = Task1_PermuteGTF.too_large_genes_list
        Array[File] contig_vcfs    = SplitVcfByContig.contig_vcf
        Array[File] contig_vcf_idx = SplitVcfByContig.contig_vcf_idx
        Array[File] anno_vcfs      = AnnotateFunctionalConsequences.anno_vcf
        Array[File] anno_vcf_idx   = AnnotateFunctionalConsequences.anno_vcf_idx
        File merged_anno_vcf       = MergeAnnotatedVcfs.merged_vcf
        File merged_anno_vcf_idx   = MergeAnnotatedVcfs.merged_vcf_idx
        File gene_data = CalculateGeneData.result
    }

}

task SplitVcfByContig {
    input {
        File vcf
        File vcf_idx
        String contig
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

        gatk --java-options "-Xmx~{java_mem_mb}m" SelectVariants \
            -V ~{vcf} \
            -L ~{contig} \
            -O ~{prefix}.~{contig}.vcf.gz
    >>>

    output {
        File contig_vcf = "~{prefix}.~{contig}.vcf.gz"
        File contig_vcf_idx = "~{prefix}.~{contig}.vcf.gz.tbi"
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

        gunzip -c ~{coding_gtf} > coding.gtf

        gatk --java-options "-Xmx~{java_mem_mb}m" SVAnnotate \
            -V ~{vcf} \
            --protein-coding-gtf coding.gtf \
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

task MergeAnnotatedVcfs {
    input {
        Array[File] vcfs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Int java_mem_mb = 1000 * ceil(0.8 * select_first([runtime_attr.mem_gb, default_attr.mem_gb]))

    command <<<
        set -euo pipefail

        INPUT_ARGS=""
        for f in ~{sep=' ' vcfs}; do
            INPUT_ARGS+=" -I ${f}"
        done

        gatk --java-options "-Xmx~{java_mem_mb}m" MergeVcfs \
            ${INPUT_ARGS} \
            -O ~{prefix}.annotated.merged.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.annotated.merged.vcf.gz"
        File merged_vcf_idx = "~{prefix}.annotated.merged.vcf.gz.tbi"
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



task CalculateGeneData {
    input {
        File script
        String prefix
        File vcf
        File vcf_idx
        File gene_info
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python ~{script} \
            --vcf ~{vcf} \
            --gene-info ~{gene_info} \
            --out ~{prefix}.tsv.gz
    >>>

    output {
        File result = "~{prefix}.tsv.gz"
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
