version 1.0

# ──────────────────────────────────────────────────────────────────────────
# ResolveGenotypeComplex.wdl
#
# Consolidated workflow that replaces the three separate workflows:
#   1. ResolveComplexVariants  (RCV)
#   2. GenotypeComplexVariants (GCV)
#   3. RefineComplexVariants   (RefCV)
#
# The logic is driven by the unified `gatk-sv-cpx` Python package with
# three subcommands:
#   gatk-sv-cpx resolve
#   gatk-sv-cpx evaluate-evidence
#   gatk-sv-cpx genotype-and-refine
#
# Workflow structure:
#   Step 0  – Build sample→batch→PE file map (once)
#   Step 1  – Uniform sharding of the input clustered VCFs
#   Step 2  – Per-shard: resolve → evaluate-evidence → genotype-and-refine
#   Step 3  – Gather finalized shard VCFs into the cohort-wide output
# ──────────────────────────────────────────────────────────────────────────

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow ResolveGenotypeComplex {
    input {
        String cohort_name

        # Clustered input VCFs (one per contig from ClusterBatch)
        Array[File] cluster_vcfs
        Array[File] cluster_vcf_indexes

        # SR lists for breakpoint overlap filtering
        Array[File] cluster_bothside_pass_lists
        Array[File] cluster_background_fail_lists

        # PE metrics (one per batch)
        Array[String] batch_name_list
        Array[File]   batch_sample_lists
        Array[File]   PE_metrics
        Array[File]   PE_metrics_indexes

        # Depth evidence (one per batch)
        Array[File] Depth_DEL_beds
        Array[File] Depth_DUP_beds

        # Depth genotyping inputs (one per batch)
        Array[File] bincov_files

        # Reference
        File contig_list
        File ref_dict

        # Sharding
        Int records_per_shard = 5000

        # Evidence thresholds
        Int min_pe_cpx = 3
        Int min_pe_ctx = 3
        Float depth_threshold = 0.5

        # Docker images
        String sv_pipeline_docker
        String sv_base_mini_docker

        # Runtime overrides
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_scatter_vcf
        RuntimeAttr? runtime_attr_resolve
        RuntimeAttr? runtime_attr_evaluate_evidence
        RuntimeAttr? runtime_attr_genotype_and_refine
        RuntimeAttr? runtime_attr_concat_vcfs
        RuntimeAttr? runtime_attr_build_pe_map
    }

    # ── Step 0: Build cohort-wide sample→batch→PE-file map ─────────────
    call BuildSampleBatchPEMap {
        input:
            batch_name_list   = batch_name_list,
            batch_sample_lists = batch_sample_lists,
            PE_metrics         = PE_metrics,
            prefix             = cohort_name,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_build_pe_map
    }

    # ── Step 1: Merge per-contig VCFs into one, then scatter uniformly ─
    call MiniTasks.ConcatVcfs as MergeClusteredVcfs {
        input:
            vcfs           = cluster_vcfs,
            vcfs_idx       = cluster_vcf_indexes,
            allow_overlaps = true,
            outfile_prefix = "~{cohort_name}.clustered_merged",
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_merge_vcfs
    }

    # Concatenate SR lists genome-wide
    call MiniTasks.CatUncompressedFiles as MergeBothsidePass {
        input:
            shards       = cluster_bothside_pass_lists,
            outfile_name = "~{cohort_name}.all.sr_bothside_pass.txt",
            sv_base_mini_docker = sv_base_mini_docker
    }
    call MiniTasks.CatUncompressedFiles as MergeBackgroundFail {
        input:
            shards       = cluster_background_fail_lists,
            outfile_name = "~{cohort_name}.all.sr_background_fail.txt",
            sv_base_mini_docker = sv_base_mini_docker
    }

    call MiniTasks.ScatterVcf {
        input:
            vcf              = MergeClusteredVcfs.concat_vcf,
            prefix           = "~{cohort_name}.shard",
            records_per_shard = records_per_shard,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_scatter_vcf
    }

    # ── Step 2: Per-shard resolution + evidence + refinement ───────────
    scatter (shard_vcf in ScatterVcf.shards) {

        # 2a. Resolve complex SV linkages
        call ResolveShard {
            input:
                input_vcf        = shard_vcf,
                bothside_pass    = MergeBothsidePass.outfile,
                background_fail  = MergeBackgroundFail.outfile,
                variant_prefix   = cohort_name,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_resolve
        }

        # 2b. Gather PE + RD evidence
        call EvaluateEvidenceShard {
            input:
                input_vcf     = ResolveShard.resolved_vcf,
                pe_metrics    = BuildSampleBatchPEMap.sample_batch_pe_map,
                depth_support = BuildSampleBatchPEMap.sample_batch_pe_map,  # placeholder
                min_pe_cpx    = min_pe_cpx,
                min_pe_ctx    = min_pe_ctx,
                depth_threshold = depth_threshold,
                prefix        = basename(shard_vcf, ".vcf.gz"),
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_evaluate_evidence
        }

        # 2c. Genotype and refine
        call GenotypeAndRefineShard {
            input:
                input_vcf     = ResolveShard.resolved_vcf,
                evidence_json = EvaluateEvidenceShard.evidence_json,
                prefix        = basename(shard_vcf, ".vcf.gz"),
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_genotype_and_refine
        }
    }

    # ── Step 3: Gather final VCF ───────────────────────────────────────
    call MiniTasks.ConcatVcfs as GatherFinalVcf {
        input:
            vcfs           = GenotypeAndRefineShard.finalized_vcf,
            vcfs_idx       = GenotypeAndRefineShard.finalized_vcf_index,
            allow_overlaps = true,
            outfile_prefix = "~{cohort_name}.cpx_final",
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat_vcfs
    }

    # ── Outputs ────────────────────────────────────────────────────────
    output {
        File cpx_final_vcf       = GatherFinalVcf.concat_vcf
        File cpx_final_vcf_index = GatherFinalVcf.concat_vcf_idx

        File bothside_pass_list  = MergeBothsidePass.outfile
        File background_fail_list = MergeBackgroundFail.outfile
    }
}


# ══════════════════════════════════════════════════════════════════════════
# Task definitions
# ══════════════════════════════════════════════════════════════════════════

task BuildSampleBatchPEMap {
    input {
        Array[String] batch_name_list
        Array[File]   batch_sample_lists
        Array[File]   PE_metrics
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail
        python3 <<CODE
import os

batch_sample_lists = ["~{sep='", "' batch_sample_lists}"]
batch_name_list = ["~{sep='", "' batch_name_list}"]
pe_metrics = ["~{sep='", "' PE_metrics}"]

with open("~{prefix}.sample_batch_pe_map.tsv", "w") as out:
    for i, batch in enumerate(batch_name_list):
        pe_file = os.path.basename(pe_metrics[i])
        with open(batch_sample_lists[i]) as inp:
            for line in inp:
                sample = line.strip()
                if sample:
                    out.write(f"{sample}\t{batch}\t{pe_file}\n")
CODE
    >>>

    output {
        File sample_batch_pe_map = "~{prefix}.sample_batch_pe_map.tsv"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}


task ResolveShard {
    input {
        File input_vcf
        File bothside_pass
        File background_fail
        String variant_prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(input_vcf, "GiB")
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(10 + input_size * 4),
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(input_vcf, ".vcf.gz")

    command <<<
        set -euo pipefail
        gatk-sv-cpx resolve \
            -i ~{input_vcf} \
            -o ~{prefix}.resolved.vcf.gz \
            --bothside-pass ~{bothside_pass} \
            --background-fail ~{background_fail} \
            --variant-prefix ~{variant_prefix}
    >>>

    output {
        File resolved_vcf       = "~{prefix}.resolved.vcf.gz"
        File resolved_vcf_index = "~{prefix}.resolved.vcf.gz.tbi"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}


task EvaluateEvidenceShard {
    input {
        File input_vcf
        File? pe_metrics
        File? depth_support
        Int min_pe_cpx
        Int min_pe_ctx
        Float depth_threshold
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(input_vcf, "GiB")
    Float pe_size = if defined(pe_metrics) then size(select_first([pe_metrics]), "GiB") else 0
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(10 + input_size * 2 + pe_size * 2),
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail
        gatk-sv-cpx evaluate-evidence \
            -i ~{input_vcf} \
            -o ~{prefix}.evidence.json \
            ~{"--pe-metrics " + pe_metrics} \
            ~{"--depth-support " + depth_support} \
            --min-pe-cpx ~{min_pe_cpx} \
            --min-pe-ctx ~{min_pe_ctx} \
            --depth-threshold ~{depth_threshold}
    >>>

    output {
        File evidence_json = "~{prefix}.evidence.json"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}


task GenotypeAndRefineShard {
    input {
        File input_vcf
        File evidence_json
        File? unresolved_ids
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(input_vcf, "GiB")
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(10 + input_size * 4),
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail
        gatk-sv-cpx genotype-and-refine \
            -i ~{input_vcf} \
            --evidence ~{evidence_json} \
            -o ~{prefix}.finalized.vcf.gz \
            ~{"--unresolved-ids " + unresolved_ids}
    >>>

    output {
        File finalized_vcf       = "~{prefix}.finalized.vcf.gz"
        File finalized_vcf_index = "~{prefix}.finalized.vcf.gz.tbi"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
