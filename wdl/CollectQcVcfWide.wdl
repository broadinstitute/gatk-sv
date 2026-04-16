version 1.0

import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow CollectQcVcfWide {
    input {
        File vcf
        File? vcf_idx
        String contig
        String prefix

        Int variants_per_shard
        Boolean create_variant_attributes = false

        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_override_preprocess_vcf
        RuntimeAttr? runtime_override_collect_sharded_vcf_stats
        RuntimeAttr? runtime_override_svtk_vcf_2_bed
        RuntimeAttr? runtime_override_scatter_vcf
        RuntimeAttr? runtime_override_merge_subvcf_stat_shards
        RuntimeAttr? runtime_override_merge_svtk_vcf_2_bed
    }

    String output_prefix = "~{prefix}.collect_qc_vcf_wide"

    call MiniTasks.ScatterVcf {
        input:
            vcf = vcf,
            vcf_index = vcf_idx,
            contig = contig,
            records_per_shard = variants_per_shard,
            prefix = "~{output_prefix}.scatter_vcf",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_override_scatter_vcf
    }
    Array[File] vcf_shards = ScatterVcf.shards

    scatter (i in range(length(vcf_shards))) {
        if (create_variant_attributes) {
            call AnnotateVariantAttributes {
                input:
                    vcf = vcf_shards[i],
                    prefix = "~{output_prefix}.annotate_attrs.shard_~{i}",
                    sv_pipeline_docker = sv_pipeline_docker
            }
        }

        call PreprocessVcf {
            input:
                vcf = select_first([AnnotateVariantAttributes.annotated_vcf, vcf_shards[i]]),
                prefix = "~{output_prefix}.preprocess.shard_~{i}",
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_override_preprocess_vcf
        }

        call CollectShardedVcfStats {
            input:
                vcf = PreprocessVcf.outvcf,
                prefix = "~{output_prefix}.collect_stats.shard_~{i}",
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_override_collect_sharded_vcf_stats
        }

        call SvtkVcf2bed {
            input:
                vcf = PreprocessVcf.outvcf,
                prefix = "~{output_prefix}.shard_~{i}",
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_override_svtk_vcf_2_bed
        }
    }

    call MiniTasks.ConcatBeds as MergeSubvcfStatShards {
        input:
            shard_bed_files = CollectShardedVcfStats.vcf_stats,
            prefix = "~{output_prefix}.VCF_sites.stats",
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_override_merge_subvcf_stat_shards
    }

    call MiniTasks.ConcatBeds as MergeSvtkVcf2bed {
        input:
            shard_bed_files = SvtkVcf2bed.vcf2bed_subworkflow_out,
            prefix = "~{output_prefix}.vcf2bed_subworkflow",
            index_output = false,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_override_merge_svtk_vcf_2_bed
    }

    output {
        File vcf_stats = MergeSubvcfStatShards.merged_bed_file
        File vcf_stats_idx = MergeSubvcfStatShards.merged_bed_idx
        File samples_list = CollectShardedVcfStats.samples_list[0]
        File vcf2bed_out = MergeSvtkVcf2bed.merged_bed_file
    }
}

task PreprocessVcf {
    input {
        File vcf
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat << 'EOF' > convert_to_symbolic.py
import pysam
import sys

vcf_in = pysam.VariantFile("~{vcf}", "r")
if "SVTYPE" not in vcf_in.header.info:
    vcf_in.header.info.add("SVTYPE", 1, "String", "Type of variant.")
if "SVLEN" not in vcf_in.header.info:
    vcf_in.header.info.add("SVLEN", 1, "Integer", "Length of variant.")
if "END" not in vcf_in.header.info:
    vcf_in.header.info.add("END", 1, "Integer", "End position of variant.")
if "ORIG_ALT" not in vcf_in.header.info:
    vcf_in.header.info.add("ORIG_ALT", 1, "String", "Original ALT nucleotide for SNVs before symbolic conversion.")
if "TRV_EXPANSION_RATIO" not in vcf_in.header.info:
    vcf_in.header.info.add("TRV_EXPANSION_RATIO", 1, "Float", "Proportion of non-ref alleles called as expansions (alt allele length > ref allele length) among non-neutral alleles for TRV variants.")
vcf_out = pysam.VariantFile("~{prefix}.vcf.gz", "w", header=vcf_in.header)

for rec in vcf_in:
    allele_type = rec.info["allele_type"].upper()
    allele_length = abs(rec.info["allele_length"]) if "allele_length" in rec.info else len(rec.ref)

    if "DEL" in allele_type:
        svtype = "DEL_SHORT" if allele_length < 50 else "DEL_SV"
    elif "INS" in allele_type:
        svtype = "INS_SHORT" if allele_length < 50 else "INS_SV"
    elif "DUP" in allele_type:
        svtype = "DUP_SHORT" if allele_length < 50 else "DUP_SV"
    elif allele_type == "TRV":
        motifs = rec.info.get("MOTIFS", None)
        if motifs:
            motif_list = list(motifs) if isinstance(motifs, (list, tuple)) else [str(motifs)]
            min_motif_len = min(len(str(m)) for m in motif_list)
        else:
            min_motif_len = 0
        svtype = "TR" if min_motif_len <= 6 else "VNTR"
    else:
        svtype = allele_type
    
    if allele_type == "SNV" and len(rec.alts) == 1 and len(rec.alts[0]) == 1:
        rec.info["ORIG_ALT"] = rec.alts[0]
    rec.info["SVTYPE"] = svtype
    rec.info["SVLEN"] = allele_length
    rec.stop = rec.pos + allele_length - 1

    # For TRV: compute expansion ratio from AL FORMAT field before GT collapsing
    if allele_type.lower() == "trv":
        expansions = 0
        contractions = 0
        for sample in rec.samples.values():
            gt = sample.get("GT", (None,))
            al_raw = sample.get("AL", None)
            if al_raw is None or gt is None:
                continue
            al = list(al_raw)
            try:
                ref_len = al[0]
            except (IndexError, TypeError):
                continue
            if ref_len is None:
                continue
            for allele_idx in gt:
                if allele_idx is None or allele_idx == 0:
                    continue
                try:
                    alt_len = al[allele_idx]
                except (IndexError, TypeError):
                    continue
                if alt_len is None:
                    continue
                if alt_len > ref_len:
                    expansions += 1
                elif alt_len < ref_len:
                    contractions += 1
        total = expansions + contractions
        if total > 0:
            rec.info["TRV_EXPANSION_RATIO"] = float(expansions) / total

    if len(rec.alts) > 1:
        for sample in rec.samples.values():
            gt = sample["GT"]
            sample["GT"] = tuple(0 if a == 0 else (1 if a is not None else None) for a in gt)
    
    rec.alts = (f"<{allele_type}>",)

    # Strip phasing from all sample genotypes so 1|0 and 0|1 are treated as 0/1
    for s_name in rec.samples:
        try:
            rec.samples[s_name].phased = False
        except (AttributeError, TypeError):
            pass

    vcf_out.write(rec)

vcf_out.close()
vcf_in.close()
EOF

        python convert_to_symbolic.py

        tabix ~{prefix}.vcf.gz
    >>>

    output {
        File outvcf = "~{prefix}.vcf.gz"
        File outvcf_index = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr runtime_default = object {
        cpu_cores: 1,
        mem_gb: 12,
        disk_gb: 4 * ceil(size(vcf, "GiB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: 12 + " GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}

task CollectShardedVcfStats {
    input {
        File vcf
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf, "GiB")

    command <<<
        set -euo pipefail
    
        /opt/sv-pipeline/scripts/vcf_qc/collectQC.vcf_wide.sh \
            ~{vcf} \
            /opt/sv-pipeline/scripts/vcf_qc/SV_colors.txt \
            collectQC_vcfwide_output/
        
        cp collectQC_vcfwide_output/data/VCF_sites.stats.bed.gz ~{prefix}.VCF_sites.stats.bed.gz
        cp collectQC_vcfwide_output/data/VCF_sites.stats.bed.gz.tbi ~{prefix}.VCF_sites.stats.bed.gz.tbi
        cp collectQC_vcfwide_output/analysis_samples.list ~{prefix}.analysis_samples.list
        
        tar -czvf \
            ~{prefix}.collectQC_vcfwide_output.tar.gz \
            collectQC_vcfwide_output
    >>>

    output {
        File vcf_stats = "~{prefix}.VCF_sites.stats.bed.gz"
        File vcf_stats_idx = "~{prefix}.VCF_sites.stats.bed.gz.tbi"
        File samples_list = "~{prefix}.analysis_samples.list"
        File vcfwide_tarball = "~{prefix}.collectQC_vcfwide_output.tar.gz"
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 1.5 + 2.0 * input_size,
        disk_gb: ceil(10.0 + 2.0 * input_size),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 0,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}

task SvtkVcf2bed {
    input {
        File vcf
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    String output_file = "~{prefix}.vcf2bed_subworkflow.bed.gz"
    Float input_size = size(vcf, "GiB")

    command <<<
        set -euo pipefail
        
        svtk vcf2bed --info ALL ~{vcf} stdout \
            | bgzip -c \
            > "~{output_file}"
    >>>

    output {
        File vcf2bed_subworkflow_out = output_file
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(10.0 + input_size * 2.0),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 0,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}

task AnnotateVariantAttributes {
    input {
        File vcf
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf, "GiB")

    command <<<
        set -euo pipefail

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' \
            ~{vcf} \
        | awk -F'\t' '{
            ref_len = length($3)
            alt_len = length($4)
            diff = alt_len - ref_len
            if (ref_len == 1 && alt_len == 1) {
                atype = "snv"
            } else if (alt_len > ref_len) {
                atype = "ins"
            } else {
                atype = "del"
            }
            print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"diff"\t"atype
        }' | bgzip -c > annot.txt.gz

        tabix -s1 -b2 -e2 annot.txt.gz

        printf '##INFO=<ID=allele_length,Number=1,Type=Integer,Description="Allele length">\n##INFO=<ID=allele_type,Number=1,Type=String,Description="Allele type">\n' > new_headers.txt

        bcftools annotate \
            -h new_headers.txt \
            -a annot.txt.gz \
            -c CHROM,POS,REF,ALT,~ID,INFO/allele_length,INFO/allele_type \
            -Oz -o ~{prefix}.vcf.gz \
            ~{vcf}

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(10.0 + input_size * 3.0),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 0,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}
