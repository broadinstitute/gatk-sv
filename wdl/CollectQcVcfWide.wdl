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

        # Split multiallelic TRV records into biallelic before processing
        bcftools norm -m-any ~{vcf} -Oz -o ~{prefix}.split.vcf.gz
        tabix -p vcf ~{prefix}.split.vcf.gz

        cat << 'EOF' > convert_to_symbolic.py
import pysam
import sys

vcf_in = pysam.VariantFile("~{prefix}.split.vcf.gz", "r")
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
if "TR_ALLELE_CLASS" not in vcf_in.header.info:
    vcf_in.header.info.add("TR_ALLELE_CLASS", 1, "String", "Sub-classification of TRV biallelic allele: TR_INS, TR_DEL, or TR_SNV.")
if "TRV_LOCUS_CARRIERS" not in vcf_in.header.info:
    vcf_in.header.info.add("TRV_LOCUS_CARRIERS", 1, "Integer", "Number of unique carriers across all biallelic alleles at this TRV locus.")
vcf_out = pysam.VariantFile("~{prefix}.vcf.gz", "w", header=vcf_in.header)

# Track consecutive TRV records from the same original multiallelic
prev_trv_key = None
trv_counter = 0
trv_buffer = []

def trim_common_bases(ref, alt):
    """Trim common prefix and suffix. Returns (trimmed_ref, trimmed_alt, prefix_len)."""
    i = 0
    while i < len(ref) and i < len(alt) and ref[i] == alt[i]:
        i += 1
    prefix_len = i
    ref = ref[i:]
    alt = alt[i:]
    while len(ref) > 0 and len(alt) > 0 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    return ref, alt, prefix_len

def flush_trv_buffer():
    """Compute locus-level carriers and write buffered TRV records."""
    if not trv_buffer:
        return
    carrier_set = set()
    for rec in trv_buffer:
        for s_name in rec.samples:
            gt = rec.samples[s_name].get("GT", (None,))
            if gt and any(a is not None and a > 0 for a in gt):
                carrier_set.add(s_name)
    n_carriers = len(carrier_set)
    for rec in trv_buffer:
        rec.info["TRV_LOCUS_CARRIERS"] = n_carriers
        vcf_out.write(rec)
    trv_buffer.clear()

for rec in vcf_in:
    allele_type = rec.info["allele_type"].upper()

    if allele_type == "TRV":
        # Handle biallelic TRV (post bcftools norm -m-any)
        current_key = (rec.chrom, rec.pos, rec.id)
        if current_key != prev_trv_key:
            flush_trv_buffer()
            trv_counter = 1
            prev_trv_key = current_key
        else:
            trv_counter += 1

        # Add suffix to ID to track which original multiallelic this came from
        rec.id = f"{rec.id}_{trv_counter}"

        # Classify TR vs VNTR by min motif length
        motifs = rec.info.get("MOTIFS", None)
        if motifs:
            motif_list = list(motifs) if isinstance(motifs, (list, tuple)) else [str(motifs)]
            min_motif_len = min(len(str(m)) for m in motif_list)
        else:
            min_motif_len = 0
        svtype = "TR" if min_motif_len <= 6 else "VNTR"

        # Trim common prefix/suffix to get the true variant portion; update position accordingly
        ref_seq = rec.ref
        alt_seq = rec.alts[0] if rec.alts else ""
        trim_ref, trim_alt, prefix_len = trim_common_bases(ref_seq, alt_seq)
        rec.pos = rec.pos + prefix_len

        # Classify allele-level detail (TR_INS/TR_DEL/TR_SNV)
        if len(trim_alt) > len(trim_ref):
            allele_class = "TR_INS"
            allele_length = len(trim_alt) - len(trim_ref)
        elif len(trim_ref) > len(trim_alt):
            allele_class = "TR_DEL"
            allele_length = len(trim_ref) - len(trim_alt)
        else:
            allele_class = "TR_SNV"
            allele_length = max(len(trim_ref), 1)

        rec.info["SVTYPE"] = svtype
        rec.info["SVLEN"] = allele_length
        rec.info["TR_ALLELE_CLASS"] = allele_class
        rec.stop = rec.pos + allele_length - 1

        rec.alts = (f"<{allele_type}>",)

        for s_name in rec.samples:
            try:
                rec.samples[s_name].phased = False
            except (AttributeError, TypeError):
                pass

        trv_buffer.append(rec.copy())
        continue

    # Non-TRV processing
    flush_trv_buffer()
    prev_trv_key = None
    allele_length = abs(rec.info["allele_length"]) if "allele_length" in rec.info else len(rec.ref)

    if "DEL" in allele_type:
        svtype = "DEL_SHORT" if allele_length < 50 else "DEL_SV"
    elif "INS" in allele_type:
        svtype = "INS_SHORT" if allele_length < 50 else "INS_SV"
    elif "DUP" in allele_type:
        svtype = "DUP_SHORT" if allele_length < 50 else "DUP_SV"
    else:
        svtype = allele_type
    
    if allele_type == "SNV" and len(rec.alts) == 1 and len(rec.alts[0]) == 1:
        rec.info["ORIG_ALT"] = rec.alts[0]
    rec.info["SVTYPE"] = svtype
    rec.info["SVLEN"] = allele_length
    rec.stop = rec.pos + allele_length - 1

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

flush_trv_buffer()
vcf_out.close()
vcf_in.close()
EOF

        python convert_to_symbolic.py
        rm -f ~{prefix}.split.vcf.gz ~{prefix}.split.vcf.gz.tbi

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
