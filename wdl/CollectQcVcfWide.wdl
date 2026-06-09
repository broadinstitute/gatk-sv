version 1.0

import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow CollectQcVcfWide {
    input {
        File vcf
        File? vcf_idx
        String contig
        String prefix

        Int variants_per_shard

        String? subset_vcf_string
        Boolean create_variant_attributes = false

        File ref_fa
        File ref_fai

        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_override_scatter_vcf
        RuntimeAttr? runtime_override_preprocess_vcf
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
        call PreprocessCollectAndConvert {
            input:
                vcf = vcf_shards[i],
                subset_flags = subset_vcf_string,
                prefix = "~{output_prefix}.shard_~{i}",
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                create_variant_attributes = create_variant_attributes,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_override_preprocess_vcf
        }
    }

    call MiniTasks.ConcatBeds as MergeSubvcfStatShards {
        input:
            shard_bed_files = PreprocessCollectAndConvert.vcf_stats,
            prefix = "~{output_prefix}.VCF_sites.stats",
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_override_merge_subvcf_stat_shards
    }

    call MiniTasks.ConcatBeds as MergeSvtkVcf2bed {
        input:
            shard_bed_files = PreprocessCollectAndConvert.vcf2bed_out,
            prefix = "~{output_prefix}.vcf2bed_subworkflow",
            index_output = false,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_override_merge_svtk_vcf_2_bed
    }

    output {
        File vcf_stats = MergeSubvcfStatShards.merged_bed_file
        File vcf_stats_idx = MergeSubvcfStatShards.merged_bed_idx
        File samples_list = PreprocessCollectAndConvert.samples_list[0]
        File vcf2bed_out = MergeSvtkVcf2bed.merged_bed_file
    }
}

task PreprocessCollectAndConvert {
    input {
        File vcf
        String? subset_flags
        String prefix
        File ref_fa
        File ref_fai
        Boolean create_variant_attributes = false
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Optionally pre-filter records on the shard (replaces the standalone SubsetVcf task)
        if ~{if defined(subset_flags) then "true" else "false"}; then
            bcftools view --no-version ~{subset_flags} -Oz -o filtered_input.vcf.gz ~{vcf}
        else
            ln -s ~{vcf} filtered_input.vcf.gz
        fi

        # Normalize multiallelics
        echo "Normalizing multiallelics"
        bcftools norm \
            -m-any \
            --check-ref w \
            --fasta-ref ~{ref_fa} \
            -Oz -o ~{prefix}.split.vcf.gz \
            filtered_input.vcf.gz

        rm -f filtered_input.vcf.gz

        # Convert to symbolic, optionally annotating allele_length / allele_type
        echo "Converting to symbolic"
        python <<'PYCODE'
import pysam

VCF_IN = "~{prefix}.split.vcf.gz"
VCF_OUT = "~{prefix}.unsorted.vcf.gz"
CREATE_VARIANT_ATTRIBUTES = "~{create_variant_attributes}".lower() == "true"

group_carrier_counts = []
vcf_in = pysam.VariantFile(VCF_IN, "r")
current_locus = None
current_carrier_set = set()
for rec in vcf_in:
    if CREATE_VARIANT_ATTRIBUTES:
        ref_len = len(rec.ref)
        first_alt = rec.alts[0] if rec.alts else ""
        alt_len = len(first_alt)
        if ref_len == 1 and alt_len == 1:
            at = "SNV"
        elif alt_len > ref_len:
            at = "INS"
        else:
            at = "DEL"
    else:
        at = rec.info["allele_type"].upper() if "allele_type" in rec.info else ""

    if at == "TRV":
        key = (rec.chrom, rec.id)
        if key != current_locus:
            if current_locus is not None:
                group_carrier_counts.append(len(current_carrier_set))
            current_locus = key
            current_carrier_set = set()
        for s_name in rec.samples:
            gt = rec.samples[s_name].get("GT", (None,))
            if gt and any(a is not None and a > 0 for a in gt):
                current_carrier_set.add(s_name)
    elif current_locus is not None:
        group_carrier_counts.append(len(current_carrier_set))
        current_locus = None
        current_carrier_set = set()

if current_locus is not None:
    group_carrier_counts.append(len(current_carrier_set))
vcf_in.close()

vcf_in = pysam.VariantFile(VCF_IN, "r")

if CREATE_VARIANT_ATTRIBUTES:
    if "allele_length" not in vcf_in.header.info:
        vcf_in.header.info.add("allele_length", 1, "Integer", "Allele length (alt_len - ref_len)")
    if "allele_type" not in vcf_in.header.info:
        vcf_in.header.info.add("allele_type", 1, "String", "Allele type (snv/ins/del)")

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
vcf_out = pysam.VariantFile(VCF_OUT, "w", header=vcf_in.header)

prev_trv_key = None
trv_counter = 0

prev_nontrv_key = None
nontrv_counter = 1

group_idx = -1
current_group_carriers = 0

for rec in vcf_in:
    if CREATE_VARIANT_ATTRIBUTES:
        ref_len = len(rec.ref)
        first_alt = rec.alts[0] if rec.alts else ""
        alt_len = len(first_alt)
        if ref_len == 1 and alt_len == 1:
            inferred_type = "snv"
        elif alt_len > ref_len:
            inferred_type = "ins"
        else:
            inferred_type = "del"
        rec.info["allele_length"] = alt_len - ref_len
        rec.info["allele_type"] = inferred_type

    allele_type = rec.info["allele_type"].upper()

    if allele_type == "TRV":
        current_key = (rec.chrom, rec.id)
        if current_key != prev_trv_key:
            group_idx += 1
            current_group_carriers = group_carrier_counts[group_idx]
            trv_counter = 1
            prev_trv_key = current_key
        else:
            trv_counter += 1
        rec.id = f"{rec.id}_{trv_counter}"

        motifs = rec.info.get("MOTIFS", None)
        if motifs:
            motif_list = list(motifs) if isinstance(motifs, (list, tuple)) else [str(motifs)]
            min_motif_len = min(len(str(m)) for m in motif_list)
        else:
            min_motif_len = 0
        svtype = "TR" if min_motif_len <= 6 else "VNTR"

        ref_seq = rec.ref
        alt_seq = rec.alts[0] if rec.alts else ""
        if len(alt_seq) > len(ref_seq):
            allele_class = "TR_INS"
            allele_length = len(alt_seq) - len(ref_seq)
        elif len(ref_seq) > len(alt_seq):
            allele_class = "TR_DEL"
            allele_length = len(ref_seq) - len(alt_seq)
        else:
            allele_class = "TR_SNV"
            allele_length = max(len(ref_seq), 1)

        rec.info["SVTYPE"] = svtype
        rec.info["SVLEN"] = allele_length
        rec.info["TR_ALLELE_CLASS"] = allele_class
        rec.info["TRV_LOCUS_CARRIERS"] = current_group_carriers
        rec.stop = rec.pos + len(rec.ref)

        rec.alts = (f"<{allele_type}>",)

        for s_name in rec.samples:
            try:
                rec.samples[s_name].phased = False
            except (AttributeError, TypeError):
                pass

        vcf_out.write(rec)
        continue

    prev_trv_key = None

    nontrv_key = (rec.chrom, rec.pos, rec.id)
    if nontrv_key == prev_nontrv_key:
        nontrv_counter += 1
        rec.id = f"{rec.id}_{nontrv_counter}"
    else:
        prev_nontrv_key = nontrv_key
        nontrv_counter = 1

    allele_length = abs(rec.info["allele_length"]) if "allele_length" in rec.info else len(rec.ref)

    if "DEL" in allele_type:
        svtype = "DEL_SHORT" if allele_length < 50 else "DEL_SV"
    elif "INS" in allele_type:
        svtype = "INS_SHORT" if allele_length < 50 else "INS_SV"
    elif "DUP" in allele_type:
        svtype = "DUP_SHORT" if allele_length < 50 else "DUP_SV"
    else:
        svtype = allele_type

    if allele_type == "SNV" and len(rec.alts[0]) == 1:
        rec.info["ORIG_ALT"] = rec.alts[0]
    rec.info["SVTYPE"] = svtype
    rec.info["SVLEN"] = allele_length
    rec.stop = rec.pos + len(rec.ref)

    rec.alts = (f"<{allele_type}>",)

    for s_name in rec.samples:
        try:
            rec.samples[s_name].phased = False
        except (AttributeError, TypeError):
            pass

    vcf_out.write(rec)

vcf_out.close()
vcf_in.close()
PYCODE

        rm -f ~{prefix}.split.vcf.gz ~{prefix}.split.vcf.gz.tbi

        # Sort
        echo "Sorting"

        mkdir -p /tmp/bcftools_sort/

        bcftools sort \
            --max-mem ~{sort_max_mem_gb}G \
            -T /tmp/bcftools_sort/ \
            -Oz -o ~{prefix}.sorted.vcf.gz \
            ~{prefix}.unsorted.vcf.gz

        rm -f ~{prefix}.unsorted.vcf.gz

        tabix ~{prefix}.sorted.vcf.gz

        # Collect stats
        echo "Collecting stats"

        /opt/sv-pipeline/scripts/vcf_qc/collectQC.vcf_wide.sh \
            ~{prefix}.sorted.vcf.gz \
            /opt/sv-pipeline/scripts/vcf_qc/SV_colors.txt \
            collectQC_vcfwide_output/
        
        cp collectQC_vcfwide_output/data/VCF_sites.stats.bed.gz ~{prefix}.VCF_sites.stats.bed.gz
        cp collectQC_vcfwide_output/data/VCF_sites.stats.bed.gz.tbi ~{prefix}.VCF_sites.stats.bed.gz.tbi
        cp collectQC_vcfwide_output/analysis_samples.list ~{prefix}.analysis_samples.list
        bgzip -c collectQC_vcfwide_output/vcf2bed_raw.bed > ~{prefix}.vcf2bed.bed.gz

        rm -f ~{prefix}.sorted.vcf.gz ~{prefix}.sorted.vcf.gz.tbi
        rm -rf collectQC_vcfwide_output/
    >>>

    output {
        File vcf_stats = "~{prefix}.VCF_sites.stats.bed.gz"
        File vcf_stats_idx = "~{prefix}.VCF_sites.stats.bed.gz.tbi"
        File samples_list = "~{prefix}.analysis_samples.list"
        File vcf2bed_out = "~{prefix}.vcf2bed.bed.gz"
    }

    RuntimeAttr runtime_default = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 5 * ceil(size(vcf, "GiB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    Int sort_max_mem_gb = floor(0.9 * select_first([runtime_override.mem_gb, runtime_default.mem_gb]))
    runtime {
        memory: 100 + " GiB"
        disks: "local-disk " + 50 + " HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}
