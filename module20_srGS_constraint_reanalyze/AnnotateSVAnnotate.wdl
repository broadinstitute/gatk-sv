version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateSVAnnotate {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int min_length
        
        File coding_gtf

        String utils_docker
        String gatk_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_convert_symbolic
        RuntimeAttr? runtime_attr_annotate_func
        RuntimeAttr? runtime_attr_concat_unannotated
        RuntimeAttr? runtime_attr_concat_annotated
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_revert_symbolic
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        call Helpers.SubsetVcfByLength {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                min_length = min_length,
                extra_args = if single_contig then "-G" else "-G --regions ~{contig}",
                prefix = "~{prefix}.~{contig}.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call Helpers.ConvertToSymbolic {
            input:
                vcf = SubsetVcfByLength.subset_vcf,
                vcf_idx = SubsetVcfByLength.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.converted",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_convert_symbolic
        }

        call AnnotateFunctionalConsequences {
            input:
                vcf = ConvertToSymbolic.processed_vcf,
                vcf_idx = ConvertToSymbolic.processed_vcf_idx,
                coding_gtf = coding_gtf,
                prefix = "~{prefix}.~{contig}.functionally_annotated",
                docker = gatk_docker,
                runtime_attr_override = runtime_attr_annotate_func
        }

        call Helpers.RevertSymbolicAlleles {
            input:
                annotated_vcf = AnnotateFunctionalConsequences.anno_vcf,
                annotated_vcf_idx = AnnotateFunctionalConsequences.anno_vcf_idx,
                original_vcf = SubsetVcfByLength.subset_vcf,
                original_vcf_idx = SubsetVcfByLength.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.reverted",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_revert_symbolic
        }

        call Helpers.ExtractVcfAnnotations {
            input:
                vcf = RevertSymbolicAlleles.reverted_vcf,
                vcf_idx = RevertSymbolicAlleles.reverted_vcf_idx,
                original_vcf = SubsetVcfByLength.subset_vcf,
                original_vcf_idx = SubsetVcfByLength.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker
        }
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs {
            input:
                tsvs = ExtractVcfAnnotations.annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.svannotate_annotations",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_annotated
        }

        call Helpers.MergeHeaderLines {
            input:
                header_files = ExtractVcfAnnotations.annotations_header,
                prefix = "~{prefix}.svannotate_annotations",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }
    }

    output {
        File annotations_tsv_svannotate = select_first([ConcatTsvs.concatenated_tsv, ExtractVcfAnnotations.annotations_tsv[0]])
        File annotations_header_svannotate = select_first([MergeHeaderLines.merged_header, ExtractVcfAnnotations.annotations_header[0]])
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
