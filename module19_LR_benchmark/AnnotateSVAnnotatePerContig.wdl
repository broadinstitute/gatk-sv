version 1.0

import "Structs.wdl"

workflow AnnotateSVAnnotatePerContig {
    input {
        File vcf
        File vcf_idx
        String contig
        String prefix
        Boolean reheader = false
        File? new_header

        Int min_length
        String? vcf_drop_fields
        
        File coding_gtf
        File noncoding_bed

        String utils_docker
        String gatk_docker
        String sv_base_mini_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_drop_fields
        RuntimeAttr? runtime_attr_convert_symbolic
        RuntimeAttr? runtime_attr_annotate_func
        RuntimeAttr? runtime_attr_concat_unannotated
        RuntimeAttr? runtime_attr_concat_annotated
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_revert_symbolic
    }

    if (reheader){
        call ReheaderVcf {
            input:
                input_vcf = ~{vcf},
                new_header = new_header,
                docker_image = sv_base_mini_docker
        }
    }

    if (defined(vcf_drop_fields)) {
        call Helpers.DropVcfFields {
            input:
                vcf = select_first([ReheaderVcf.output_vcf, vcf]),
                vcf_idx = select_first([ReheaderVcf.output_vcf_index, vcf_idx]),
                drop_fields = select_first([vcf_drop_fields]),
                prefix = "~{prefix}.~{contig}.dropped",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_drop_fields
        }
    }

    call Helpers.ConvertToSymbolic {
        input:
            vcf = select_first([DropVcfFields.dropped_vcf, SubsetVcfAnnotated.subset_vcf]),
            vcf_idx = select_first([DropVcfFields.dropped_vcf_idx, SubsetVcfAnnotated.subset_vcf_idx]),
            prefix = "~{prefix}.~{contig}.converted",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_convert_symbolic
    }

    call AnnotateFunctionalConsequences {
            input:
                vcf = ConvertToSymbolic.processed_vcf,
                vcf_idx = ConvertToSymbolic.processed_vcf_idx,
                noncoding_bed = noncoding_bed,
                coding_gtf = coding_gtf,
                prefix = "~{prefix}.~{contig}.functionally_annotated",
                docker = gatk_docker,
                runtime_attr_override = runtime_attr_annotate_func
    }

    call Helpers.RevertSymbolicAlleles {
        input:
            annotated_vcf = AnnotateFunctionalConsequences.anno_vcf,
            annotated_vcf_idx = AnnotateFunctionalConsequences.anno_vcf_idx,
            original_vcf = SubsetVcfAnnotated.subset_vcf,
            original_vcf_idx = SubsetVcfAnnotated.subset_vcf_idx,
            prefix = "~{prefix}.~{contig}.reverted",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_revert_symbolic
    }

    call Helpers.ExtractVcfAnnotations as ExtractAnnotations {
        input:
            vcf = RevertSymbolicAlleles.reverted_vcf,
            vcf_idx = RevertSymbolicAlleles.reverted_vcf_idx,
            original_vcf = SubsetVcfAnnotated.subset_vcf,
            original_vcf_idx = SubsetVcfAnnotated.subset_vcf_idx,
            prefix = "~{prefix}.~{contig}",
            docker = utils_docker
    }


    call Helpers.ConcatTsvs as MergeTsvs {
        input:
            tsvs = ExtractAnnotations.annotations_tsv,
            prefix = prefix + ".svannotate_annotations",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_annotated
    }

    call Helpers.MergeHeaderLines as MergeHeaders {
        input:
            header_files = ExtractAnnotations.annotations_header,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge
    }

    output {
        File annotations_tsv_svannotate = MergeTsvs.concatenated_tsv
        File annotations_header_svannotate = MergeHeaders.merged_header
    }
}

task AnnotateFunctionalConsequences {
    input {
        File vcf
        File vcf_idx
        File noncoding_bed
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
            --non-coding-bed ~{noncoding_bed} \
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


task ReheaderVcf {
  input {
    File input_vcf
    File? new_header
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: ceil(10 + size(input_vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  String prefix = basename(input_vcf, ".vcf.gz")

  command <<<
    set -euo pipefail

    bcftools reheader \
      -h ~{new_header} \
      -o ~{prefix}.reheadered.vcf.gz \
      ~{input_vcf}

    bcftools index ~{prefix}.reheadered.vcf.gz
  >>>

  output {
    File output_vcf = "~{prefix}.reheadered.vcf.gz"
    File output_vcf_index = "~{prefix}.reheadered.vcf.gz.csi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

