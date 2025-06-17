version 1.0

import "Structs.wdl"

workflow AnnotateFunctionalConsequences {
  input {
    File vcf
    File? vcf_index
    String prefix

    File? protein_coding_gtf
    File? noncoding_bed
    Int? promoter_window
    Int? max_breakend_as_cnv_length
    String? additional_args
    Int? min_annotation_size

    String gatk_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_svannotate
    RuntimeAttr? runtime_attr_subset_vcf
    RuntimeAttr? runtime_attr_concat_vcf
  }

  if (defined(min_annotation_size)) {
    call SubsetVcfBySvlen {
      input:
        vcf = vcf,
        min_svlen = select_first([min_annotation_size]),
        prefix = prefix,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_subset_vcf
    }
  }

  call SVAnnotate {
    input:
      vcf = select_first([SubsetVcfBySvlen.subset_vcf, vcf]),
      vcf_index = SubsetVcfBySvlen.subset_vcf_index,
      prefix = prefix,
      protein_coding_gtf = protein_coding_gtf,
      noncoding_bed = noncoding_bed,
      promoter_window = promoter_window,
      max_breakend_as_cnv_length = max_breakend_as_cnv_length,
      additional_args = additional_args,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_svannotate
  }

  if (defined(min_annotation_size)) {
    call ConcatAnnotatedVcf {
      input:
        original_vcf = vcf,
        annotated_vcf = SVAnnotate.annotated_vcf,
        min_svlen = select_first([min_annotation_size]),
        prefix = prefix,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_concat_vcf
    }
  }

  output {
    File annotated_vcf = select_first([ConcatAnnotatedVcf.final_vcf, SVAnnotate.annotated_vcf])
    File annotated_vcf_index = select_first([ConcatAnnotatedVcf.final_vcf_index, SVAnnotate.annotated_vcf_index])
  }
}

task SubsetVcfBySvlen {
  input {
    File vcf
    Int min_svlen
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: ceil(10 + size(vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String subset_vcf_name = "~{prefix}.subset.vcf.gz"

  command <<<
    set -euo pipefail

    if [ ! -f "~{vcf}.tbi" ]; then
      tabix -p vcf ~{vcf}
    fi

    bcftools view -i "INFO/SVLEN >= ~{min_svlen}" -O z -o ~{subset_vcf_name} ~{vcf}
    tabix -p vcf ~{subset_vcf_name}
  >>>

  output {
    File subset_vcf = subset_vcf_name
    File subset_vcf_index = "~{subset_vcf_name}.tbi"
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

task SVAnnotate {
  input {
    File vcf
    File? vcf_index
    String prefix

    File? protein_coding_gtf
    File? noncoding_bed
    Int? promoter_window
    Int? max_breakend_as_cnv_length
    String? additional_args

    String gatk_docker
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

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.7)

  String outfile = "~{prefix}.annotated.vcf.gz"

  output {
    File annotated_vcf = "~{outfile}"
    File annotated_vcf_index = "~{outfile}.tbi"
  }
  command <<<
    set -euo pipefail

    if [ ! -f "~{vcf}.tbi" ]; then
      tabix -p vcf ~{vcf}
    fi

    gatk --java-options "-Xmx~{java_mem_mb}m" SVAnnotate \
      -V ~{vcf} \
      -O ~{outfile} \
      ~{"--protein-coding-gtf " + protein_coding_gtf} \
      ~{"--non-coding-bed " + noncoding_bed} \
      ~{"--promoter-window-length " + promoter_window} \
      ~{"--max-breakend-as-cnv-length " + max_breakend_as_cnv_length} \
      ~{additional_args}

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ConcatAnnotatedVcf {
  input {
    File original_vcf
    File annotated_vcf
    Int min_svlen
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: ceil(10 + size([original_vcf, annotated_vcf], "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String final_vcf_name = "~{prefix}.final_annotated.vcf.gz"

  command <<<
    set -euo pipefail

    if [ ! -f "~{original_vcf}.tbi" ]; then
      tabix -p vcf ~{original_vcf}
    fi

    bcftools view -e "INFO/SVLEN >= ~{min_svlen}" -O z -o unannotated.vcf.gz ~{original_vcf}
    tabix -p vcf unannotated.vcf.gz

    bcftools concat -a -O z -o ~{final_vcf_name} ~{annotated_vcf} unannotated.vcf.gz
    tabix -p vcf ~{final_vcf_name}
  >>>

  output {
    File final_vcf = final_vcf_name
    File final_vcf_index = "~{final_vcf_name}.tbi"
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