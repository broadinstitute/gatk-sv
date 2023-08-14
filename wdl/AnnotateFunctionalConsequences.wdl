version 1.0

import "Structs.wdl"

workflow AnnotateFunctionalConsequences {
  input {
    File vcf
    File? vcf_index
    String prefix

    File protein_coding_gtf
    File? noncoding_bed
    Int? promoter_window
    Int? max_breakend_as_cnv_length
    String? additional_args

    String gatk_docker
    RuntimeAttr? runtime_attr_svannotate
  }

  call SVAnnotate {
    input:
      vcf = vcf,
      vcf_index = vcf_index,
      prefix = prefix,
      protein_coding_gtf = protein_coding_gtf,
      noncoding_bed = noncoding_bed,
      promoter_window = promoter_window,
      max_breakend_as_cnv_length = max_breakend_as_cnv_length,
      additional_args = additional_args,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_svannotate
  }

  output {
    File annotated_vcf = SVAnnotate.annotated_vcf
    File annotated_vcf_index = SVAnnotate.annotated_vcf_index
  }
}

task SVAnnotate {
  input {
    File vcf
    File? vcf_index
    String prefix

    File protein_coding_gtf
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

    # check index is in expected location. if not, tabix
    if [ ! -f "~{vcf}.tbi" ]; then
      tabix -p vcf ~{vcf}
    fi

    gatk --java-options "-Xmx~{java_mem_mb}m" SVAnnotate \
      -V ~{vcf} \
      -O ~{outfile} \
      --protein-coding-gtf ~{protein_coding_gtf} \
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