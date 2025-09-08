version 1.0

import "Structs.wdl"

workflow FederateCohorts {
  input {
    File vcf_a
    File vcf_b

    String prefix_a
    String prefix_b
    String additional_args

    File scores
    String prefix

    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    Float? java_mem_fraction
    String gatk_docker
  }

  call SVFederate {
    input:
      vcf_a=vcf_a,
      vcf_b=vcf_b,
      scores=scores,
      prefix=prefix,
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_fasta_fai,
      reference_dict=reference_dict,
      java_mem_fraction=java_mem_fraction,
      gatk_docker=gatk_docker
  }

  output {
    File federated_vcf = SVFederate.out
    File federated_vcf_index = SVFederate.out_index
  }
}

task SVFederate {
  input {
    File vcf_a
    File vcf_b

    String prefix_a
    String prefix_b

    String prefix
    File scores
    String additional_args

    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    Float? java_mem_fraction
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    vcf_a: {
                 localization_optional: true
               }
    vcf_b:  {
                 localization_optional: true
               }
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 16,
                               disk_gb: ceil(10 + size(vcf_a, "GB") * 4 + size(vcf_b, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{prefix}.vcf.gz"
    File out_index = "~{prefix}.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

    function getJavaMem() {
    # get JVM memory in MiB by getting total memory from /proc/meminfo
    # and multiplying by java_mem_fraction
      cat /proc/meminfo \
        | awk -v MEM_FIELD="$1" '{
          f[substr($1, 1, length($1)-1)] = $2
        } END {
          printf "%dM", f[MEM_FIELD] * ~{default="0.85" java_mem_fraction} / 1024
        }'
    }
    JVM_MAX_MEM=$(getJavaMem MemTotal)
    echo "JVM memory: $JVM_MAX_MEM"

    gatk --java-options "-Xmx${JVM_MAX_MEM}" SVFederate \
      -A ~{vcf_a} \
      -B ~{vcf_b} \
      --prefix-A ~{prefix_a} \
      --prefix-B ~{prefix_b} \
      --sv-pairs ~{scores} \
      -R ~{reference_fasta} \
      ~{additional_args} \
      -O ~{prefix}.vcf.gz

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
