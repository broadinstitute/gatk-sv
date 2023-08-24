## This WDL pipeline implements mobile element calling with DeepMEI

version 1.0

import "Structs.wdl"

workflow DeepMei {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_id
    String? reference_id
    String deepmei_docker
    RuntimeAttr? runtime_attr_override
  }

  call RunDeepMei {
    input:
      bam_or_cram_file = bam_or_cram_file,
      bam_or_cram_index = bam_or_cram_index,
      sample_id = sample_id,
      reference_id = reference_id,
      deepmei_docker = deepmei_docker,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File deepmei_vcf = RunDeepMei.vcf
    File deepmei_vcf_index = RunDeepMei.vcf_index
  }
}

task RunDeepMei {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_id
    String reference_id = "38"
    String deepmei_docker
    RuntimeAttr? runtime_attr_override
  }

  Float disk_overhead = 100.0
  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + disk_overhead)

  RuntimeAttr default_attr = object {
                               cpu_cores: 32,
                               mem_gb: 15.0,
                               disk_gb: vm_disk_size,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File vcf = "~{sample_id}.deepmei.vcf.gz"
    File vcf_index = "~{sample_id}.deepmei.vcf.gz.tbi"
  }
  command <<<
    set -euxo pipefail

    export PATH=/root/miniconda3/bin:$PATH
    CRAM_NAME=$(basename ~{bam_or_cram_file})
    OUTDIR=${PWD}/output/
    mkdir ${OUTDIR}

    bash /root/DeepMEI/DeepMEI_model/model_test_batch.sh \
      -i ~{bam_or_cram_file} \
      -r ~{reference_id} \
      -w ${PWD}/output/

    mv ${OUTDIR}/DeepMEI_output/${CRAM_NAME}/${CRAM_NAME}.vcf ~{sample_id}.deepmei.vcf
    bgzip ~{sample_id}.deep_mei.vcf
    tabix ~{sample_id}.deep_mei.vcf.gz
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: deepmei_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

