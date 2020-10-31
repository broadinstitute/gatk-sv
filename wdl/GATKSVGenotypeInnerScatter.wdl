version 1.0

import "Structs.wdl"

workflow GATKSVGenotypeInnerScatter {
  input {
    Array[File] vcfs
    File sample_mean_depth_file

    Int records_per_shard = 2000
    Int discrete_samples = 1000

    # Model parameters
    Float? eps_pe
    Float? eps_sr1
    Float? eps_sr2

    String genotyping_gatk_docker

    # cpu or cuda
    String device_train = "cpu"
    String device_genotype = "cpu"

    String gpu_type = "nvidia-tesla-k80"
    String nvidia_driver_version = "450.80.02"

    RuntimeAttr? runtime_attr_train
    RuntimeAttr? runtime_attr_infer
  }

  scatter (i in range(length(vcfs))) {
    String model_name = basename(vcfs[i], ".vcf.gz")
    call SVTrainGenotyping {
      input:
        vcf = vcfs[i],
        sample_mean_depth_file = sample_mean_depth_file,
        model_name = model_name,
        eps_pe = eps_pe,
        eps_sr1 = eps_sr1,
        eps_sr2 = eps_sr2,
        gatk_docker = genotyping_gatk_docker,
        device = device_train,
        gpu_type = gpu_type,
        nvidia_driver_version = nvidia_driver_version,
        runtime_attr_override = runtime_attr_train
    }
    call SVGenotype {
      input:
        vcf = vcfs[i],
        model_tar = SVTrainGenotyping.out,
        model_name = model_name,
        discrete_samples = discrete_samples,
        output_vcf_filename = "~{model_name}.genotyped.vcf.gz",
        gatk_docker = genotyping_gatk_docker,
        device = device_genotype,
        gpu_type = gpu_type,
        nvidia_driver_version = nvidia_driver_version,
        runtime_attr_override = runtime_attr_infer
    }
  }

  output {
    Array[File] out = SVGenotype.out
    Array[File] out_index = SVGenotype.out_index
  }
}

task SVTrainGenotyping {
  input {
    File vcf
    File sample_mean_depth_file
    String model_name
    String gatk_docker
    String device
    Float? eps_pe
    Float? eps_sr1
    Float? eps_sr2
    Int? max_iter
    String gpu_type
    String nvidia_driver_version
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2.0,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{model_name}.sv_genotype_model.tar.gz"
  }
  command <<<
    set -euo pipefail
    mkdir svmodel
    tabix ~{vcf}

    gatk --java-options -Xmx~{java_mem_mb}M SVTrainGenotyping \
      --variant ~{vcf} \
      --coverage-file ~{sample_mean_depth_file} \
      --output-name ~{model_name} \
      --output-dir svmodel \
      --device ~{device} \
      --jit \
      ~{"--max-iter " + max_iter} \
      ~{"--eps-pe " + eps_pe} \
      ~{"--eps-sr1 " + eps_sr1} \
      ~{"--eps-sr2 " + eps_sr2}

    tar czf ~{model_name}.sv_genotype_model.tar.gz svmodel/*
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    #gpuCount: 1
    #gpuType: gpu_type
    #nvidiaDriverVersion: nvidia_driver_version
    #zones: "us-east1-d us-east1-c us-central1-a us-central1-c us-west1-b"
  }
}

task SVGenotype {
  input {
    File vcf
    File model_tar
    Int discrete_samples
    String model_name
    String output_vcf_filename
    String gatk_docker
    String device
    String gpu_type
    String nvidia_driver_version
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4.0,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{output_vcf_filename}"
    File out_index = "~{output_vcf_filename}.tbi"
  }
  command <<<

    set -eo pipefail
    mkdir svmodel
    tar xzf ~{model_tar} svmodel/
    tabix ~{vcf}

    gatk --java-options -Xmx~{java_mem_mb}M SVGenotype \
      -V ~{vcf} \
      --output ~{output_vcf_filename} \
      --discrete-samples ~{discrete_samples} \
      --model-name ~{model_name} \
      --model-dir svmodel \
      --device ~{device} \
      --jit
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    #gpuCount: 1
    #gpuType: gpu_type
    #nvidiaDriverVersion: nvidia_driver_version
    #zones: "us-east1-d us-east1-c us-central1-a us-central1-c us-west1-b"
  }
}
