version 1.0
import "Structs.wdl"

workflow ExtractPerSampleVCFs {
    input {
        File vcf
        File vcf_index
        Array[String]? sample_list   # optional
        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_get_samples
        RuntimeAttr? runtime_attr_extract_sample
    }

    # If no list is provided, query all sample names
    call GetSamples {
        input:
            vcf = vcf,
            vcf_index = vcf_index,
            sample_list = sample_list,
            docker_image = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_get_samples
    }

    scatter (s in GetSamples.sample_ids) {
        call ExtractSample {
            input:
                vcf = vcf,
                vcf_index = vcf_index,
                sample_id = s,
                docker_image = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_extract_sample
        }
    }

    output {
        Array[File] per_sample_vcfs = ExtractSample.per_sample_vcf
        Array[File] per_sample_tbis = ExtractSample.per_sample_tbi
    }
}


task GetSamples {
    input {
        File vcf
        File vcf_index
        Array[String]? sample_list
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        if [[ -z "${sample_list}" ]]; then
            # extract sample names from VCF
            bcftools query -l ~{vcf} > samples.txt
        else
            # write provided list to file
            printf "%s\n" ~{sep="\n" sample_list} > samples.txt
        fi
    >>>

    output {
        Array[String] sample_ids = read_lines("samples.txt")
    }

    RuntimeAttr default_attr = object {
      cpu_cores: 1, 
      mem_gb: 15, 
      disk_gb: 30,
      boot_disk_gb: 10,
      preemptible_tries: 1,
      max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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


task ExtractSample {
    input {
        File vcf
        File vcf_index
        String sample_id
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        bcftools view -s ~{sample_id} -Oz ~{vcf} > ~{sample_id}.vcf.gz
        bcftools index -t ~{sample_id}.vcf.gz
    >>>

    output {
        File per_sample_vcf = "~{sample_id}.vcf.gz"
        File per_sample_tbi = "~{sample_id}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
      cpu_cores: 1, 
      mem_gb: 5, 
      disk_gb: 15,
      boot_disk_gb: 10,
      preemptible_tries: 1,
      max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

