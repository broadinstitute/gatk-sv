version 1.0

import "Structs.wdl"
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks

workflow RemovePanelCreationErrorSites {
  input {
    File input_vcf
    File input_vcfs_idx
    File get_panel_creation_error_sites_py
    String sv_pipeline_base_docker
  }

  call RemovePanelCreationError {
        input:
          input_vcf = input_vcf,
          input_vcf_idx = input_vcfs_idx,
          get_panel_creation_error_sites_py = get_panel_creation_error_sites_py,
          docker_image = sv_pipeline_base_docker
  }

  output {
    File panel_error_removed_vcf = RemovePanelCreationError.output_vcf
    File panel_error_removed_idx = RemovePanelCreationError.output_vcf_idx
  }
}



task RemovePanelCreationError {
  input {
    File input_vcf
    File input_vcf_idx
    File get_panel_creation_error_sites_py
    String docker_image
    File? monitoring_script
    RuntimeAttr? runtime_attr_override
  }

  
  String prefix = basename(input_vcf, ".vcf.gz")
  command <<<
    set -euxo pipefail

    # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
    touch monitoring.log
    if [ -s ~{monitoring_script} ]; then
        bash ~{monitoring_script} > monitoring.log &
    fi


    python ~{get_panel_creation_error_sites_py} ~{input_vcf} 31 1 > remove_sites.txt

    bcftools view -T ^remove_sites.txt -S ^remove_samples.txt ~{input_vcf} -Oz -o "~{prefix}.remove_error_sites.vcf.gz"

    bcftools index -t "~{prefix}.remove_error_sites.vcf.gz"

  >>>

  output {
    File output_vcf = "~{prefix}.remove_error_sites.vcf.gz"
    File output_vcf_idx = "~{prefix}.remove_error_sites.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: ceil(size(input_vcf,"GiB")*2)+10,
    disk_gb: ceil(size(input_vcf,"GiB")*2)+20,
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


