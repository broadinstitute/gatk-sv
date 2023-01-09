#### Copyright (C) 2022 Ryan Collins & the Talkowski Laboratory
####
#### Workflow to compute and filter no-call GT rates from a GATK-SV VCF


version 1.0

import "Structs.wdl"
import "Utils.wdl" as Utils
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow ApplyNoCallRateCutoffs {
  input {
    File vcf
    File? vcf_idx
    Array[String] svtype_list
    Array[String] ncr_filter_field
    Array[Float] NCR_cff_list
    String? global_ncr_filter_field
    Float? global_max_ncr
    Boolean? exclude_CTX

    String sv_pipeline_base_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_index_vcf
    RuntimeAttr? runtime_attr_split_vcf_by_type
    RuntimeAttr? runtime_attr_apply_ncr_filter
    RuntimeAttr? runtime_attr_concat_filtered_vcfs
    RuntimeAttr? runtime_attr_exclude_type_from_vcf
  }
  
  if ( !defined(vcf_idx)){
    call IndexVcf{
      input:
        vcf = vcf,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_index_vcf
    }
  }

  File vcf_index = select_first([IndexVcf.output_vcf_idx, vcf_idx])

  scatter (i in range(length(svtype_list))){
    call SplitVcfByType{
      input:
        vcf = vcf,
        vcf_idx = vcf_index,
        svtype = svtype_list[i],
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_split_vcf_by_type
    }

    call ApplyNCRFilter {
      input:
        vcf=SplitVcfByType.svtype_vcf,
        global_max_ncr=NCR_cff_list[i],
        global_ncr_filter_field=ncr_filter_field[i],
        sv_pipeline_base_docker=sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_apply_ncr_filter
    }
  }

  call MiniTasks.ConcatVcfs as ConcatFilteredVcfsByType {
    input:
      vcfs = ApplyNCRFilter.filtered_vcf,
      vcfs_idx = ApplyNCRFilter.filtered_vcf_idx,
      outfile_prefix=basename(vcf, ".vcf.gz") + ".NCR_filtered.specific_SVTYPE",
      allow_overlaps = true,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_filtered_vcfs
  }

  if (defined(exclude_CTX)){
    call ExcludeTypeAndCtxFromVcf{
      input:
        vcf = vcf,
        vcf_idx = vcf_index,
        svtype_list = svtype_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_exclude_type_from_vcf
    }
  }

  if (!defined(exclude_CTX)){
    call ExcludeTypeFromVcf{
      input:
        vcf = vcf,
        vcf_idx = vcf_index,
        svtype_list = svtype_list,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_exclude_type_from_vcf
    }
 }

  File other_vcf = select_first([ExcludeTypeAndCtxFromVcf.oth_wo_ctx_vcf, ExcludeTypeFromVcf.oth_vcf])
  File other_vcf_idx = select_first([ExcludeTypeAndCtxFromVcf.oth_wo_ctx_vcf_idx, ExcludeTypeFromVcf.oth_vcf_idx])


  File global_ncr_filter_field_oth = select_first([global_ncr_filter_field, "NCR"])

  call ApplyNCRFilter as apply_ncr_filter_oth{
      input:
        vcf = other_vcf,
        global_max_ncr = global_max_ncr,
        global_ncr_filter_field = global_ncr_filter_field_oth,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_apply_ncr_filter
  }

  call MiniTasks.ConcatVcfs {
    input:
      vcfs=[ConcatFilteredVcfsByType.concat_vcf, apply_ncr_filter_oth.filtered_vcf],
      vcfs_idx = [ConcatFilteredVcfsByType.concat_vcf_idx, apply_ncr_filter_oth.filtered_vcf_idx],
      outfile_prefix=basename(vcf, ".vcf.gz") + ".NCR_filtered",
      allow_overlaps = true,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_filtered_vcfs
  }

  output{
    File ncr_filtered_vcf = ConcatVcfs.concat_vcf
    File ncr_filtered_vcf_idx = ConcatVcfs.concat_vcf_idx
  }
}




task SplitVcfByType{
    input{
        File vcf
        File vcf_idx
        String svtype
        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: ceil(10.0 + size(vcf, "GB") * 2.0),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String prefix = basename(vcf,"vcf.gz")
    command<<<
        set -euo pipefail
        python <<CODE

        import os
        import pysam
        fin=pysam.VariantFile("~{vcf}")
        fo=pysam.VariantFile("~{prefix}.~{svtype}.vcf.gz", 'w', header = fin.header)
        for record in fin:
          if record.info['SVTYPE'] == "~{svtype}":
            fo.write(record)
        fin.close()
        fo.close()
        CODE

        tabix -p vcf "~{prefix}.~{svtype}.vcf.gz"

    >>>

    output{
        File svtype_vcf = "~{prefix}.~{svtype}.vcf.gz"
        File svtype_vcf_idx = "~{prefix}.~{svtype}.vcf.gz.tbi"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_base_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ExcludeTypeFromVcf{
    input{
        File vcf
        File vcf_idx
        Array[String] svtype_list
        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: ceil(10.0 + size(vcf, "GB") * 2.0),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String prefix = basename(vcf,"vcf.gz")
    command<<<
        set -euo pipefail
        python <<CODE

        import os
        import pysam
        fin=pysam.VariantFile("~{vcf}")
        fo=pysam.VariantFile("~{prefix}.OTH.vcf.gz", 'w', header = fin.header)
        for record in fin:
          if not record.info['SVTYPE'] in ~{svtype_list}:
            fo.write(record)
        fin.close()
        fo.close()
        CODE

        tabix -p vcf "~{prefix}.OTH.vcf.gz"

    >>>

    output{
        File oth_vcf = "~{prefix}.OTH.vcf.gz"
        File oth_vcf_idx = "~{prefix}.OTH.vcf.gz.tbi"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_base_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ExcludeTypeAndCtxFromVcf{
    input{
        File vcf
        File vcf_idx
        Array[String] svtype_list
        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: ceil(10.0 + size(vcf, "GB") * 2.0),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String prefix = basename(vcf,"vcf.gz")
    command<<<
        set -euo pipefail


        python <<CODE

        import os
        import pysam
        fa = open("~{write_lines(svtype_list)}")
        svtype_list = []
        for line in fa:
          pin=line.strip().split()
          svtype_list+=pin
        fa.close()

        fin=pysam.VariantFile("~{vcf}")
        fo=pysam.VariantFile("~{prefix}.OTH.vcf.gz", 'w', header = fin.header)
        for record in fin:
          if not record.info['SVTYPE'] in svtype_list and not record.info['SVTYPE']=="CTX":
            fo.write(record)
        fin.close()
        fo.close()
        CODE

        tabix -p vcf "~{prefix}.OTH.vcf.gz"

    >>>

    output{
        File oth_wo_ctx_vcf = "~{prefix}.OTH.vcf.gz"
        File oth_wo_ctx_vcf_idx = "~{prefix}.OTH.vcf.gz.tbi"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_base_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

# Apply filter status to records in VCF based on NCR
task ApplyNCRFilter {
  input {
    File vcf
    File? vcf_idx
    Float? global_max_ncr
    String global_ncr_filter_field
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(vcf, ".vcf.gz")

  RuntimeAttr default_attr = object {
      cpu_cores: 1, 
      mem_gb: 5, 
      disk_gb: ceil(10.0 + size(vcf, "GB") * 2.0),
      boot_disk_gb: 10,
      preemptible_tries: 1,
      max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  command <<<
    set -eu -o pipefail

    if [ ~{defined(vcf_idx)} == "false" ]; then
      tabix -f ~{vcf}
    fi

    script_options=""
    if [ ~{defined(global_max_ncr)} == "true" ]; then
      script_options="$script_options --global-max-ncr ~{global_max_ncr}"
    fi

    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/nocall_rate_filter.py \
      --verbose \
      --global-filter-on ~{global_ncr_filter_field} \
      $script_options \
      ~{vcf} \
      ~{prefix}.NCR_filtered.vcf.gz

    tabix -p vcf "~{prefix}.NCR_filtered.vcf.gz"
  >>>

  output {
    File filtered_vcf = "~{prefix}.NCR_filtered.vcf.gz"
    File filtered_vcf_idx = "~{prefix}.NCR_filtered.vcf.gz.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_base_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

task IndexVcf{
    input{
        File vcf
        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: ceil(10.0 + size(vcf, "GB") * 2),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(vcf,".vcf.gz")
    command<<<
        set -euo pipefail
        gsutil cp ~{vcf} ./
        tabix -p vcf "~{prefix}.vcf.gz"

    >>>

    output{
        File output_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_base_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
