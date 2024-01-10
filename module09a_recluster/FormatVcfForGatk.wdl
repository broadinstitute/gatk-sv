version 1.0

import "Structs.wdl"
import "TasksClusterBatch.wdl" as tasks_cluster
import "TasksMakeCohortVcf.wdl" as tasks

workflow FormatVcfForGatk {
  input {
    File vcf
    String prefix
    File ped_file
    Int records_per_shard = 40000

    File contig_list
    File? contigs_header  # Replaces vcf contig dictionary if provided
    String? formatter_args

    String? chr_x
    String? chr_y

    File? svtk_to_gatk_script  # For debugging

    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_create_ploidy
    RuntimeAttr? runtime_attr_scatter
    RuntimeAttr? runtime_attr_format
    RuntimeAttr? runtime_override_concat
    RuntimeAttr? runtime_override_preconcat_step1
    RuntimeAttr? runtime_override_hail_merge_step1
    RuntimeAttr? runtime_override_fix_header_step1
  }

  call tasks_cluster.CreatePloidyTableFromPed {
    input:
      ped_file=ped_file,
      contig_list=contig_list,
      retain_female_chr_y=false,
      chr_x=chr_x,
      chr_y=chr_y,
      output_prefix="~{prefix}.ploidy",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_create_ploidy
  }

  call tasks.ScatterVcf {
    input:
      vcf=vcf,
      records_per_shard = records_per_shard,
      prefix = "~{prefix}.scatter_vcf",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_scatter
  }

  scatter ( i in range(length(ScatterVcf.shards)) ) {
    call FormatVcf {
      input:
        vcf=ScatterVcf.shards[i],
        ploidy_table=CreatePloidyTableFromPed.out,
        args=formatter_args,
        output_prefix="~{prefix}.format.shard_~{i}",
        contigs_header=contigs_header,
        script=svtk_to_gatk_script,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_format
    }
  }

  Boolean shards_unsorted = defined(contigs_header)
  call tasks.ConcatVcfs {
    input:
      vcfs=FormatVcf.out,
      vcfs_idx=FormatVcf.out_index,
      naive=!shards_unsorted,
      allow_overlaps=shards_unsorted,
      outfile_prefix="~{prefix}.gatk_formatted",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat
  }

  output {
    File gatk_formatted_vcf = ConcatVcfs.concat_vcf
    File gatk_formatted_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}

task FormatVcf {
  input {
    File vcf
    File ploidy_table
    File? script
    String? args
    File? contigs_header  # Overwrites contig dictionary, in case they are out of order
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(50 + size(vcf, "GB") * 3),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{output_prefix}.vcf.gz"
    File out_index = "~{output_prefix}.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

    # Convert format
    python ~{default="/opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py" script} \
      --vcf ~{vcf} \
      --out tmp.vcf.gz \
      --ploidy-table ~{ploidy_table} \
      ~{args}

    if ~{defined(contigs_header)}; then
      bcftools view --no-version -h tmp.vcf.gz > original_header.vcf
      grep -v "^##contig=" original_header.vcf | grep -v "^#CHROM" > header.vcf
      cat ~{contigs_header} >> header.vcf
      grep "^#CHROM" original_header.vcf >> header.vcf
      bcftools reheader -h header.vcf tmp.vcf.gz | bcftools sort -Oz -o ~{output_prefix}.vcf.gz
    else
      mv tmp.vcf.gz ~{output_prefix}.vcf.gz
    fi

    tabix ~{output_prefix}.vcf.gz
  >>>
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