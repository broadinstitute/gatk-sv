version 1.0

import "ShardedCountSVs.wdl" as count
import "TasksMakeCohortVcf.wdl" as tasks

workflow GetVcfStats {
  input {
    String prefix
    Array[File] vcfs  # sharded by contig, in order

    Boolean count_per_sample = true

    String sites_query_string = "%ID\t%INFO/SVTYPE\t%FILTER\n"
    String bcftools_preprocessing_options = "-i 'FILTER=\"PASS\" || FILTER=\"MULTIALLELIC\"'" # for preprocessing prior to per-sample SV counting

    File contigs_list

    String sv_pipeline_docker
    String sv_base_mini_docker

  }

  Array[String] contigs = read_lines(contigs_list)

  scatter ( i in range(length(vcfs)) ) {
    call GetSitesInfo {
      input:
        vcf = vcfs[i],
        prefix = "~{prefix}.~{contigs[i]}",
        sites_query_string = sites_query_string,
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  call tasks.ZcatCompressedFiles {
    input:
      shards=GetSitesInfo.sites_table,
      outfile_name="~{prefix}.all_contigs.sites_info.tsv.gz",
      sv_base_mini_docker=sv_base_mini_docker
  }

  if (count_per_sample) {
    scatter ( i in range(length(vcfs)) ) {
      call count.ShardedCountSVs {
        input:
          vcf = vcfs[i],
          bcftools_preprocessing_options = bcftools_preprocessing_options,
          prefix = "~{prefix}.~{contigs[i]}",
          sv_pipeline_docker = sv_pipeline_docker,
          sv_base_mini_docker = sv_base_mini_docker
      }
    }
    call count.SumSVCounts {
      input:
        counts = ShardedCountSVs.sv_counts,
        prefix = "~{prefix}.all_contigs",
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  output {
    File? sv_counts = SumSVCounts.sv_counts_summed
    File sites_info = ZcatCompressedFiles.outfile
  }
}


task GetSitesInfo {
  input {
    File vcf
    String prefix
    String sites_query_string
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: 20 + ceil(size([vcf], "GiB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    bcftools query -f '~{sites_query_string}' ~{vcf} | bgzip -c > ~{prefix}.sites_info.tsv.gz
  >>>

  output {
    File sites_table = "~{prefix}.sites_info.tsv.gz"
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
