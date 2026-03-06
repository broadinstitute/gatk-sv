version 1.0

import "Structs.wdl"
import "FormatVcfForGatk.wdl" as format_vcf
import "TasksClusterBatch.wdl" as tasks_cluster
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow MergeBatchSites {
  input {
    String cohort             # Cohort name or project prefix for all cohort-level outputs
    Array[File] ploidy_tables
    Array[File] depth_vcfs    # Filtered depth VCFs across batches
    Array[File] pesr_vcfs     # Filtered PESR VCFs across batches

    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    String sv_base_mini_docker
    String sv_pipeline_docker
    String gatk_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_svtk_to_gatk_vcf
    RuntimeAttr? runtime_attr_merge_pesr
    RuntimeAttr? runtime_attr_merge_depth
    RuntimeAttr? runtime_override_concat
    RuntimeAttr? runtime_attr_convert_bnds
  }

  scatter (i in range(length(depth_vcfs))) {
    call format_vcf.FormatVcf as FormatDepth {
      input:
        vcf=depth_vcfs[i],
        ploidy_table=ploidy_tables[i],
        output_prefix=basename(depth_vcfs[i]) + ".reformatted",
        args="--add-sr-pos",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_svtk_to_gatk_vcf
    }
  }

  scatter (i in range(length(pesr_vcfs))) {
    call ConvertBreakends {
      input:
        vcf=pesr_vcfs[i],
        prefix=basename(pesr_vcfs[i]) + ".bnd_converted",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_convert_bnds
    }
    call format_vcf.FormatVcf as FormatPesr {
      input:
        vcf=ConvertBreakends.converted_vcf,
        ploidy_table=ploidy_tables[i],
        output_prefix=basename(pesr_vcfs[i]) + ".reformatted",
        args="--add-sr-pos",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_svtk_to_gatk_vcf
    }
  }

  call tasks_cluster.SVCluster as MergePesr {
    input:
      vcfs=FormatPesr.out,
      output_prefix="~{cohort}.merge_batch_sites.pesr",
      sites_only=true,
      omit_members=true,
      algorithm="SINGLE_LINKAGE",
      depth_sample_overlap=0,
      depth_interval_overlap=1,
      depth_size_similarity=1,
      depth_breakend_window=0,
      mixed_sample_overlap=0,
      mixed_interval_overlap=1,
      mixed_size_similarity=1,
      mixed_breakend_window=0,
      pesr_sample_overlap=0,
      pesr_interval_overlap=1,
      pesr_size_similarity=1,
      pesr_breakend_window=0,
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_fasta_fai,
      reference_dict=reference_dict,
      java_mem_fraction=java_mem_fraction,
      gatk_docker=gatk_docker,
      runtime_attr_override=runtime_attr_merge_pesr
  }

  call tasks_cluster.SVCluster as MergeDepth {
    input:
      vcfs=FormatDepth.out,
      output_prefix="~{cohort}.merge_batch_sites.depth",
      sites_only=true,
      omit_members=true,
      algorithm="SINGLE_LINKAGE",
      depth_sample_overlap=0,
      depth_interval_overlap=1,
      depth_size_similarity=1,
      depth_breakend_window=0,
      mixed_sample_overlap=0,
      mixed_interval_overlap=1,
      mixed_size_similarity=1,
      mixed_breakend_window=0,
      pesr_sample_overlap=0,
      pesr_interval_overlap=1,
      pesr_size_similarity=1,
      pesr_breakend_window=0,
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_fasta_fai,
      reference_dict=reference_dict,
      java_mem_fraction=java_mem_fraction,
      gatk_docker=gatk_docker,
      runtime_attr_override=runtime_attr_merge_depth
  }

  call tasks_cohort.ConcatVcfs {
    input:
      vcfs = [MergePesr.out, MergeDepth.out],
      vcfs_idx = [MergePesr.out_index, MergeDepth.out_index],
      allow_overlaps = true,
      outfile_prefix = cohort + ".merge_batch_sites",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_concat
  }

  output {
    File merge_batch_sites_vcf = ConcatVcfs.concat_vcf
    File merge_batch_sites_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}


task ConvertBreakends {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: ceil(10 + 3 * size(vcf, "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File converted_vcf = "${prefix}.vcf.gz"
  }
  command <<<

    set -euo pipefail

    python3 <<CODE
import gzip
with gzip.open("~{vcf}", 'rt') as f, open("~{prefix}.vcf", 'w') as out:
    for line in f:
        if line.startswith('#'):
            out.write(line)
            continue
        columns = line.split('\t', 8)
        alt = columns[4]
        if "]" in alt or "[" in alt:
            columns[4] = "<BND>"
        out.write("\t".join(columns))
CODE
bgzip ~{prefix}.vcf

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
