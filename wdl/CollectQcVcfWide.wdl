version 1.0


import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow CollectQcVcfWide {
  input {
    Array[File] vcfs
    String contig
    Int sv_per_shard
    String prefix

    File ref_fa
    File ref_fai

    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_override_preprocess_vcf
    RuntimeAttr? runtime_override_collect_sharded_vcf_stats
    RuntimeAttr? runtime_override_svtk_vcf_2_bed
    RuntimeAttr? runtime_override_scatter_vcf
    RuntimeAttr? runtime_override_merge_subvcf_stat_shards
    RuntimeAttr? runtime_override_merge_svtk_vcf_2_bed
  }

  String output_prefix = "~{prefix}.collect_qc_vcf_wide"

  scatter ( vcf in vcfs ) {
    call MiniTasks.ScatterVcf {
      input:
        vcf=vcf,
        vcf_index=vcf + ".tbi",
        contig=contig,
        records_per_shard=sv_per_shard,
        prefix="~{output_prefix}.scatter_vcf",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_scatter_vcf
    }
  }
  Array[File] vcf_shards = flatten(ScatterVcf.shards)

  scatter (i in range(length(vcf_shards))) {
    call PreprocessVcf {
      input:
        vcf=vcf_shards[i],
        ref_fa=ref_fa,
        ref_fai=ref_fai,
        prefix="~{output_prefix}.preprocess.shard_~{i}",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_preprocess_vcf
    }

    call CollectShardedVcfStats {
      input:
        vcf=PreprocessVcf.outvcf,
        prefix="~{output_prefix}.collect_stats.shard_~{i}",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_collect_sharded_vcf_stats
      }

    call SvtkVcf2bed {
      input:
        vcf=PreprocessVcf.outvcf,
        prefix="~{output_prefix}.shard_~{i}",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_svtk_vcf_2_bed
      }
    }

  call MiniTasks.ConcatBeds as MergeSubvcfStatShards {
    input:
      shard_bed_files=CollectShardedVcfStats.vcf_stats,
      prefix="~{output_prefix}.VCF_sites.stats",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_merge_subvcf_stat_shards
  }

  call MiniTasks.ConcatBeds as MergeSvtkVcf2bed {
    input:
      shard_bed_files=SvtkVcf2bed.vcf2bed_subworkflow_out,
      prefix="~{output_prefix}.vcf2bed_subworkflow",
      index_output=false,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_merge_svtk_vcf_2_bed
  }

  output {
    File vcf_stats=MergeSubvcfStatShards.merged_bed_file
    File vcf_stats_idx=MergeSubvcfStatShards.merged_bed_idx
    File samples_list=CollectShardedVcfStats.samples_list[0]
    File vcf2bed_out=MergeSvtkVcf2bed.merged_bed_file
  }
}

task PreprocessVcf {
  input {
    File vcf
    File ref_fa
    File ref_fai
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euo pipefail

    bcftools norm \
      -m -any \
      -f ~{ref_fa} \
      -Oz -o unsorted.vcf.gz \
      ~{vcf}
    
    bcftools sort \
      -Oz -o split.vcf.gz
      unsorted.vcf.gz
    
    tabix split.vcf.gz

    cat << 'EOF' > convert_to_symbolic.py
import pysam
import sys

vcf_in = pysam.VariantFile("split.vcf.gz", "r")
vcf_in.header.info.add("SVTYPE", 1, "String", "Type of variant.")
vcf_in.header.info.add("SVLEN", 1, "Integer", "Length of variant.")
vcf_in.header.info.add("END", 1, "Integer", "End position of variant.")

vcf_out = pysam.VariantFile("~{prefix}.vcf.gz", "w", header=vcf_in.header)

for rec in vcf_in:
    allele_type = rec.info["allele_type"].upper()
    allele_type = allele_type[0] if isinstance(allele_type, tuple) else allele_type
    allele_length = abs(len(rec.alts[0]) - len(rec.ref))

    rec.info["SVTYPE"] = allele_type
    rec.info["SVLEN"] = allele_length
    rec.info["END"] = rec.pos + len(ref_len) - 1
    rec.alts = (f"<{svtype}>",)

    vcf_out.write(rec)

vcf_out.close()
vcf_in.close()
EOF

    python convert_to_symbolic.py
    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File outvcf = "~{prefix}.vcf.gz"
    File outvcf_index = "~{prefix}.vcf.gz.tbi"
  }

  RuntimeAttr runtime_default = object {
    mem_gb: 8,
    disk_gb: 5 * ceil(size(vcf, "GiB")) + 10,
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 0,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }
}

task CollectShardedVcfStats {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GiB")
  RuntimeAttr runtime_default = object {
    mem_gb: 1.5 + 2.0 * input_size,
    disk_gb: ceil(10.0 + 2.0 * input_size),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 0,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
  
    /opt/sv-pipeline/scripts/vcf_qc/collectQC.vcf_wide.sh \
      ~{vcf} \
      /opt/sv-pipeline/scripts/vcf_qc/SV_colors.txt \
      collectQC_vcfwide_output/
    
    cp collectQC_vcfwide_output/data/VCF_sites.stats.bed.gz \
      ~{prefix}.VCF_sites.stats.bed.gz
    cp collectQC_vcfwide_output/data/VCF_sites.stats.bed.gz.tbi \
      ~{prefix}.VCF_sites.stats.bed.gz.tbi
    cp collectQC_vcfwide_output/analysis_samples.list \
      ~{prefix}.analysis_samples.list
    tar -czvf ~{prefix}.collectQC_vcfwide_output.tar.gz \
      collectQC_vcfwide_output
  >>>

  output {
    File vcf_stats = "~{prefix}.VCF_sites.stats.bed.gz"
    File vcf_stats_idx = "~{prefix}.VCF_sites.stats.bed.gz.tbi"
    File samples_list = "~{prefix}.analysis_samples.list"
    File vcfwide_tarball = "~{prefix}.collectQC_vcfwide_output.tar.gz"
  }
}

task SvtkVcf2bed {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_file = "~{prefix}.vcf2bed_subworkflow.bed.gz"
  Float input_size = size(vcf, "GiB")
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(10.0 + input_size * 2.0),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 0,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    
    svtk vcf2bed --info ALL ~{vcf} stdout \
      | bgzip -c \
      > "~{output_file}"
  >>>

  output {
    File vcf2bed_subworkflow_out = output_file
  }
}
