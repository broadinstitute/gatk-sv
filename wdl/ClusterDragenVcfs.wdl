version 1.0

import "Structs.wdl"
import "CombineBatches.wdl" as Combine
import "FormatVcfForGatk.wdl" as GatkFormatting
import "TasksClusterBatch.wdl" as ClusterTasks

workflow ClusterDragenVcfs {
  input {
    String cohort_name
    File ped_file

    # Raw DRAGEN inputs, one VCF per batch/input set
    Array[File] dragen_sv_vcfs
    Array[File] dragen_cnv_vcfs

    # Use the same primary-contig fai for ploidy generation and svtk standardization.
    File contig_list
    Int min_size = 50

    # Grouped clustering resources. Provide the downloaded workspace resources here.
    File clustering_config_part1
    File stratification_config_part1
    File clustering_config_part2
    File stratification_config_part2
    Array[String] track_names
    Array[File] track_bed_files

    File reference_fasta
    File reference_fasta_fai
    File reference_dict
    String? chr_x
    String? chr_y

    Float grouped_sample_overlap = 0
    Float? java_mem_fraction

    # Use a GATK image built from the current master branch.
    String gatk_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_create_ploidy
    RuntimeAttr? runtime_attr_standardize_sv
    RuntimeAttr? runtime_attr_standardize_cnv
    RuntimeAttr? runtime_attr_format_sv
    RuntimeAttr? runtime_attr_format_cnv
    RuntimeAttr? runtime_attr_rewrite_config_part1
    RuntimeAttr? runtime_attr_rewrite_config_part2
    RuntimeAttr? runtime_attr_join_vcfs
    RuntimeAttr? runtime_attr_cluster_sites
    RuntimeAttr? runtime_attr_recluster_part1
    RuntimeAttr? runtime_attr_recluster_part2
  }

  call ClusterTasks.CreatePloidyTableFromPed {
    input:
      ped_file=ped_file,
      contig_list=contig_list,
      retain_female_chr_y=true,
      chr_x=chr_x,
      chr_y=chr_y,
      output_prefix="~{cohort_name}.ploidy",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_create_ploidy
  }

  scatter (i in range(length(dragen_sv_vcfs))) {
    call StandardizeVcf as StandardizeDragenSv {
      input:
        vcf=dragen_sv_vcfs[i],
        caller="dragen",
        output_prefix="~{cohort_name}.dragen_sv.~{i}.standardized",
        contig_list=contig_list,
        min_size=min_size,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_standardize_sv
    }

    call GatkFormatting.FormatVcf as FormatDragenSv {
      input:
        vcf=StandardizeDragenSv.out,
        ploidy_table=CreatePloidyTableFromPed.out,
        output_prefix="~{cohort_name}.dragen_sv.~{i}.gatk",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_format_sv
    }
  }

  scatter (i in range(length(dragen_cnv_vcfs))) {
    # DRAGEN CNV is currently registered as "depth" in this checkout.
    call StandardizeVcf as StandardizeDragenCnv {
      input:
        vcf=dragen_cnv_vcfs[i],
        caller="depth",
        output_prefix="~{cohort_name}.dragen_cnv.~{i}.standardized",
        contig_list=contig_list,
        min_size=min_size,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_standardize_cnv
    }

    call GatkFormatting.FormatVcf as FormatDragenCnv {
      input:
        vcf=StandardizeDragenCnv.out,
        ploidy_table=CreatePloidyTableFromPed.out,
        output_prefix="~{cohort_name}.dragen_cnv.~{i}.gatk",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_format_cnv
    }
  }

  call RewriteGroupedClusteringConfig as RewriteClusteringConfigPart1 {
    input:
      config=clustering_config_part1,
      sample_overlap=grouped_sample_overlap,
      output_prefix="~{cohort_name}.cluster_config.part1",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_rewrite_config_part1
  }

  call RewriteGroupedClusteringConfig as RewriteClusteringConfigPart2 {
    input:
      config=clustering_config_part2,
      sample_overlap=grouped_sample_overlap,
      output_prefix="~{cohort_name}.cluster_config.part2",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_rewrite_config_part2
  }

  call ClusterTasks.SVCluster as JoinVcfs {
    input:
      vcfs=flatten([FormatDragenSv.out, FormatDragenCnv.out]),
      ploidy_table=CreatePloidyTableFromPed.out,
      output_prefix="~{cohort_name}.join_vcfs",
      fast_mode=false,
      pesr_sample_overlap=0,
      pesr_interval_overlap=1,
      pesr_breakend_window=0,
      depth_sample_overlap=0,
      depth_interval_overlap=1,
      depth_breakend_window=0,
      mixed_sample_overlap=0,
      mixed_interval_overlap=1,
      mixed_breakend_window=0,
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_fasta_fai,
      reference_dict=reference_dict,
      java_mem_fraction=java_mem_fraction,
      gatk_docker=gatk_docker,
      runtime_attr_override=runtime_attr_join_vcfs
  }

  call ClusterTasks.SVCluster as ClusterSites {
    input:
      vcfs=[JoinVcfs.out],
      ploidy_table=CreatePloidyTableFromPed.out,
      output_prefix="~{cohort_name}.cluster_sites",
      fast_mode=false,
      breakpoint_summary_strategy="REPRESENTATIVE",
      pesr_sample_overlap=0,
      pesr_interval_overlap=0.1,
      pesr_breakend_window=300,
      depth_sample_overlap=0,
      depth_interval_overlap=0.5,
      depth_breakend_window=500000,
      mixed_sample_overlap=0,
      mixed_interval_overlap=0.5,
      mixed_breakend_window=1000000,
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_fasta_fai,
      reference_dict=reference_dict,
      java_mem_fraction=java_mem_fraction,
      variant_prefix="~{cohort_name}_",
      gatk_docker=gatk_docker,
      runtime_attr_override=runtime_attr_cluster_sites
  }

  call Combine.GroupedSVClusterTask as GroupedSVClusterPart1 {
    input:
      vcf=ClusterSites.out,
      ploidy_table=CreatePloidyTableFromPed.out,
      output_prefix="~{cohort_name}.recluster_part_1",
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_fasta_fai,
      reference_dict=reference_dict,
      clustering_config=RewriteClusteringConfigPart1.out,
      stratification_config=stratification_config_part1,
      track_bed_files=track_bed_files,
      track_names=track_names,
      breakpoint_summary_strategy="REPRESENTATIVE",
      java_mem_fraction=java_mem_fraction,
      gatk_docker=gatk_docker,
      runtime_attr_override=runtime_attr_recluster_part1
  }

  call Combine.GroupedSVClusterTask as GroupedSVClusterPart2 {
    input:
      vcf=GroupedSVClusterPart1.out,
      ploidy_table=CreatePloidyTableFromPed.out,
      output_prefix="~{cohort_name}.recluster_part_2",
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_fasta_fai,
      reference_dict=reference_dict,
      clustering_config=RewriteClusteringConfigPart2.out,
      stratification_config=stratification_config_part2,
      track_bed_files=track_bed_files,
      track_names=track_names,
      breakpoint_summary_strategy="REPRESENTATIVE",
      java_mem_fraction=java_mem_fraction,
      gatk_docker=gatk_docker,
      runtime_attr_override=runtime_attr_recluster_part2
  }

  output {
    File clustered_vcf = GroupedSVClusterPart2.out
    File clustered_vcf_index = GroupedSVClusterPart2.out_index
  }
}

task StandardizeVcf {
  input {
    File vcf
    String caller
    String output_prefix
    File contig_list
    Int min_size = 50
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + size(vcf, "GB") * 3),
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

    svtk standardize \
      --contigs ~{contig_list} \
      --min-size ~{min_size} \
      ~{vcf} \
      ~{output_prefix}.vcf.gz \
      ~{caller}

    tabix -p vcf ~{output_prefix}.vcf.gz
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

task RewriteGroupedClusteringConfig {
  input {
    File config
    Float sample_overlap = 0
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{output_prefix}.tsv"
  }
  command <<<
    set -euo pipefail

    python <<'PY'
import csv

input_path = "~{config}"
output_path = "~{output_prefix}.tsv"
sample_overlap = "~{sample_overlap}"

with open(input_path, newline="") as src, open(output_path, "w", newline="") as dst:
    writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
    reader = csv.reader(src, delimiter="\t")
    header = None
    sample_overlap_idx = None
    for row in reader:
      if not row:
        writer.writerow(row)
        continue
      if row[0].startswith("#") and header is None:
        writer.writerow(row)
        continue
      if header is None:
        header = row
        writer.writerow(header)
        sample_overlap_idx = header.index("SAMPLE_OVERLAP")
        continue
        if not row:
            writer.writerow(row)
            continue
        row[sample_overlap_idx] = sample_overlap
        writer.writerow(row)

  if header is None:
    raise ValueError("Failed to find clustering config header")
PY
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