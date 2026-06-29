version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

import "GenotypeCpxCnvs.wdl" as GenotypeCpx
import "TasksMakeCohortVcf.wdl" as MiniTasks

# Workflow to perform depth-based genotyping for a single vcf shard scattered 
# across batches on predicted CPX CNVs
workflow ScatterCpxGenotyping {
  input {
    File bin_exclude
    File vcf
    Int records_per_shard
    File batches_file
    File coverage_files_file
    File rd_depth_sep_cutoff_files_file
    File ped_file
    File median_coverage_files_file
    Int n_per_split_small
    Int n_per_split_large
    Int n_rd_test_bins
    Int? min_ddup_thresh
    String prefix
    String contig
    File ref_dict

    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_split_vcf_to_genotype
    RuntimeAttr? runtime_override_concat_cpx_cnv_vcfs

    RuntimeAttr? runtime_override_preconcat
    RuntimeAttr? runtime_override_fix_header

    # overrides for GenotypeCpx
    RuntimeAttr? runtime_override_ids_from_median
    RuntimeAttr? runtime_override_get_cpx_cnv_intervals
    RuntimeAttr? runtime_override_parse_genotypes
    RuntimeAttr? runtime_override_merge_melted_gts
    RuntimeAttr? runtime_override_split_bed_by_size
    RuntimeAttr? runtime_override_rd_genotype
    RuntimeAttr? runtime_override_concat_melted_genotypes
  }

  String contig_prefix = prefix + "." + contig

  # Shard VCF into even slices
  call SelectCandidatesAndScatter {
    input:
      vcf=vcf,
      prefix=contig_prefix,
      records_per_shard=records_per_shard,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_split_vcf_to_genotype
  }

  # Scatter genotyping over shards
  scatter ( shard in SelectCandidatesAndScatter.shards ) {
    # Run genotyping
    call GenotypeCpx.GenotypeCpxCnvs as GenotypeShard {
      input:
        bin_exclude=bin_exclude,
        vcf=shard,
        batches_file=batches_file,
        coverage_files_file=coverage_files_file,
        rd_depth_sep_cutoff_files_file=rd_depth_sep_cutoff_files_file,
        ped_file=ped_file,
        median_coverage_files_file=median_coverage_files_file,
        n_per_split_large=n_per_split_large,
        n_per_split_small=n_per_split_small,
        n_rd_test_bins=n_rd_test_bins,
        min_ddup_thresh=min_ddup_thresh,
        prefix=prefix,
        ped_file=ped_file,
        contig=contig,
        ref_dict=ref_dict,
        linux_docker=linux_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_override_ids_from_median=runtime_override_ids_from_median,
        runtime_override_get_cpx_cnv_intervals=runtime_override_get_cpx_cnv_intervals,
        runtime_override_parse_genotypes=runtime_override_parse_genotypes,
        runtime_override_merge_melted_gts=runtime_override_merge_melted_gts,
        runtime_override_split_bed_by_size=runtime_override_split_bed_by_size,
        runtime_override_rd_genotype=runtime_override_rd_genotype,
        runtime_override_concat_melted_genotypes=runtime_override_concat_melted_genotypes
    }
  }

  call MiniTasks.ConcatVcfs as ConcatCpxCnvVcfs {
    input:
      vcfs=flatten([GenotypeShard.cpx_depth_gt_resolved_vcf, [SelectCandidatesAndScatter.non_candidates_vcf]]),
      vcfs_idx=flatten([GenotypeShard.cpx_depth_gt_resolved_vcf_idx, [SelectCandidatesAndScatter.non_candidates_vcf_index]]),
      allow_overlaps=true,
      outfile_prefix="~{prefix}.regenotyped",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_cpx_cnv_vcfs
  }

  # Output merged VCF
  output {
    File cpx_depth_gt_resolved_vcf = ConcatCpxCnvVcfs.concat_vcf
    File cpx_depth_gt_resolved_vcf_idx = ConcatCpxCnvVcfs.concat_vcf_idx
  }
 }

 # Note: requires docker with updated bcftools
task SelectCandidatesAndScatter {
  input {
    File vcf
    File? vcf_index
    String prefix
    Int records_per_shard
    Int? threads = 1
    String? contig
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
                                  cpu_cores: 2,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    # select candidate CPX
    zcat ~{vcf} | awk '
    {
        if (/^#/ || (/CPX/ && !/UNRESOLVED/) || (/INVERSION_SINGLE_ENDER/ && /UNRESOLVED/)) {
            print | "bgzip -c > candidates.vcf.gz"
        }

        if (/^#/ || (!/CPX/ && !(/INVERSION_SINGLE_ENDER/ && /UNRESOLVED/)) || (/CPX/ && /UNRESOLVED/)) {
            print | "bgzip -c > noncandidates.vcf.gz"
        }
    }'
    tabix candidates.vcf.gz

    # in case the file is empty create an empty shard
    bcftools view -h candidates.vcf.gz | bgzip -c > "~{prefix}.0.vcf.gz"

    # scatter the candidates vcf
    bcftools +scatter candidates.vcf.gz -o . -O z -p "~{prefix}". --threads ~{threads} -n ~{records_per_shard} ~{"-r " + contig}

    ls "~{prefix}".*.vcf.gz | sort -k1,1V > vcfs.list
    i=0
    while read VCF; do
      shard_no=`printf %06d $i`
      mv "$VCF" "~{prefix}.shard_${shard_no}.vcf.gz"
      i=$((i+1))
    done < vcfs.list

    # rename noncandidates vcf AFTER handling shard naming to prevent lumping together
    mv noncandidates.vcf.gz ~{prefix}.noncandidates.vcf.gz
    tabix ~{prefix}.noncandidates.vcf.gz
  >>>
  output {
    Array[File] shards = glob("~{prefix}.shard_*.vcf.gz")
    File non_candidates_vcf = "~{prefix}.noncandidates.vcf.gz"
    File non_candidates_vcf_index = "~{prefix}.noncandidates.vcf.gz.tbi"
  }
}
