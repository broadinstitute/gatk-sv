version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "ClusterSingleChromosome.wdl" as VcfClusterTasks

# Workflow to run parallelized vcf clustering for a single chromosome
workflow VcfClusterSingleChrom {
  input {
    Array[File] vcfs
    Int num_samples
    String prefix
    String evidence_type
    String cohort_name
    Int dist
    Float frac
    Float sample_overlap
    File? exclude_list
    Array[String] batches
    Int sv_size
    Array[String] sv_types
    String contig
    Int localize_shard_size
    Boolean subset_sr_lists
    File bothside_pass
    File background_fail
    File empty_file

    File hail_script
    String project

    String sv_pipeline_docker
    String sv_base_mini_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_localize_vcfs
    RuntimeAttr? runtime_override_join_vcfs
    RuntimeAttr? runtime_override_fix_multiallelic
    RuntimeAttr? runtime_override_fix_ev_tags

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_subset_bothside_pass
    RuntimeAttr? runtime_override_subset_background_fail

    # overrides for VcfClusterTasks
    RuntimeAttr? runtime_override_shard_clusters
    RuntimeAttr? runtime_override_shard_vids
    RuntimeAttr? runtime_override_subset_sv_type
    RuntimeAttr? runtime_override_shard_vcf_precluster
    RuntimeAttr? runtime_override_pull_vcf_shard
    RuntimeAttr? runtime_override_svtk_vcf_cluster
    RuntimeAttr? runtime_override_get_vcf_header_with_members_info_line
    RuntimeAttr? runtime_override_concat_vcf_cluster
    RuntimeAttr? runtime_override_concat_svtypes
    RuntimeAttr? runtime_override_concat_sharded_cluster
    RuntimeAttr? runtime_override_make_sites_only
    RuntimeAttr? runtime_override_sort_merged_vcf

    RuntimeAttr? runtime_override_preconcat_sharded_cluster
    RuntimeAttr? runtime_override_hail_merge_sharded_cluster
    RuntimeAttr? runtime_override_fix_header_sharded_cluster
  }

  scatter (i in range(length(vcfs))) {
    call LocalizeContigVcfs {
      input:
        vcf=vcfs[i],
        vcf_index = vcfs[i] + ".tbi",
        shard_size = localize_shard_size,
        contig=contig,
        prefix=prefix + "." + contig + "." + batches[i],
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_localize_vcfs
    }
  }
  Array[Array[File]] sharded_vcfs_ = transpose(LocalizeContigVcfs.out)

  scatter (i in range(length(sharded_vcfs_))) {
    call JoinVcfs {
      input:
        vcfs=sharded_vcfs_[i],
        contig=contig,
        prefix=prefix + "." + i,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_join_vcfs
    }
    call FixMultiallelicRecords {
      input:
        joined_vcf=JoinVcfs.out,
        batch_contig_vcfs=sharded_vcfs_[i],
        contig=contig,
        prefix=prefix + "." + i,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_fix_multiallelic
    }
    call FixEvidenceTags {
      input:
        vcf=FixMultiallelicRecords.out,
        contig=contig,
        prefix=prefix + "." + i,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_fix_ev_tags
    }
  }

  call MiniTasks.ConcatVcfs {
    input:
      vcfs=FixEvidenceTags.out,
      vcfs_idx=FixEvidenceTags.out_index,
      naive=true,
      outfile_prefix="~{prefix}.precluster_concat",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_vcf_cluster
  }

  #Run vcfcluster per chromosome
  call VcfClusterTasks.ClusterSingleChrom {
    input:
      vcf=ConcatVcfs.concat_vcf,
      vcf_index=ConcatVcfs.concat_vcf_idx,
      num_samples=num_samples,
      contig=contig,
      cohort_name=cohort_name,
      evidence_type=evidence_type,
      prefix=prefix,
      dist=dist,
      frac=frac,
      sample_overlap=sample_overlap,
      exclude_list=exclude_list,
      sv_size=sv_size,
      sv_types=sv_types,
      empty_file=empty_file,
      hail_script=hail_script,
      project=project,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_override_subset_sv_type=runtime_override_subset_sv_type,
      runtime_override_shard_clusters=runtime_override_shard_clusters,
      runtime_override_shard_vids=runtime_override_shard_vids,
      runtime_override_pull_vcf_shard=runtime_override_pull_vcf_shard,
      runtime_override_svtk_vcf_cluster=runtime_override_svtk_vcf_cluster,
      runtime_override_get_vcf_header_with_members_info_line=runtime_override_get_vcf_header_with_members_info_line,
      runtime_override_concat_svtypes=runtime_override_concat_svtypes,
      runtime_override_concat_sharded_cluster=runtime_override_concat_sharded_cluster,
      runtime_override_make_sites_only=runtime_override_make_sites_only,
      runtime_override_sort_merged_vcf=runtime_override_sort_merged_vcf,
      runtime_override_preconcat_sharded_cluster=runtime_override_preconcat_sharded_cluster,
      runtime_override_hail_merge_sharded_cluster=runtime_override_hail_merge_sharded_cluster,
      runtime_override_fix_header_sharded_cluster=runtime_override_fix_header_sharded_cluster
  }

  if(subset_sr_lists) {
    #Subset bothside_pass & background_fail to chromosome of interest
    call SubsetVariantList as SubsetBothsidePass {
      input:
        vid_list=bothside_pass,
        vid_col=2,
        vcf=ConcatVcfs.concat_vcf,
        outfile_name="~{prefix}.pass.VIDs.list",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_subset_bothside_pass
    }
    call SubsetVariantList as SubsetBackgroundFail {
      input:
        vid_list=background_fail,
        vid_col=1,
        vcf=ConcatVcfs.concat_vcf,
        outfile_name="~{prefix}.fail.VIDs.list",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_subset_background_fail
    }
  }

  output {
    Array[File] clustered_vcfs = ClusterSingleChrom.clustered_vcfs
    Array[File] clustered_vcf_indexes = ClusterSingleChrom.clustered_vcf_indexes
    File filtered_bothside_pass = select_first([SubsetBothsidePass.filtered_vid_list, empty_file])
    File filtered_background_fail = select_first([SubsetBackgroundFail.filtered_vid_list, empty_file])
  }
}

task LocalizeContigVcfs {
  input {
    File vcf
    File vcf_index
    Int shard_size
    String contig
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10 + size(vcf, "GiB") * 1.5),
                                  cpu_cores: 1,
                                  preemptible_tries: 1,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    # See Issue #52 "Use GATK to retrieve VCF records in JoinContigFromRemoteVcfs"
    # https://github.com/broadinstitute/gatk-sv/issues/52

    tabix -h "~{vcf}" "~{contig}" \
      | sed "s/AN=[0-9]*;//g" \
      | sed "s/AC=[0-9]*;//g" \
      | bgzip \
      > contig.vcf.gz
    tabix contig.vcf.gz

    python3 <<CODE
    import pysam

    vcf = pysam.VariantFile("contig.vcf.gz")
    SHARD_SIZE = ~{shard_size}

    i = 0
    shard = 0
    vcf_out = None
    for record in vcf:
      if i % SHARD_SIZE == 0:
        if vcf_out is not None:
          vcf_out.close()
        path = f"~{prefix}.{shard:06d}.vcf.gz"
        vcf_out = pysam.VariantFile(path, mode='w', header=vcf.header)
        shard += 1
      vcf_out.write(record)
      i += 1

    if vcf_out is not None:
      vcf_out.close()
    CODE
  >>>

  output {
    Array[File] out = glob("~{prefix}.*.vcf.gz")
  }
}

# Merge contig vcfs across batches
task JoinVcfs {
  input {
    Array[File] vcfs
    String contig
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcfs, "GiB")
  Float input_size_ratio = 3.0
  Float base_disk_gb = 10.0
  RuntimeAttr runtime_default = object {
                                  mem_gb: 1.0,
                                  disk_gb: ceil(base_disk_gb + input_size * input_size_ratio),
                                  cpu_cores: 1,
                                  preemptible_tries: 1,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    python3 <<CODE | bgzip > ~{prefix}.~{contig}.joined.vcf.gz
    import sys
    import gzip

    fl = open("~{write_lines(vcfs)}")
    files = [gzip.open(f.strip(), 'rb') for f in fl.readlines()]
    lines_zip = zip(*files)

    for linesb in lines_zip:
      lines = [l.decode('utf-8') for l in linesb]
      ex = lines[0]
      if ex.startswith('##'):
        sys.stdout.write(ex)
      else:
        sys.stdout.write(ex.strip())
        if len(lines) > 1:
          sys.stdout.write('\t')
          out_lines = [l.strip().split('\t', 9)[-1] for l in lines[1:]]
          sys.stdout.write("\t".join(out_lines))
        sys.stdout.write('\n')
    CODE
    tabix ~{prefix}.~{contig}.joined.vcf.gz
  >>>

  output {
    File out = "~{prefix}.~{contig}.joined.vcf.gz"
    File out_index = "~{prefix}.~{contig}.joined.vcf.gz.tbi"
  }
}

# Add in max CN state to multiallelics
task FixMultiallelicRecords {
  input {
    File joined_vcf
    Array[File] batch_contig_vcfs
    String contig
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(joined_vcf, "GiB") * 2 + size(batch_contig_vcfs, "GiB")
  Float input_size_fraction = 2.0
  Float base_disk_gb = 10.0
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(base_disk_gb + input_size * input_size_fraction),
                                  cpu_cores: 1,
                                  preemptible_tries: 1,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail
    /opt/sv-pipeline/04_variant_resolution/scripts/make_concordant_multiallelic_alts.py \
      ~{joined_vcf} \
      ~{write_lines(batch_contig_vcfs)} \
      ~{prefix}.~{contig}.fixed_multiallelics.vcf.gz
    tabix ~{prefix}.~{contig}.fixed_multiallelics.vcf.gz
  >>>

  output {
    File out = "~{prefix}.~{contig}.fixed_multiallelics.vcf.gz"
    File out_index = "~{prefix}.~{contig}.fixed_multiallelics.vcf.gz.tbi"
  }
}

# Convert EV field from String to Integer
task FixEvidenceTags {
  input {
    File vcf
    String contig
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GiB")
  Float input_size_ratio = 2.0
  Float base_disk_gb = 10.0
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(base_disk_gb + input_size * input_size_ratio),
                                  cpu_cores: 1,
                                  preemptible_tries: 1,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail
    zcat ~{vcf} \
      | sed -e 's/:RD,PE,SR/:7/g' \
      | sed -e 's/:PE,SR/:6/g' \
      | sed -e 's/:RD,SR/:5/g' \
      | sed -e 's/:RD,PE/:3/g' \
      | sed -e 's/:PE\t/:2\t/g' -e 's/:SR\t/:4\t/g' -e 's/:RD\t/:1\t/g' \
      | sed -e 's/ID=EV,Number=.,Type=String/ID=EV,Number=1,Type=Integer/g' \
      | bgzip \
      > ~{prefix}.~{contig}.unclustered.vcf.gz
    tabix ~{prefix}.~{contig}.unclustered.vcf.gz
  >>>

  output {
    File out = "~{prefix}.~{contig}.unclustered.vcf.gz"
    File out_index = "~{prefix}.~{contig}.unclustered.vcf.gz.tbi"
  }
}

# Find intersection of Variant IDs from vid_list with those present in vcf, return as filtered_vid_list
task SubsetVariantList {
  input {
    File vid_list
    Int vid_col
    File vcf
    String outfile_name
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + size(vid_list, "GB") * 2.0 + size(vcf, "GB")),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    zgrep -v "^#" ~{vcf} | cut -f3 > valid_vids.list
    awk -F'\t' -v OFS='\t' 'ARGIND==1{inFileA[$1]; next} {if ($~{vid_col} in inFileA) print }' valid_vids.list ~{vid_list} \
      > ~{outfile_name}
  >>>

  output {
    File filtered_vid_list = outfile_name
  }
}