version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks


workflow ReviseSVtypeINStoMEIperContig {
  input{
    File vcf
    File vcf_idx
    String prefix
    String contig

    Int max_shards_per_chrom_step1
    Int min_records_per_shard_step1
    Boolean concat_shards = true

    String sv_base_mini_docker
    String sv_pipeline_base_docker
    String sv_pipeline_updates_docker

    RuntimeAttr? runtime_override_split_vcf_to_clean
    RuntimeAttr? runtime_attr_ReviseSVtypeMEI
    RuntimeAttr? runtime_override_combine_step_1_vcfs
  }


  call MiniTasks.ScatterVcf as SplitVcfReviseSVtypeMEI {
      input:
        vcf=vcf,
        prefix="~{prefix}.~{contig}",
        records_per_shard=min_records_per_shard_step1,
        sv_pipeline_docker=sv_pipeline_updates_docker,
        runtime_attr_override=runtime_override_split_vcf_to_clean
  }

  Array[Pair[File, File]] vcf_shards = zip(SplitVcfReviseSVtypeMEI.shards, SplitVcfReviseSVtypeMEI.shards_idx)
  scatter (vcf_shard in vcf_shards) {
      call ReviseSVtypeMEI {
        input:
          vcf = vcf_shard.left,
          vcf_idx = vcf_shard.right,
          sv_pipeline_base_docker = sv_pipeline_base_docker,
          prefix = basename(vcf_shard.left, ".vcf.gz") + ".SVtypeRevised",
          runtime_attr_override = runtime_attr_ReviseSVtypeMEI
      }
  }

  if (concat_shards) {
    call MiniTasks.ConcatVcfs as CombineStep1Vcfs {
      input:
        vcfs=ReviseSVtypeMEI.updated_vcf,
        vcfs_idx=ReviseSVtypeMEI.updated_vcf_idx,
        naive=true,
        outfile_prefix="~{prefix}.~{contig}.SVtypeRevisedINStoMEI",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_combine_step_1_vcfs
    }
  }

  output{
      File? updated_vcf = CombineStep1Vcfs.concat_vcf
      File? updated_vcf_idx = CombineStep1Vcfs.concat_vcf_idx
      Array[File] updated_vcf_shards = ReviseSVtypeMEI.updated_vcf
      Array[File] updated_vcf_shard_idxs = ReviseSVtypeMEI.updated_vcf_idx
  }
}


# Revise svtype of MEIs to SVTYPE=MEI
task ReviseSVtypeMEI {
  input {
    File vcf
    File vcf_idx
    String prefix
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2, 
    disk_gb: 10 + (3 * ceil(size([vcf], "GB"))),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  
  command <<<
    set -eu -o pipefail

    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/revise_MEI_svtypes.py \
      ~{vcf} \
      ~{prefix}.vcf.gz
    tabix -p vcf -f ~{prefix}.vcf.gz
  >>>

  output{
    File updated_vcf = "~{prefix}.vcf.gz"
    File updated_vcf_idx = "~{prefix}.vcf.gz.tbi"
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

