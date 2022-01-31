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

    String sv_base_mini_docker

    RuntimeAttr? runtime_override_split_vcf_to_clean
    RuntimeAttr? runtime_attr_ReviseSVtypeMEI
    RuntimeAttr? runtime_override_combine_step_1_vcfs
  }


  call MiniTasks.SplitVcf as SplitVcfReviseSVtypeMEI {
      input:
        vcf=vcf,
        contig=contig,
        prefix="~{prefix}.~{contig}.shard_",
        n_shards=max_shards_per_chrom_step1,
        min_vars_per_shard=min_records_per_shard_step1,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_split_vcf_to_clean
  }

  scatter (vcf_shard in SplitVcfReviseSVtypeMEI.vcf_shards) {
      call ReviseSVtypeMEI{
        input:
          vcf = vcf_shard,
          sv_base_mini_docker = sv_base_mini_docker,
          prefix = "~{prefix}.~{contig}.SVtypeRevised.shard_",
          runtime_attr_override = runtime_attr_ReviseSVtypeMEI
      }
  }

  call MiniTasks.ConcatVcfs as CombineStep1Vcfs {
      input:
        vcfs=ReviseSVtypeMEI.updated_vcf,
        vcfs_idx=ReviseSVtypeMEI.updated_vcf_idx,
        naive=true,
        outfile_prefix="~{prefix}.~{contig}.SVtypeRevisedINStoMEI",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_combine_step_1_vcfs
  }

  output{
      File updated_vcf = CombineStep1Vcfs.concat_vcf
      File updated_vcf_idx = CombineStep1Vcfs.concat_vcf_idx
  }
}



# revise svtype of MEIs to SVTYPE=MEI
task ReviseSVtypeMEI{
  input{
    File vcf
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  
  command <<<
    zcat ~{vcf} | grep '#' > ~{prefix}.vcf
    zcat ~{vcf} | grep -v '#' | grep "INS:ME" | sed -e "s/SVTYPE=INS/SVTYPE=MEI/" >> ~{prefix}.vcf
    zcat ~{vcf} | grep -v '#' | grep -v "INS:ME"  >> ~{prefix}.vcf
    mkdir tmp
    vcf-sort -t tmp/ ~{prefix}.vcf | bgzip > ~{prefix}.vcf.gz
    tabix -p vcf ~{prefix}.vcf.gz
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
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}





