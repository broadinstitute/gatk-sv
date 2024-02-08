version 1.0


import "CalcAF.wdl" as calcAF
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow prune_and_add_vafs {
  input {
    File vcf
    File vcf_idx
    String prefix
    String sv_pipeline_docker

    File? sample_pop_assignments  #Two-column file with sample ID & pop assignment. "." for pop will ignore sample
    File? prune_list              #List of samples to be excluded from the output vcf
    File? famfile                 #Used for M/F AF calculations
    File? par_bed                 #Used to mark hemizygous males on chrX/Y
    Int sv_per_shard
    File contiglist
    String? drop_empty_records
    
  }
  Array[Array[String]] contigs=read_tsv(contiglist)
  #Iterate over chromosomes
  scatter (contig in contigs) {
    #Prune VCF
    call PruneVcf {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        contig=contig[0],
        prune_list=prune_list,
        prefix=prefix,
        sv_pipeline_docker=sv_pipeline_docker
    }
    #Compute AC, AN, and AF per population & sex combination
    call calcAF.CalcAF as getAFs {
      input:
        vcf=PruneVcf.pruned_vcf,
        vcf_idx=PruneVcf.pruned_vcf_idx,
        sv_per_shard=sv_per_shard,
        prefix=prefix,
        sample_pop_assignments=sample_pop_assignments,
        famfile=famfile,
        par_bed=par_bed,
        drop_empty_records=drop_empty_records,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  #Merge pruned VCFs with allele info
  call MiniTasks.ConcatVcfs as concat_vcfs {
    input:
      vcfs=getAFs.vcf_wAFs,
      naive=true,
      outfile_prefix="~{prefix}.pruned_wAFs",
      sv_base_mini_docker=sv_pipeline_docker
  }

  output {
    File output_vcf = concat_vcfs.concat_vcf
    File output_vcf_idx = concat_vcfs.concat_vcf_idx
  }
}


#Shard vcf into single chromosome shards & drop pruned samples
task PruneVcf {
  input{
    File vcf
    File vcf_idx
    String contig
    File? prune_list
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 250,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    #Tabix chromosome of interest
    tabix -h ~{vcf} ~{contig} | bgzip -c > ~{contig}.vcf.gz
    #Get column indexes corresponding to samples to drop, if any exist
    if ~{defined(prune_list)}; then
      dropidx=$( zfgrep "#" ~{contig}.vcf.gz | fgrep -v "##" \
                 | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print NR, $1 }' \
                 | fgrep -wf ~{prune_list} | cut -f1 | paste -s -d, )
      zcat ~{contig}.vcf.gz \
      | cut --complement -f"$dropidx" \
      | bgzip -c \
      > "~{prefix}.~{contig}.pruned.vcf.gz"
    else
      cp "~{contig}.vcf.gz" "~{prefix}.~{contig}.pruned.vcf.gz"
    fi
    tabix -f "~{prefix}.~{contig}.pruned.vcf.gz"
  >>>

  output {
    File pruned_vcf = "~{prefix}.~{contig}.pruned.vcf.gz"
    File pruned_vcf_idx = "~{prefix}.~{contig}.pruned.vcf.gz.tbi"
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
