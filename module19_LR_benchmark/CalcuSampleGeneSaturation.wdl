version 1.0
import "Structs.wdl"


workflow CalcuSampleGeneSaturation {

  input {
    File sample_list
    File bed_file
    File gene_loeuf_file
    Int shard_size = 100
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_split_samples
    RuntimeAttr? runtime_attr_process_loeuf
    RuntimeAttr? runtime_attr_process_shard
    RuntimeAttr? runtime_attr_combine_shard

  }

  call SplitSamples {
    input:
      sample_list = sample_list,
      shard_size = shard_size,
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_split_samples
  }

  call ProcessLOEUF {
    input:
      gene_loeuf_file = gene_loeuf_file,
      output_gene_name = "gene_loeuf_processed.tsv",
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_process_loeuf
  }

  scatter (shard in SplitSamples.sample_shards) {
    call ProcessShard {
      input:
        sample_list = shard,
        input_bed_file = bed_file,
        gene_loeuf = ProcessLOEUF.output_gene, 
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_process_shard
    }
  }

  call CombineShardOutputs {
    input:
      shard_outputs = ProcessShard.sample_gene_per_shard,
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_combine_shard
  }

  output {
    File processed_gene_loeuf = ProcessLOEUF.output_gene
    File combined_sample_gene = CombineShardOutputs.combined_output
  }
}

###############################################################################
# Task 1: Split sample list into shards
###############################################################################
task SplitSamples {

  input {
    File sample_list
    Int shard_size
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(sample_list, "GiB") * 2),
    disk_gb: 15 + ceil(size(sample_list, "GiB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    cut -f1 ~{sample_list} | split -l ~{shard_size} - shard_

    for f in shard_*; do
      echo -e "sample_id\n$(cat $f)" > ${f}.tsv
      rm $f
    done
  >>>

  output {
    Array[File] sample_shards = glob("shard_*.tsv")
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

###############################################################################
# Task 2: Process LOEUF file
###############################################################################
task ProcessLOEUF {

  input {
    File gene_loeuf_file
    String output_gene_name
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(gene_loeuf_file, "GiB") * 2),
    disk_gb: 15 + ceil(size(gene_loeuf_file, "GiB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    Rscript - << 'EOF'
    dat=read.table("~{gene_loeuf_file}", header=TRUE, comment.char="")
    dat=dat[!is.na(dat$LOEUF),]
    dat=dat[order(dat$LOEUF),]
    dat$LOEUF_tile = as.integer((0:(nrow(dat)-1)) / nrow(dat) * 10)
    write.table( dat[,c('gene_name','LOEUF','LOEUF_tile')], "~{output_gene_name}", quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE )
    EOF
  >>>

  output {
    File output_gene = "~{output_gene_name}"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

###############################################################################
# Task 3: Process each shard
###############################################################################
task ProcessShard {

  input {
    File input_bed_file
    File gene_loeuf
    File sample_list
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(input_bed_file, "GiB") * 2),
    disk_gb: 15 + ceil(size(input_bed_file, "GiB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    Rscript - << 'EOF'
    dat=read.table("~{input_bed_file}", header=FALSE)
    aou_sample=read.table("~{sample_list}", header=TRUE)
    loeuf = read.table("~{gene_loeuf}", header=T)

    sample_gene = data.frame(sample_id=aou_sample[,1], stringsAsFactors=FALSE)

    sample_gene$genes = sapply(sample_gene$sample_id, function(x) {
      paste(dat[grepl(x, dat[,6]),][,8], collapse=',')
    })

    sample_gene$gene_count = sapply(sample_gene$genes, function(x) {
      if (x == "") return(0)
      length(unique(strsplit(x, ',')[[1]]))
    })

    for(i in sort(unique(loeuf$LOEUF_tile))){

        gene_list = loeuf[loeuf$LOEUF_tile==i,]$gene
        sample_gene[,ncol(sample_gene)+1] = sapply(sample_gene[,2], function(x){paste(unique(gene_list[gene_list%in%strsplit(as.character(x),',')[[1]]]), collapse=',') })
        colnames(sample_gene)[ncol(sample_gene)] = paste('gene_list.LOEUF_decile', i, sep='_')
        sample_gene[,ncol(sample_gene)+1] = sapply(sample_gene[,2], function(x){length(unique(gene_list[gene_list%in%strsplit(as.character(x),',')[[1]]])) })
        colnames(sample_gene)[ncol(sample_gene)] = paste('gene_count.LOEUF_decile', i, sep='_')
    }

    write.table(sample_gene, "sample_gene_per_shard.tsv", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE )
    
    EOF
  >>>

  output {
    File sample_gene_per_shard = "sample_gene_per_shard.tsv"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

###############################################################################
# Task 4: Combine shard outputs
###############################################################################
task CombineShardOutputs {

  input {
    Array[File] shard_outputs
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(shard_outputs[0], "GiB") * 2),
    disk_gb: 15 + ceil(size(shard_outputs[0], "GiB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    cat ~{sep=' ' shard_outputs} > combined_sample_gene.tsv
  >>>

  output {
    File combined_output = "combined_sample_gene.tsv"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
