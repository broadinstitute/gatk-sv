version 1.0

workflow CalcuSampleGeneStats {
  input {
    File sample_gene_gnomad
    File sample_pop_gnomad
    File sample_gene_aou
    File sample_pop_aou
    Array[Int] shard_indices
    String sv_pipeline_base_docker
  }

  scatter (i in shard_indices) {
    call RunSampleGeneStats {
      input:
        sample_gene_gnomad = sample_gene_gnomad,
        sample_pop_gnomad  = sample_pop_gnomad,
        sample_gene_aou    = sample_gene_aou,
        sample_pop_aou     = sample_pop_aou,
        shard_index        = i,
        docker_image       = sv_pipeline_base_docker
    }
  }

  call ConcatTSV as concat_table {
    input:
      tsv_files =  RunSampleGeneStats.out_table,
      output_name = "gnomAD_AoU.sample_gene_table",
      docker_image       = sv_pipeline_base_docker
  }

  call ConcatTSV as concat_stat{
    input:
      tsv_files =  RunSampleGeneStats.stat_table,
      output_name = "gnomAD_AoU.sample_gene_stat",
      docker_image       = sv_pipeline_base_docker
  }

  output {
    File shard_tables = concat_table.combined_tsv
    File shard_stats  = concat_stat.combined_tsv
  }
}

task RunSampleGeneStats {
  input {
    File sample_gene_gnomad
    File sample_pop_gnomad
    File sample_gene_aou
    File sample_pop_aou
    Int shard_index
    String docker_image
  }

  command <<<
    set -euo pipefail

    Rscript -e '

    merge_unique_items <- function(x) {
      x %>%
        strsplit(",") %>%          # split each element by comma
        unlist() %>%               # flatten into one vector
        trimws() %>%               # remove extra whitespace
        unique()                   # keep only unique items
    }

    read_in_stat<-function(sample_pop, sample_gene){
      d1=read.table(sample_pop)
      d2=read.table(sample_gene, sep='\t')
      d2[,2] = toupper(d2[,2])
      dat=merge(d1,d2, by='V1')
      return(dat)
    }

    merge_sample_gene_stat.across_gnomad_aou<-function(sample_pop_gnomad, sample_pop_aou, sample_gene_gnomad, sample_gene_aou){
      dat = read_in_stat(sample_pop_gnomad, sample_gene_gnomad)
      dat2 = read_in_stat(sample_pop_aou, sample_gene_aou)

      dat[dat[,1]%in%dat2[,1],][,1] = paste(dat[dat[,1]%in%dat2[,1],][,1],'a', sep='_')

      out=rbind(dat, dat2)
      out[,2] =toupper(out[,2])
      out[out[,2]%in%c('ASJ','FIN','NFE') ,][,2] = 'EUR'
      out[,2] = factor(out[,2], levels = c('EUR','AFR','AMR','EAS','SAS','AMI','MID','OTH'))
      out = out[order(out[,4], decreasing = T),]
      out = out[order(out[,2]),]
      return(out)
    }

    calcu_sample_gene_stat<-function(out, sample_count_range){
      library(dplyr)
      stat = data.frame(sample_count_range)
      rec = 0
      for(i in sample_count_range){
          rec = rec+1
          tmp = out[1:i,]
          print(c(i,nrow(tmp)))
          gene_list = merge_unique_items(c(tmp[,3]))
          stat[rec,2] = length(gene_list)
      }
      return(stat)
    }


    i <- ~{shard_index}

    out <- merge_sample_gene_stat.across_gnomad_aou( "~{sample_pop_gnomad}", 
                                                      "~{sample_pop_aou}", 
                                                      "~{sample_gene_gnomad}", 
                                                      "~{sample_gene_aou}")

    stat <- calcu_sample_gene_stat( out[out[,4] > 0, ], c(((i - 1) * 500 + 1):(i * 500)) )

    write.table( out, paste0("sample_vs_lof_genes.shard_", i, ".tsv"), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE )

    write.table(stat, paste0("sample_vs_lof_genes.shard_", i, ".stat"), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE )

    write.table( data.frame(table(out[,2])), "sample_vs_pop.stat", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE )
    
    '

  >>>

  output {
    File out_table  = "sample_vs_lof_genes.shard_~{shard_index}.tsv"
    File stat_table = "sample_vs_lof_genes.shard_~{shard_index}.stat"
    File pop_stat   = "sample_vs_pop.stat"
  }

  runtime {
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk 10 HDD"
    docker: docker_image
  }
}

task ConcatTSV {
  input {
    Array[File] tsv_files
    String output_name
    String docker_image
  }

  command <<<
    set -euo pipefail

    # concatenate all shard TSVs
    cat ~{sep=' ' tsv_files} > ~{output_name}
  >>>

  output {
    File combined_tsv = "~{output_name}"
  }

  runtime {
    cpu: 1
    memory: "2 GiB"
    disks: "local-disk 5 HDD"
    docker: docker_image
  }
}


