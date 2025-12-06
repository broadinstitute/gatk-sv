
version 1.0

import "Structs.wdl"
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks

workflow IntegrateUpsetStat {
    input {
        String     sample
        Array[File] vcf_list
        String sv_base_mini_docker
        String sv_pipeline_base_docker
    }

  call CalcuUpsetStat1{
      input:
          sample = sample,
          vcf_list = vcf_list,
          docker_image = sv_pipeline_base_docker
    }

  call CalcuUpsetStat2{
      input:
          sample = sample,
          vcf_list = vcf_list,
          docker_image = sv_pipeline_base_docker
    }

  call CalcuUpsetStat3{
      input:
          sample = sample,
          vcf_list = vcf_list,
          docker_image = sv_pipeline_base_docker
    }

  call CalcuUpsetStat4{
      input:
          sample = sample,
          vcf_list = vcf_list,
          docker_image = sv_pipeline_base_docker
    }

  call CalcuUpsetStat5{
      input:
          sample = sample,
          vcf_list = vcf_list,
          docker_image = sv_pipeline_base_docker
    }

  call CalcuUpsetStat6{
      input:
          sample = sample,
          vcf_list = vcf_list,
          docker_image = sv_pipeline_base_docker
    }

  call CalcuUpsetStat7{
      input:
          sample = sample,
          vcf_list = vcf_list,
          docker_image = sv_pipeline_base_docker
    }

  call CalcuUpsetStat8{
      input:
          sample = sample,
          vcf_list = vcf_list,
          docker_image = sv_pipeline_base_docker
    }

  output{
      File stat1 = CalcuUpsetStat1.inte_stat
      File stat2 = CalcuUpsetStat2.inte_stat
      File stat3 = CalcuUpsetStat3.inte_stat
      File stat4 = CalcuUpsetStat4.inte_stat
      File stat5 = CalcuUpsetStat5.inte_stat
      File stat6 = CalcuUpsetStat6.inte_stat
      File stat7 = CalcuUpsetStat7.inte_stat
      File stat8 = CalcuUpsetStat8.inte_stat
  }
 }



task CalcuUpsetStat1 {
    input {
        String sample
        Array[File] vcf_list
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        while read path; do
            echo "Downloading ${path}"
            gsutil cp "${path}" ./
        done < ~{vcf_list}

        Rscript -e '
          
        #!/usr/bin/env Rscript
              

        calcu_benchmark_stat<-function(dv_vs_pav.fp, dv_vs_pav.tp, dv_vs_mc_release.fp, dv_vs_mc_release.tp, dv_vs_mc_liftover.fp, dv_vs_mc_liftover.tp, benchmark_callset_1, benchmark_callset_2, benchmark_callset_3){
            all_variant_sites = unique(rbind(dv_vs_pav.fp[,c(1:3)],dv_vs_pav.tp[,c(1:3)],
                                    dv_vs_mc_release.fp[,c(1:3)] ,dv_vs_mc_release.tp[,c(1:3)],
                                    dv_vs_mc_liftover.fp[,c(1:3)],dv_vs_mc_liftover.tp[,c(1:3)]))

            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_pav.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_1, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_mc_release.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_2, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_mc_liftover.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_3, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = apply(all_variant_sites[,c(4:6)],1,function(x){paste(unique(x[x!="uniq"]),collapse = ",")})
            all_variant_sites[,ncol(all_variant_sites)+1] = sapply(all_variant_sites[,3], function(x){strsplit(as.character(x),"_")[[1]][4]})
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "type"
            all_variant_sites[,ncol(all_variant_sites)+1] = sapply(all_variant_sites[,3], function(x){strsplit(as.character(x),"_")[[1]][5]})
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "length"
            all_variant_sites[,ncol(all_variant_sites)+1] = "SNV"
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "cate"
            all_variant_sites[all_variant_sites$type == "DEL" & all_variant_sites$length<50,]$cate = "DEL_under50bp"
            all_variant_sites[all_variant_sites$type == "INS" & all_variant_sites$length<50,]$cate = "INS_under50bp"
            all_variant_sites[all_variant_sites$type == "DEL" & all_variant_sites$length>49,]$cate = "DEL_over50bp"
            all_variant_sites[all_variant_sites$type == "INS" & all_variant_sites$length>49,]$cate = "INS_over50bp"

            out = table(all_variant_sites[,c(10,7)])
            return(out)
        }

        sample = ~{sample}
        query = paste(sample, "deepvariant.g.alt_cleaned", sep=".")
        ref1 = paste("hgsv_3", sample, sep=".")
        ref2 = paste("hprc_y2_release", sample, sep=".")
        ref3 = paste("hprc_y2_mc_pangenie", sample, sep=".")
        name_query = "deepvariant"
        name1 = "PAV"
        name2 = "MC_release"
        name3 = "MC_liftover"

        dv_vs_pav.fp = read.table(paste(sample, query, "vs", ref1,"fp_query.vcf.gz", sep="."))
        dv_vs_pav.tp = read.table(paste(sample, query, "vs", ref1,"tp_query.vcf.gz", sep="."))

        dv_vs_mc_release.fp = read.table(paste(sample, query,"vs", ref2,"fp_query.vcf.gz", sep="."))
        dv_vs_mc_release.tp = read.table(paste(sample, query,"vs", ref2,"tp_query.vcf.gz", sep="."))

        dv_vs_mc_liftover.fp = read.table(paste(sample, query,"vs", ref3,"fp_query.vcf.gz", sep="."))
        dv_vs_mc_liftover.tp = read.table(paste(sample, query,"vs", ref3,"tp_query.vcf.gz", sep="."))

        dv_stat = calcu_benchmark_stat(dv_vs_pav.fp, dv_vs_pav.tp, dv_vs_mc_release.fp, dv_vs_mc_release.tp, dv_vs_mc_liftover.fp, dv_vs_mc_liftover.tp, name1, name2, name3)
        write.table(dv_stat, paste(sample, paste(name_query,name1, name2, name3, sep=".vs."), "stat", sep="."), quote=F, sep="\t", col.names=T, row.names=T)

        '


    >>>

    output {
        File inte_stat = "~{sample}.deepvariant.vs.PAV.vs.MC_release.vs.MC_liftover.stat"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(vcf_list, "GiB")*2),
        disk_gb: 15 + ceil(size(vcf_list, "GiB")*2),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

task CalcuUpsetStat2 {
    input {
        String sample
        Array[File] vcf_list
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        while read path; do
            echo "Downloading ${path}"
            gsutil cp "${path}" ./
        done < ~{vcf_list}

        Rscript -e '
          
        #!/usr/bin/env Rscript
              

        calcu_benchmark_stat<-function(dv_vs_pav.fp, dv_vs_pav.tp, dv_vs_mc_release.fp, dv_vs_mc_release.tp, dv_vs_mc_liftover.fp, dv_vs_mc_liftover.tp, benchmark_callset_1, benchmark_callset_2, benchmark_callset_3){
            all_variant_sites = unique(rbind(dv_vs_pav.fp[,c(1:3)],dv_vs_pav.tp[,c(1:3)],
                                    dv_vs_mc_release.fp[,c(1:3)] ,dv_vs_mc_release.tp[,c(1:3)],
                                    dv_vs_mc_liftover.fp[,c(1:3)],dv_vs_mc_liftover.tp[,c(1:3)]))

            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_pav.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_1, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_mc_release.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_2, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_mc_liftover.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_3, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = apply(all_variant_sites[,c(4:6)],1,function(x){paste(unique(x[x!="uniq"]),collapse = ",")})
            all_variant_sites[,ncol(all_variant_sites)+1] = sapply(all_variant_sites[,3], function(x){strsplit(as.character(x),"_")[[1]][4]})
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "type"
            all_variant_sites[,ncol(all_variant_sites)+1] = sapply(all_variant_sites[,3], function(x){strsplit(as.character(x),"_")[[1]][5]})
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "length"
            all_variant_sites[,ncol(all_variant_sites)+1] = "SNV"
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "cate"
            all_variant_sites[all_variant_sites$type == "DEL" & all_variant_sites$length<50,]$cate = "DEL_under50bp"
            all_variant_sites[all_variant_sites$type == "INS" & all_variant_sites$length<50,]$cate = "INS_under50bp"
            all_variant_sites[all_variant_sites$type == "DEL" & all_variant_sites$length>49,]$cate = "DEL_over50bp"
            all_variant_sites[all_variant_sites$type == "INS" & all_variant_sites$length>49,]$cate = "INS_over50bp"

            out = table(all_variant_sites[,c(10,7)])
            return(out)
        }

        sample = ~{sample}

        query = paste("hgsv_3", sample, sep=".")
        ref1 = paste(sample, "deepvariant.g.alt_cleaned", sep=".")
        ref2 = paste("hprc_y2_release", sample, sep=".")
        ref3 = paste("hprc_y2_mc_pangenie", sample, sep=".")
        name_query = "PAV"
        name1 = "deep_variant"
        name2 = "MC_release"
        name3 = "MC_liftover"
        pav_vs_dv.fp = read.table(paste(sample, query,"vs", ref1,"fp_ref.vcf.gz", sep="."))
        pav_vs_dv.tp = read.table(paste(sample, query,"vs", ref1,"tp_ref.vcf.gz", sep="."))

        pav_vs_mc_release.fp = read.table(paste(sample, query,"vs", ref2,"fp_ref.vcf.gz", sep="."))
        pav_vs_mc_release.tp = read.table(paste(sample, query,"vs", ref2,"tp_ref.vcf.gz", sep="."))

        pav_vs_mc_liftover.fp = read.table(paste(sample, query,"vs", ref3,"fp_ref.vcf.gz", sep="."))
        pav_vs_mc_liftover.tp = read.table(paste(sample, query,"vs", ref3,"tp_ref.vcf.gz", sep="."))

        pav_stat = calcu_benchmark_stat(pav_vs_dv.fp, pav_vs_dv.tp, 
                                        pav_vs_mc_release.fp, pav_vs_mc_release.tp,
                                        pav_vs_mc_liftover.fp, pav_vs_mc_liftover.tp,
                                                name1, name2, name3)
        write.table(dv_stat, paste(sample, paste(name_query,name1, name2, name3, sep=".vs."), "stat", sep="."), quote=F, sep="\t", col.names=T, row.names=T)

        '


    >>>

    output {
        File inte_stat = "~{sample}.PAV.vs.deep_variant.vs.MC_release.vs.MC_liftover.stat"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(vcf_list, "GiB")*2),
        disk_gb: 15 + ceil(size(vcf_list, "GiB")*2),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

task CalcuUpsetStat3 {
    input {
        String sample
        Array[File] vcf_list
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        while read path; do
            echo "Downloading ${path}"
            gsutil cp "${path}" ./
        done < ~{vcf_list}

        Rscript -e '
          
        #!/usr/bin/env Rscript
              

        calcu_benchmark_stat<-function(dv_vs_pav.fp, dv_vs_pav.tp, dv_vs_mc_release.fp, dv_vs_mc_release.tp, dv_vs_mc_liftover.fp, dv_vs_mc_liftover.tp, benchmark_callset_1, benchmark_callset_2, benchmark_callset_3){
            all_variant_sites = unique(rbind(dv_vs_pav.fp[,c(1:3)],dv_vs_pav.tp[,c(1:3)],
                                    dv_vs_mc_release.fp[,c(1:3)] ,dv_vs_mc_release.tp[,c(1:3)],
                                    dv_vs_mc_liftover.fp[,c(1:3)],dv_vs_mc_liftover.tp[,c(1:3)]))

            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_pav.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_1, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_mc_release.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_2, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_mc_liftover.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_3, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = apply(all_variant_sites[,c(4:6)],1,function(x){paste(unique(x[x!="uniq"]),collapse = ",")})
            all_variant_sites[,ncol(all_variant_sites)+1] = sapply(all_variant_sites[,3], function(x){strsplit(as.character(x),"_")[[1]][4]})
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "type"
            all_variant_sites[,ncol(all_variant_sites)+1] = sapply(all_variant_sites[,3], function(x){strsplit(as.character(x),"_")[[1]][5]})
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "length"
            all_variant_sites[,ncol(all_variant_sites)+1] = "SNV"
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "cate"
            all_variant_sites[all_variant_sites$type == "DEL" & all_variant_sites$length<50,]$cate = "DEL_under50bp"
            all_variant_sites[all_variant_sites$type == "INS" & all_variant_sites$length<50,]$cate = "INS_under50bp"
            all_variant_sites[all_variant_sites$type == "DEL" & all_variant_sites$length>49,]$cate = "DEL_over50bp"
            all_variant_sites[all_variant_sites$type == "INS" & all_variant_sites$length>49,]$cate = "INS_over50bp"

            out = table(all_variant_sites[,c(10,7)])
            return(out)
        }

        sample = ~{sample}

        query = paste("hprc_y2_release", sample, sep=".")
        ref1 = paste(sample, "deepvariant.g.alt_cleaned", sep=".")
        ref2 = paste("hgsv_3", sample, sep=".")
        ref3 = paste("hprc_y2_mc_pangenie", sample, sep=".")
        name_query = "MC_release"
        name1 = "deep_variant"
        name2 = "PAV"
        name3 = "MC_liftover"

        mc_release_vs_dv.fp = read.table(paste(sample, query,"vs", ref1,"fp_ref.vcf.gz", sep="."))
        mc_release_vs_dv.tp = read.table(paste(sample, query,"vs", ref1,"tp_ref.vcf.gz", sep="."))

        mc_release_vs_pav.fp = read.table(paste(sample, query,"vs", ref2,"fp_query.vcf.gz", sep="."))
        mc_release_vs_pav.tp = read.table(paste(sample, query,"vs", ref2,"tp_query.vcf.gz", sep="."))

        mc_release_vs_mc_liftover.fp = read.table(paste(sample, query,"vs", ref3,"fp_query.vcf.gz", sep="."))
        mc_release_vs_mc_liftover.tp = read.table(paste(sample, query,"vs", ref3,"tp_query.vcf.gz", sep="."))

        mc_release_stat = calcu_benchmark_stat(mc_release_vs_dv.fp, mc_release_vs_dv.tp, 
                                        mc_release_vs_pav.fp, mc_release_vs_pav.tp,
                                        mc_release_vs_mc_liftover.fp, mc_release_vs_mc_liftover.tp,
                                                name1, name2, name3)
        write.table(dv_stat, paste(sample, paste(name_query,name1, name2, name3, sep=".vs."), "stat", sep="."), quote=F, sep="\t", col.names=T, row.names=T)


        '


    >>>

    output {
        File inte_stat = "~{sample}.MC_release.vs.deep_variant.vs.PAV.vs.MC_liftover.stat"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(vcf_list, "GiB")*2),
        disk_gb: 15 + ceil(size(vcf_list, "GiB")*2),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

task CalcuUpsetStat4 {
    input {
        String sample
        Array[File] vcf_list
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        while read path; do
            echo "Downloading ${path}"
            gsutil cp "${path}" ./
        done < ~{vcf_list}

        Rscript -e '
          
        #!/usr/bin/env Rscript
              

        calcu_benchmark_stat<-function(dv_vs_pav.fp, dv_vs_pav.tp, dv_vs_mc_release.fp, dv_vs_mc_release.tp, dv_vs_mc_liftover.fp, dv_vs_mc_liftover.tp, benchmark_callset_1, benchmark_callset_2, benchmark_callset_3){
            all_variant_sites = unique(rbind(dv_vs_pav.fp[,c(1:3)],dv_vs_pav.tp[,c(1:3)],
                                    dv_vs_mc_release.fp[,c(1:3)] ,dv_vs_mc_release.tp[,c(1:3)],
                                    dv_vs_mc_liftover.fp[,c(1:3)],dv_vs_mc_liftover.tp[,c(1:3)]))

            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_pav.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_1, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_mc_release.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_2, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_mc_liftover.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_3, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = apply(all_variant_sites[,c(4:6)],1,function(x){paste(unique(x[x!="uniq"]),collapse = ",")})
            all_variant_sites[,ncol(all_variant_sites)+1] = sapply(all_variant_sites[,3], function(x){strsplit(as.character(x),"_")[[1]][4]})
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "type"
            all_variant_sites[,ncol(all_variant_sites)+1] = sapply(all_variant_sites[,3], function(x){strsplit(as.character(x),"_")[[1]][5]})
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "length"
            all_variant_sites[,ncol(all_variant_sites)+1] = "SNV"
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "cate"
            all_variant_sites[all_variant_sites$type == "DEL" & all_variant_sites$length<50,]$cate = "DEL_under50bp"
            all_variant_sites[all_variant_sites$type == "INS" & all_variant_sites$length<50,]$cate = "INS_under50bp"
            all_variant_sites[all_variant_sites$type == "DEL" & all_variant_sites$length>49,]$cate = "DEL_over50bp"
            all_variant_sites[all_variant_sites$type == "INS" & all_variant_sites$length>49,]$cate = "INS_over50bp"

            out = table(all_variant_sites[,c(10,7)])
            return(out)
        }

        sample = ~{sample}


        query = paste("hprc_y2_mc_pangenie", sample, sep=".")
        ref1 = paste(sample, "deepvariant.g.alt_cleaned", sep=".")
        ref2 = paste("hgsv_3", sample, sep=".")
        ref3 = paste("hprc_y2_release", sample, sep=".")
        name_query = "MC_liftover"
        name1 = "deep_variant"
        name2 = "PAV"
        name3 = "MC_release"
        mc_liftover_vs_dv.fp = read.table(paste(sample, query,"vs", ref1,"fp_ref.vcf.gz", sep="."))
        mc_liftover_vs_dv.tp = read.table(paste(sample, query,"vs", ref1,"tp_ref.vcf.gz", sep="."))

        mc_liftover_vs_pav.fp = read.table(paste(sample, query,"vs", ref2,"fp_query.vcf.gz", sep="."))
        mc_liftover_vs_pav.tp = read.table(paste(sample, query,"vs", ref2,"tp_query.vcf.gz", sep="."))

        mc_liftover_vs_mc_release.fp = read.table(paste(sample, query,"vs", ref3,"fp_ref.vcf.gz", sep="."))
        mc_liftover_vs_mc_release.tp = read.table(paste(sample, query,"vs", ref3,"tp_ref.vcf.gz", sep="."))

        mc_liftover_stat = calcu_benchmark_stat(mc_liftover_vs_dv.fp, mc_liftover_vs_dv.tp, 
                                                mc_liftover_vs_pav.fp, mc_liftover_vs_pav.tp,
                                                mc_liftover_vs_mc_release.fp, mc_liftover_vs_mc_release.tp,
                                                name1, name2, name3)
        write.table(dv_stat, paste(sample, paste(name_query,name1, name2, name3, sep=".vs."), "stat", sep="."), quote=F, sep="\t", col.names=T, row.names=T)


        '


    >>>

    output {
        File inte_stat = "~{sample}.MC_liftover.vs.deep_variant.vs.PAV.vs.MC_release.stat"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(vcf_list, "GiB")*2),
        disk_gb: 15 + ceil(size(vcf_list, "GiB")*2),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

task CalcuUpsetStat5 {
    input {
        String sample
        Array[File] vcf_list
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        while read path; do
            echo "Downloading ${path}"
            gsutil cp "${path}" ./
        done < ~{vcf_list}

        Rscript -e '
          
        #!/usr/bin/env Rscript
              

        calcu_benchmark_stat<-function(dv_vs_pav.fp, dv_vs_pav.tp, dv_vs_mc_release.fp, dv_vs_mc_release.tp, dv_vs_mc_liftover.fp, dv_vs_mc_liftover.tp, benchmark_callset_1, benchmark_callset_2, benchmark_callset_3){
            all_variant_sites = unique(rbind(dv_vs_pav.fp[,c(1:3)],dv_vs_pav.tp[,c(1:3)],
                                    dv_vs_mc_release.fp[,c(1:3)] ,dv_vs_mc_release.tp[,c(1:3)],
                                    dv_vs_mc_liftover.fp[,c(1:3)],dv_vs_mc_liftover.tp[,c(1:3)]))

            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_pav.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_1, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_mc_release.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_2, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_mc_liftover.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_3, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = apply(all_variant_sites[,c(4:6)],1,function(x){paste(unique(x[x!="uniq"]),collapse = ",")})
            all_variant_sites[,ncol(all_variant_sites)+1] = sapply(all_variant_sites[,3], function(x){strsplit(as.character(x),"_")[[1]][4]})
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "type"
            all_variant_sites[,ncol(all_variant_sites)+1] = sapply(all_variant_sites[,3], function(x){strsplit(as.character(x),"_")[[1]][5]})
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "length"
            all_variant_sites[,ncol(all_variant_sites)+1] = "SNV"
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "cate"
            all_variant_sites[all_variant_sites$type == "DEL" & all_variant_sites$length<50,]$cate = "DEL_under50bp"
            all_variant_sites[all_variant_sites$type == "INS" & all_variant_sites$length<50,]$cate = "INS_under50bp"
            all_variant_sites[all_variant_sites$type == "DEL" & all_variant_sites$length>49,]$cate = "DEL_over50bp"
            all_variant_sites[all_variant_sites$type == "INS" & all_variant_sites$length>49,]$cate = "INS_over50bp"

            out = table(all_variant_sites[,c(10,7)])
            return(out)
        }

        sample = ~{sample}

        query = paste("Jiadong_integration", sample, sep=".")
        ref1 = paste("hgsv_3", sample, sep=".")
        ref2 = paste("hprc_y2_release", sample, sep=".")
        ref3 = paste("hprc_y2_mc_pangenie", sample, sep=".")
        name_query = "Jiadong_integration"
        name1 = "PAV"
        name2 = "MC_release"
        name3 = "MC_liftover"

        jiadong_vs_pav.fp = read.table(paste(sample, query,"vs", ref1,"fp_query.vcf.gz", sep="."))
        jiadong_vs_pav.tp = read.table(paste(sample, query,"vs", ref1,"tp_query.vcf.gz", sep="."))

        jiadong_vs_mc_release.fp = read.table(paste(sample, query,"vs", ref2,"fp_query.vcf.gz", sep="."))
        jiadong_vs_mc_release.tp = read.table(paste(sample, query,"vs", ref2,"tp_query.vcf.gz", sep="."))

        jiadong_vs_mc_liftover.fp = read.table(paste(sample, query,"vs", ref3,"fp_query.vcf.gz", sep="."))
        jiadong_vs_mc_liftover.tp = read.table(paste(sample, query,"vs", ref3,"tp_query.vcf.gz", sep="."))

        jiadong_stat = calcu_benchmark_stat(jiadong_vs_pav.fp, jiadong_vs_pav.tp, 
                                       jiadong_vs_mc_release.fp, jiadong_vs_mc_release.tp, 
                                       jiadong_vs_mc_liftover.fp, jiadong_vs_mc_liftover.tp, 
                                                name1, name2, name3)
        write.table(dv_stat, paste(sample, paste(name_query,name1, name2, name3, sep=".vs."), "stat", sep="."), quote=F, sep="\t", col.names=T, row.names=T)


        '


    >>>

    output {
        File inte_stat = "~{sample}.Jiadong_integration.vs.PAV.vs.MC_release.vs.MC_liftover.stat"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(vcf_list, "GiB")*2),
        disk_gb: 15 + ceil(size(vcf_list, "GiB")*2),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

task CalcuUpsetStat6 {
    input {
        String sample
        Array[File] vcf_list
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        while read path; do
            echo "Downloading ${path}"
            gsutil cp "${path}" ./
        done < ~{vcf_list}

        Rscript -e '
          
        #!/usr/bin/env Rscript
              

        calcu_benchmark_stat<-function(dv_vs_pav.fp, dv_vs_pav.tp, dv_vs_mc_release.fp, dv_vs_mc_release.tp, dv_vs_mc_liftover.fp, dv_vs_mc_liftover.tp, benchmark_callset_1, benchmark_callset_2, benchmark_callset_3){
            all_variant_sites = unique(rbind(dv_vs_pav.fp[,c(1:3)],dv_vs_pav.tp[,c(1:3)],
                                    dv_vs_mc_release.fp[,c(1:3)] ,dv_vs_mc_release.tp[,c(1:3)],
                                    dv_vs_mc_liftover.fp[,c(1:3)],dv_vs_mc_liftover.tp[,c(1:3)]))

            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_pav.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_1, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_mc_release.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_2, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_mc_liftover.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_3, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = apply(all_variant_sites[,c(4:6)],1,function(x){paste(unique(x[x!="uniq"]),collapse = ",")})
            all_variant_sites[,ncol(all_variant_sites)+1] = sapply(all_variant_sites[,3], function(x){strsplit(as.character(x),"_")[[1]][4]})
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "type"
            all_variant_sites[,ncol(all_variant_sites)+1] = sapply(all_variant_sites[,3], function(x){strsplit(as.character(x),"_")[[1]][5]})
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "length"
            all_variant_sites[,ncol(all_variant_sites)+1] = "SNV"
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "cate"
            all_variant_sites[all_variant_sites$type == "DEL" & all_variant_sites$length<50,]$cate = "DEL_under50bp"
            all_variant_sites[all_variant_sites$type == "INS" & all_variant_sites$length<50,]$cate = "INS_under50bp"
            all_variant_sites[all_variant_sites$type == "DEL" & all_variant_sites$length>49,]$cate = "DEL_over50bp"
            all_variant_sites[all_variant_sites$type == "INS" & all_variant_sites$length>49,]$cate = "INS_over50bp"

            out = table(all_variant_sites[,c(10,7)])
            return(out)
        }

        sample = ~{sample}

        query = paste("hgsv_3", sample, sep=".")
        ref1 = paste("Jiadong_integration", sample, sep=".")
        ref2 = paste("hprc_y2_release", sample, sep=".")
        ref3 = paste("hprc_y2_mc_pangenie", sample, sep=".")
        name_query = "PAV"
        name1 = "Jiadong_integration"
        name2 = "MC_release"
        name3 = "MC_liftover"
        pav_vs_jiadong.fp = read.table(paste(sample, query,"vs", ref1,"fp_ref.vcf.gz", sep="."))
        pav_vs_jiadong.tp = read.table(paste(sample, query,"vs", ref1,"tp_ref.vcf.gz", sep="."))

        pav_vs_mc_release.fp = read.table(paste(sample, query,"vs", ref2,"fp_ref.vcf.gz", sep="."))
        pav_vs_mc_release.tp = read.table(paste(sample, query,"vs", ref2,"tp_ref.vcf.gz", sep="."))

        pav_vs_mc_liftover.fp = read.table(paste(sample, query,"vs", ref3,"fp_ref.vcf.gz", sep="."))
        pav_vs_mc_liftover.tp = read.table(paste(sample, query,"vs", ref3,"tp_ref.vcf.gz", sep="."))

        pav_stat = calcu_benchmark_stat(pav_vs_jiadong.fp, pav_vs_jiadong.tp, 
                                        pav_vs_mc_release.fp, pav_vs_mc_release.tp,
                                        pav_vs_mc_liftover.fp, pav_vs_mc_liftover.tp,
                                                name1, name2, name3)
        write.table(dv_stat, paste(sample, paste(name_query,name1, name2, name3, sep=".vs."), "stat", sep="."), quote=F, sep="\t", col.names=T, row.names=T)


        '


    >>>

    output {
        File inte_stat = "~{sample}.PAV.vs.Jiadong_integration.vs.MC_release.vs.MC_liftover.stat"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(vcf_list, "GiB")*2),
        disk_gb: 15 + ceil(size(vcf_list, "GiB")*2),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

task CalcuUpsetStat7 {
    input {
        String sample
        Array[File] vcf_list
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        while read path; do
            echo "Downloading ${path}"
            gsutil cp "${path}" ./
        done < ~{vcf_list}

        Rscript -e '
          
        #!/usr/bin/env Rscript
              

        calcu_benchmark_stat<-function(dv_vs_pav.fp, dv_vs_pav.tp, dv_vs_mc_release.fp, dv_vs_mc_release.tp, dv_vs_mc_liftover.fp, dv_vs_mc_liftover.tp, benchmark_callset_1, benchmark_callset_2, benchmark_callset_3){
            all_variant_sites = unique(rbind(dv_vs_pav.fp[,c(1:3)],dv_vs_pav.tp[,c(1:3)],
                                    dv_vs_mc_release.fp[,c(1:3)] ,dv_vs_mc_release.tp[,c(1:3)],
                                    dv_vs_mc_liftover.fp[,c(1:3)],dv_vs_mc_liftover.tp[,c(1:3)]))

            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_pav.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_1, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_mc_release.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_2, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_mc_liftover.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_3, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = apply(all_variant_sites[,c(4:6)],1,function(x){paste(unique(x[x!="uniq"]),collapse = ",")})
            all_variant_sites[,ncol(all_variant_sites)+1] = sapply(all_variant_sites[,3], function(x){strsplit(as.character(x),"_")[[1]][4]})
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "type"
            all_variant_sites[,ncol(all_variant_sites)+1] = sapply(all_variant_sites[,3], function(x){strsplit(as.character(x),"_")[[1]][5]})
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "length"
            all_variant_sites[,ncol(all_variant_sites)+1] = "SNV"
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "cate"
            all_variant_sites[all_variant_sites$type == "DEL" & all_variant_sites$length<50,]$cate = "DEL_under50bp"
            all_variant_sites[all_variant_sites$type == "INS" & all_variant_sites$length<50,]$cate = "INS_under50bp"
            all_variant_sites[all_variant_sites$type == "DEL" & all_variant_sites$length>49,]$cate = "DEL_over50bp"
            all_variant_sites[all_variant_sites$type == "INS" & all_variant_sites$length>49,]$cate = "INS_over50bp"

            out = table(all_variant_sites[,c(10,7)])
            return(out)
        }

        sample = ~{sample}

        query = paste("hprc_y2_release", sample, sep=".")
        ref1 = paste("Jiadong_integration", sample, sep=".")
        ref2 = paste("hgsv_3", sample, sep=".")
        ref3 = paste("hprc_y2_mc_pangenie", sample, sep=".")
        name_query = "MC_release"
        name1 = "Jiadong_integration"
        name2 = "PAV"
        name3 = "MC_liftover"

        mc_release_vs_jiadong.fp = read.table(paste(sample, query,"vs", ref1,"fp_ref.vcf.gz", sep="."))
        mc_release_vs_jiadong.tp = read.table(paste(sample, query,"vs", ref1,"tp_ref.vcf.gz", sep="."))

        mc_release_vs_pav.fp = read.table(paste(sample, query,"vs", ref2,"fp_query.vcf.gz", sep="."))
        mc_release_vs_pav.tp = read.table(paste(sample, query,"vs", ref2,"tp_query.vcf.gz", sep="."))

        mc_release_vs_mc_liftover.fp = read.table(paste(sample, query,"vs", ref3,"fp_query.vcf.gz", sep="."))
        mc_release_vs_mc_liftover.tp = read.table(paste(sample, query,"vs", ref3,"tp_query.vcf.gz", sep="."))

        mc_release_stat = calcu_benchmark_stat(mc_release_vs_jiadong.fp, mc_release_vs_jiadong.tp, 
                                        mc_release_vs_pav.fp, mc_release_vs_pav.tp,
                                        mc_release_vs_mc_liftover.fp, mc_release_vs_mc_liftover.tp,
                                                name1, name2, name3)
        write.table(dv_stat, paste(sample, paste(name_query,name1, name2, name3, sep=".vs."), "stat", sep="."), quote=F, sep="\t", col.names=T, row.names=T)

        '
    >>>

    output {
        File inte_stat = "~{sample}.MC_release.vs.Jiadong_integration.vs.PAV.vs.MC_liftover.stat"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(vcf_list, "GiB")*2),
        disk_gb: 15 + ceil(size(vcf_list, "GiB")*2),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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

task CalcuUpsetStat8 {
    input {
        String sample
        Array[File] vcf_list
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        while read path; do
            echo "Downloading ${path}"
            gsutil cp "${path}" ./
        done < ~{vcf_list}

        Rscript -e '
          
        #!/usr/bin/env Rscript
              

        calcu_benchmark_stat<-function(dv_vs_pav.fp, dv_vs_pav.tp, dv_vs_mc_release.fp, dv_vs_mc_release.tp, dv_vs_mc_liftover.fp, dv_vs_mc_liftover.tp, benchmark_callset_1, benchmark_callset_2, benchmark_callset_3){
            all_variant_sites = unique(rbind(dv_vs_pav.fp[,c(1:3)],dv_vs_pav.tp[,c(1:3)],
                                    dv_vs_mc_release.fp[,c(1:3)] ,dv_vs_mc_release.tp[,c(1:3)],
                                    dv_vs_mc_liftover.fp[,c(1:3)],dv_vs_mc_liftover.tp[,c(1:3)]))

            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_pav.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_1, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_mc_release.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_2, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = "uniq"
            all_variant_sites[all_variant_sites[,3]%in%dv_vs_mc_liftover.tp[,3],][,ncol(all_variant_sites)] = paste("vs", benchmark_callset_3, sep="_")
            all_variant_sites[,ncol(all_variant_sites)+1] = apply(all_variant_sites[,c(4:6)],1,function(x){paste(unique(x[x!="uniq"]),collapse = ",")})
            all_variant_sites[,ncol(all_variant_sites)+1] = sapply(all_variant_sites[,3], function(x){strsplit(as.character(x),"_")[[1]][4]})
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "type"
            all_variant_sites[,ncol(all_variant_sites)+1] = sapply(all_variant_sites[,3], function(x){strsplit(as.character(x),"_")[[1]][5]})
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "length"
            all_variant_sites[,ncol(all_variant_sites)+1] = "SNV"
            colnames(all_variant_sites)[ncol(all_variant_sites)] = "cate"
            all_variant_sites[all_variant_sites$type == "DEL" & all_variant_sites$length<50,]$cate = "DEL_under50bp"
            all_variant_sites[all_variant_sites$type == "INS" & all_variant_sites$length<50,]$cate = "INS_under50bp"
            all_variant_sites[all_variant_sites$type == "DEL" & all_variant_sites$length>49,]$cate = "DEL_over50bp"
            all_variant_sites[all_variant_sites$type == "INS" & all_variant_sites$length>49,]$cate = "INS_over50bp"

            out = table(all_variant_sites[,c(10,7)])
            return(out)
        }

        sample = ~{sample}

        query = paste("hprc_y2_mc_pangenie", sample, sep=".")
        ref1 = paste("Jiadong_integration", sample, sep=".")
        ref2 = paste("hgsv_3", sample, sep=".")
        ref3 = paste("hprc_y2_release", sample, sep=".")
        name_query = "MC_liftover"
        name1 = "Jiadong_integration"
        name2 = "PAV"
        name3 = "MC_release"
        mc_liftover_vs_jiadong.fp = read.table(paste(sample, query,"vs", ref1,"fp_ref.vcf.gz", sep="."))
        mc_liftover_vs_jiadong.tp = read.table(paste(sample, query,"vs", ref1,"tp_ref.vcf.gz", sep="."))

        mc_liftover_vs_pav.fp = read.table(paste(sample, query,"vs", ref2,"fp_query.vcf.gz", sep="."))
        mc_liftover_vs_pav.tp = read.table(paste(sample, query,"vs", ref2,"tp_query.vcf.gz", sep="."))

        mc_liftover_vs_mc_release.fp = read.table(paste(sample, query,"vs", ref3,"fp_ref.vcf.gz", sep="."))
        mc_liftover_vs_mc_release.tp = read.table(paste(sample, query,"vs", ref3,"tp_ref.vcf.gz", sep="."))

        mc_liftover_stat = calcu_benchmark_stat(mc_liftover_vs_jiadong.fp, mc_liftover_vs_jiadong.tp, 
                                                mc_liftover_vs_pav.fp, mc_liftover_vs_pav.tp,
                                                mc_liftover_vs_mc_release.fp, mc_liftover_vs_mc_release.tp,
                                                name1, name2, name3)
        write.table(dv_stat, paste(sample, paste(name_query,name1, name2, name3, sep=".vs."), "stat", sep="."), quote=F, sep="\t", col.names=T, row.names=T)

        '
    >>>

    output {
        File inte_stat = "~{sample}.MC_liftover.vs.Jiadong_integration.vs.PAV.vs.MC_release.stat"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(vcf_list, "GiB")*2),
        disk_gb: 15 + ceil(size(vcf_list, "GiB")*2),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

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
















