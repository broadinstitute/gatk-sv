version 1.0

import "Structs.wdl"
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks

workflow PlotSaturation {
  input {
    File vcf_file
    File vcf_idx
    File sample_list
    Array[String] contig_list

    File anno_script_bash
    File anno_script_helper_R
    File repeat_mask
    File simple_repeats
    File segmental_duplicates

    String sv_base_mini_docker
    String sv_pipeline_base_docker
    File? monitoring_script

    RuntimeAttr? runtime_attr_extract_chrom_vcf
    RuntimeAttr? runtime_attr_split_variants_by_size
    RuntimeAttr? runtime_attr_extract_snv_sites
    RuntimeAttr? runtime_attr_extract_sm_indel_sites
    RuntimeAttr? runtime_attr_extract_lg_indel_sites
    RuntimeAttr? runtime_attr_extract_sv_sites
    RuntimeAttr? runtime_attr_annotate_genomic_context_snv
    RuntimeAttr? runtime_attr_annotate_genomic_context_sm_indel
    RuntimeAttr? runtime_attr_annotate_genomic_context_lg_indel
    RuntimeAttr? runtime_attr_annotate_genomic_context_sv
    RuntimeAttr? runtime_attr_calcu_satu_table_snv
    RuntimeAttr? runtime_attr_calcu_satu_table_sm_indel
    RuntimeAttr? runtime_attr_calcu_satu_table_lg_indel
    RuntimeAttr? runtime_attr_calcu_satu_table_sv

  }

  scatter(contig in contig_list){
    call  LongReadGenotypeTasks.ExtractChromosomeVcf {
      input:
        input_vcf = vcf_file,
        input_vcf_idx = vcf_idx,
        chromosome = contig,
        docker_image = sv_pipeline_base_docker,
        monitoring_script = monitoring_script,
        runtime_attr_override = runtime_attr_extract_chrom_vcf
    }

    call LongReadGenotypeTasks.SplitVariantsBySize{
      input:
        input_vcf = ExtractChromosomeVcf.output_vcf,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_split_variants_by_size
    }

    cal LongReadGenotypeTasks.ExtractVariantSites as extract_snv_sites{
      input:
        input_vcf = SplitVariantsBySize.snv_vcf,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_extract_snv_sites
    }

    call LongReadGenotypeTasks.AnnotateGenomicContext as annotate_genomic_context_snv{
      input:
        variant_sites = extract_snv_sites.variant_sites,
        anno_script_bash = anno_script_bash,
        anno_script_Rscript = anno_script_helper_R,
        repeat_mask = repeat_mask,
        simple_repeats = simple_repeats,
        segmental_duplicates = segmental_duplicates,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_annotate_genomic_context_snv
    }

    call LongReadGenotypeTasks.ExtractVariantSites as extract_sm_indel_sites{
      input:
        input_vcf = SplitVariantsBySize.indel_1_30_vcf,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_extract_sm_indel_sites
    }

    call LongReadGenotypeTasks.AnnotateGenomicContext as annotate_genomic_context_sm_indel{
      input:
        variant_sites = extract_sm_indel_sites.variant_sites,
        anno_script_bash = anno_script_bash,
        anno_script_Rscript = anno_script_helper_R,
        repeat_mask = repeat_mask,
        simple_repeats = simple_repeats,
        segmental_duplicates = segmental_duplicates,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_annotate_genomic_context_sm_indel
    }

    call LongReadGenotypeTasks.ExtractVariantSites as extract_lg_indel_sites{
      input:
        input_vcf = SplitVariantsBySize.indel_31_50_vcf,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_extract_lg_indel_sites
    }

    call LongReadGenotypeTasks.AnnotateGenomicContext as annotate_genomic_context_lg_indel{
      input:
        variant_sites = extract_lg_indel_sites.variant_sites,
        anno_script_bash = anno_script_bash,
        anno_script_Rscript = anno_script_helper_R,
        repeat_mask = repeat_mask,
        simple_repeats = simple_repeats,
        segmental_duplicates = segmental_duplicates,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_annotate_genomic_context_lg_indel
    }

    call LongReadGenotypeTasks.ExtractVariantSites as extract_sv_sites{
      input:
        input_vcf = SplitVariantsBySize.sv_vcf,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_extract_sv_sites
    }

    call LongReadGenotypeTasks.AnnotateGenomicContext as annotate_genomic_context_sv{
      input:
        variant_sites = extract_sv_sites.variant_sites,
        anno_script_bash = anno_script_bash,
        anno_script_Rscript = anno_script_helper_R,
        repeat_mask = repeat_mask,
        simple_repeats = simple_repeats,
        segmental_duplicates = segmental_duplicates,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_annotate_genomic_context_sv
    }

    call CalcuSaturationTable as calcu_satu_table_snv{
      input:
        vcf_file   = extract_snv_sites.updated_vcf,
        SVID_GC = annotate_genomic_context_snv.variant_anno,
        sample_list = sample_list,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_calcu_satu_table_snv
    }

    call CalcuSaturationTable as calcu_satu_table_sm_indel{
      input:
        vcf_file   = extract_sm_indel_sites.updated_vcf,
        SVID_GC = annotate_genomic_context_sm_indel.variant_anno,
        sample_list = sample_list,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_calcu_satu_table_sm_indel
    }

    call CalcuSaturationTable as calcu_satu_table_lg_indel{
      input:
        vcf_file   = extract_lg_indel_sites.updated_vcf,
        SVID_GC = annotate_genomic_context_lg_indel.variant_anno,
        sample_list = sample_list,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_calcu_satu_table_lg_indel
    }

    call CalcuSaturationTable as calcu_satu_table_sv{
      input:
        vcf_file   = extract_sv_sites.updated_vcf,
        SVID_GC = annotate_genomic_context_sv.variant_anno,
        sample_list = sample_list,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_calcu_satu_table_sv
    }

    call IntegrateSaturateTables as inte_satu_tables{
      input:
        snv_table = calcu_satu_table_snv.Satu_Stat,
        sm_indel_table = calcu_satu_table_sm_indel.Satu_Stat,
        lg_indel_table = calcu_satu_table_lg_indel.Satu_Stat,
        sv_table = calcu_satu_table_sv.Satu_Stat,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_inte_satu_tables
    }
    }

    output{
      File Array[satu_stat] = inte_satu_tables.integrated_stat
    }
  }



task CalcuSaturationTable {
  input {
    File vcf_file         # Input VCF file (bgzipped or plain)
    File SVID_GC
    File sample_list
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(vcf_file,'.vcf.gz')

  command <<<
    set -euxo pipefail


    #use python script to extract IDs of carriers for each variants
    python3 <<CODE

    import sys
    import pysam

    def vcf_to_bed(input_vcf, output_bed):
        vcf = pysam.VariantFile(input_vcf, "r")
        with open(output_bed, "w") as out:
            # Write header
            header_fields = ["#CHROM", "START", "END", "SVID", "REF", "ALT", "AC", "AF", "AN", "FILTER", "ALT_CARRIERS"]
            out.write("\t".join(header_fields) + "\n")
            for record in vcf:
                chrom = record.chrom
                start = record.pos - 1  # BED is 0-based
                end = start + len(record.ref)
                id = record.id
                ref = record.ref
                alts = record.alts if record.alts else ["."]
                alt = ",".join(alts)
                ac = ",".join(str(x) for x in record.info.get("AC", ["."]))
                af = ",".join(f"{x:.5f}" for x in record.info.get("AF", [0.0]))
                an = record.info.get("AN", ".")
                filt = ";".join(record.filter.keys()) if record.filter else "PASS"
                # Identify samples with ALT allele (i.e., genotype includes non-zero)
                alt_carriers = []
                for sample in vcf.header.samples:
                    gt = record.samples[sample].get("GT")
                    if gt and any(allele not in (None, 0) for allele in gt):
                        alt_carriers.append(sample)
                alt_carrier_str = ",".join(alt_carriers) if alt_carriers else "."
                out.write("\t".join(map(str, [
                    chrom, start, end, id, ref, alt, ac, af, an, filt, alt_carrier_str
                ])) + "\n")

    vcf_to_bed("~{vcf_file}", "~{prefix}.bed")

    CODE
    


    # use R script to calculate the saturation table
    Rscript -e '

    calcu_satu_table<-function(dat){
        samples=read.table("~{sample_list}")
        samples[,2] = 0
        samples[,3] = 0

        #calculate the count of new variants carried by each individual
        for(i in 1:nrow(samples)){
          dat[,col_count+1] = sapply(dat$ALT_CARRIERS, function(x){length(strsplit(x,",")[[1]][strsplit(x,",")[[1]]%in%c(samples[i,1])])})
          samples[i,2] = nrow(dat[dat[,col_count+1] == 1, ])
          samples[i,3] = sum(dat[1:i, 2])
          dat = dat[dat[,col_count+1] == 0,]
        }
        colnames(samples)  = c("Samp","New_Variants","Total_Variants")

        return(samples)
    }

    input_bed = "~{prefix}.bed"

    dat=read.table(input_bed, sep="\t", comment.char="", header=T)
    dat = dat[dat$FILTER=="PASS",]
    col_count = ncol(dat)

    svid_gc = read.table("~{SVID_GC}")
    colnames(svid_gc) = c("SVID", "Genomic_Context")

    satu_stat_us = calcu_satu_table(dat[dat$SVID%in%svid_gc[svid_gc[,2]=="US",][,1],])
    satu_stat_rm = calcu_satu_table(dat[dat$SVID%in%svid_gc[svid_gc[,2]=="RM",][,1],])
    satu_stat_sd = calcu_satu_table(dat[dat$SVID%in%svid_gc[svid_gc[,2]=="SD",][,1],])
    satu_stat_sr = calcu_satu_table(dat[dat$SVID%in%svid_gc[svid_gc[,2]=="SR",][,1],])

    satu_stat = merge(satu_stat_us[,c(1,3)], satu_stat_rm[,c(1,3)],  by = "Samp") 
    satu_stat = merge(satu_stat, satu_stat_sd[,c(1,3)],  by = "Samp") 
    satu_stat = merge(satu_stat, satu_stat_sr[,c(1,3)],  by = "Samp") 
    colnames(satu_stat) = c("samp","US","RM","SD","SR")

    write.table(satu_stat, "~{prefix}.satu.stat", quote=F, sep="\t", col.names=T, row.names=F)
    
    '

  >>>

  output {
    File Satu_Stat = "~{prefix}.satu.stat"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 20 + ceil(size(vcf_file,"GiB")*2),
    disk_gb: 30 + ceil(size(vcf_file,"GiB")*3),
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

task IntegrateSaturateTables {
  input {
    File snv_table
    File sm_indel_table
    File lg_indel_table
    File sv_table
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(snv_table,'.satu.stat')

  command <<<
    set -euxo pipefail


    # use R script to calculate the saturation table
    Rscript -e '

    snv_table = read.table("~{snv_table}", header = T)
    colnames(snv_table)[c(2:ncol(snv_table))] = paste("snv", colnames(snv_table)[c(2:ncol(snv_table))], sep="_") 

    sm_indel_table = read.table("~{sm_indel_table}", header = T)
    colnames(sm_indel_table)[c(2:ncol(sm_indel_table))] = paste("sm_indel", colnames(sm_indel_table)[c(2:ncol(sm_indel_table))], sep="_") 

    lg_indel_table = read.table("~{lg_indel_table}", header = T)
    colnames(lg_indel_table)[c(2:ncol(lg_indel_table))] = paste("lg_indel", colnames(lg_indel_table)[c(2:ncol(lg_indel_table))], sep="_") 

    sv_table = read.table("~{sv_table}", header = T)
    colnames(sv_table)[c(2:ncol(sv_table))] = paste("sv", colnames(sv_table)[c(2:ncol(sv_table))], sep="_") 


    out = merge(snv_table, sm_indel_table, by="samp")
    out = merge(out, lg_indel_table, by="samp")
    out = merge(out, sv_table, by="samp")

    write.table(out, "~{prefix}.integrated.satu.stat", quote=F, sep="\t", col.names=T, row.names=F)
    
    '

  >>>

  output {
    File integrated_stat = "~{prefix}.integrated.satu.stat"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 20 + ceil(size(vcf_file,"GiB")*2),
    disk_gb: 30 + ceil(size(vcf_file,"GiB")*3),
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

























