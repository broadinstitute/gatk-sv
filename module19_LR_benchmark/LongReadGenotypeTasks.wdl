version 1.0

import "Structs.wdl"

task AddDummyGT {
  input {
    File sites_file         # Input VCF file (bgzipped or plain)
    File sites_idx
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(sites_file,'.sites.gz')

  command <<<
    set -euxo pipefail

    #take out vcf header
    bcftools view -h ~{sites_file} > header.tmp
    #revise to add sample column
    head -n $(($(wc -l < header.tmp) - 1)) header.tmp > header.vcf
    tail -n 1 header.tmp | sed -e 's/$/\tFORMAT\tDUMMY/' >> header.vcf
    bcftools view -H ~{sites_file} | sed -e 's/$/\tGT\t0|1/' >> header.vcf
    bgzip header.vcf 
    mv header.vcf.gz "~{prefix}.with_dummy_gt.vcf.gz"
    tabix -p vcf  "~{prefix}.with_dummy_gt.vcf.gz"

  >>>

  output {
    File vcf_file = "~{prefix}.with_dummy_gt.vcf.gz"
    File vcf_idx = "~{prefix}.with_dummy_gt.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(sites_file, "GiB")*2),
    disk_gb: 15 + ceil(size(sites_file, "GiB")*2),
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

task AnnotateGenomicContext {
  input {
    File variant_sites
    File anno_script_Rscript
    File anno_script_bash
    File repeat_mask
    File segmental_duplicates
    File simple_repeats
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(variant_sites, '.bed')
  command <<<
    set -euxo pipefail

    bash ~{anno_script_bash} ~{variant_sites} ~{prefix}.SVID_GC \
    --rm ~{repeat_mask} \
    --sd ~{segmental_duplicates} \
    --sr ~{simple_repeats}

  >>>

  output {
    File variant_anno = "~{prefix}.SVID_GC"
  }


  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: ceil(size(variant_sites, "GiB"))*10,
    disk_gb: 10 + ceil(size(variant_sites, "GiB"))*10,
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

task AddGenomicContextToVcfPython {
  input {
    File vcf_file         # Input VCF file (bgzipped or plain)
    File svid_annotation  # 2-column TSV: SVID <tab> annotation
    String docker_image
    RuntimeAttr? runtime_attr_override
  }


  String prefix = basename(vcf_file,'.vcf.gz')

  command <<<

    python3 <<CODE
    import pysam
    import os
    from collections import defaultdict

    def determine_svtype(ref, alt):
        if len(ref) == 1 and len(alt) == 1:
            return "SNV"
        elif len(ref) < len(alt):
            return "INS"
        elif len(ref) > len(alt):
            return "DEL"
        else:
            return "OTH"

    def determine_svlen(ref, alt):
        return str(abs(len(alt) - len(ref)))

    # Load SVID -> annotation mapping
    svid_anno = {}

    print('reading in SVID and their genomic context ... ')
    with open("~{svid_annotation}", "r") as f:
        for line in f:
            if line.strip():
                svid, annotation = line.strip().split('\t')
                svid_anno[svid] = annotation

    print('complete')

    # Prepare output VCF writers per annotation
    vcf_in = pysam.VariantFile("~{vcf_file}", "r")

    header = vcf_in.header
    header.info.add("GC", number=1, type="String", description="Genomic context of the variant")

    vcf_out = pysam.VariantFile("~{prefix}.GC_anno.vcf.gz", 'w', header = header)

    for rec in vcf_in.fetch():
        svid = rec.id
        rec_new = rec
        rec_new.info["GC"] = svid_anno[svid]
        vcf_out.write(rec_new)

    vcf_out.close()

    CODE
  >>>

  output {
    File annotated_vcf = "~{prefix}.GC_anno.vcf.gz"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: 15 + ceil(size(vcf_file, "GiB")*2),
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

task AddGenomicContextToVcfR {
  input {
    File vcf_file         # Input VCF file (bgzipped or plain)
    File svid_annotation  # 2-column TSV: SVID <tab> annotation
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(vcf_file,'.vcf.gz')

  command <<<
    set -euxo pipefail

    # use R script to add GC to the vcf
    Rscript -e '


    read_or_empty <- function(file_path) {

      lines <- readLines(file_path, warn = FALSE)
      if(any(!grepl("^#", lines))){
        return(read.table(file_path, header = FALSE, sep = "", stringsAsFactors = FALSE))
      } else{
        empty_df <- data.frame(matrix(ncol = 10, nrow = 0))
        colnames(empty_df) <- paste0("V", 1:10)
        return(empty_df)

      }
    }


    svid_gc <- read.table("~{svid_annotation}", header = TRUE)
    vcf_in <- read_or_empty("~{vcf_file}")
    colnames(vcf_in)[3] = "SVID"
    svid_gc = unique(svid_gc)
    vcf_out <- merge(vcf_in, svid_gc, by="SVID")
    vcf_out[,8] = apply(vcf_out[,c("V8","GC")], 1, function(x){paste(x[1], paste("GC", x[2], sep="="), sep=";")})
    vcf_out_v2 = vcf_out[,c(2,3,1,4:(ncol(vcf_out)-1))]

    write.table(vcf_out_v2, file = "~{prefix}.GC_anno.vcf", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  
    '


    #use python script to add GC to vcf header
    python3 <<CODE

    import pysam
    import os

    vcf_in = pysam.VariantFile("~{vcf_file}", "r")
    header = vcf_in.header
    header.info.add("GC", number=1, type="String", description="Genomic context of the variant")
    vcf_out = pysam.VariantFile("~{prefix}.header.vcf.gz", 'w', header = header)
    vcf_in.close()
    vcf_out.close()

    CODE
    

    cat <(zcat ~{prefix}.header.vcf.gz) ~{prefix}.GC_anno.vcf | bgzip > ~{prefix}.GC_anno.vcf.gz
    bcftools sort -Oz -o ~{prefix}.GC_anno.sorted.vcf.gz ~{prefix}.GC_anno.vcf.gz
    tabix -p vcf ~{prefix}.GC_anno.sorted.vcf.gz

  >>>


  output {
    File test_vcf = "~{prefix}.GC_anno.vcf"
    File annotated_vcf = "~{prefix}.GC_anno.sorted.vcf.gz"
    File annotated_vcf_idx = "~{prefix}.GC_anno.sorted.vcf.gz.tbi"
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

task BenchmarkSNVs{
  input{
    File comp_vcf
    File base_vcf
    String prefix
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(comp_vcf,"GiB") + size(base_vcf, "GiB"))*2,
    disk_gb: 15 + ceil(size(comp_vcf,"GiB") + size(base_vcf, "GiB"))*2,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String comp_pbaseix = basename(comp_vcf, ".vcf.gz")
  String base_pbaseix = basename(base_vcf, ".vcf.gz")

  command <<<
    set -euxo pipefail

    # use R script to add GC to the vcf
    Rscript -e '


    read_or_empty <- function(file1) {
      if (!file.exists(file1)) {
        stop("File does not exist.")
      }
      
      # Check if file has any lines
      lines <- readLines(file1, warn = FALSE)
      line_count <- length(lines[!grepl("^#", lines)])
      
      if (line_count > 0) {
        df <- read.table(file1, header = FALSE, stringsAsFactors = FALSE)
      } else {
        df <- data.frame(
          V1 = character(),
          V2 = character(),
          V3 = character(),
          V4 = character(),
          V5 = character(),
          stringsAsFactors = FALSE
        )
      }
      
      return(df)
    }

    comp = read_or_empty("~{comp_vcf}")
    base = read_or_empty("~{base_vcf}")

    colnames(comp)[3] = "SVID_comp"
    colnames(base)[3] = "SVID_truth"

    dat=merge(comp[,c(1:5)], base[,c(1:5)], by=c("V1","V2","V4","V5"))

    fp_comp = comp[!comp$SVID_comp%in%dat$SVID_comp, ]
    tp_comp = comp[comp$SVID_comp%in%dat$SVID_comp, ]
    fn_base = base[!base$SVID_truth%in%dat$SVID_truth, ]
    tp_base = base[base$SVID_truth%in%dat$SVID_truth, ]

    write.table(fp_comp, "fp_comp.vcf" , quote=F, sep="\t", col.names=F, row.names=F)
    write.table(tp_comp, "tp_comp.vcf" , quote=F, sep="\t", col.names=F, row.names=F)
    write.table(fn_base, "fn_base.vcf" , quote=F, sep="\t", col.names=F, row.names=F)
    write.table(tp_base, "tp_base.vcf" , quote=F, sep="\t", col.names=F, row.names=F)
    '

    #tabix input vcf
    tabix -p vcf ~{comp_vcf}
    tabix -p vcf ~{base_vcf}

    #extract vcf headers
    bcftools view -h ~{comp_vcf} > comp.header
    bcftools view -h ~{base_vcf} > base.header

    #generate output vcf
    cat comp.header fp_comp.vcf | bgzip > "~{prefix}.fp_comp.vcf.gz"
    cat comp.header tp_comp.vcf | bgzip > "~{prefix}.tp_comp.vcf.gz"
    cat base.header fn_base.vcf | bgzip > "~{prefix}.fn_base.vcf.gz"
    cat base.header tp_base.vcf | bgzip > "~{prefix}.tp_base.vcf.gz"

    >>>


  output {
    File fp_vcf = "~{prefix}.fp_comp.vcf.gz"
    File fn_vcf = "~{prefix}.fn_base.vcf.gz"
    File tp_comp_vcf = "~{prefix}.tp_comp.vcf.gz"
    File tp_base_vcf = "~{prefix}.tp_base.vcf.gz"
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

task BenchmarkSVs{
  input{
    File comp_vcf
    File base_vcf
    File benchmark_bash
    File banchmark_helper_R
    String prefix
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(comp_vcf,"GiB") + size(base_vcf, "GiB"))*2,
    disk_gb: 15 + ceil(size(comp_vcf,"GiB") + size(base_vcf, "GiB"))*2,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euxo pipefail

    #index the input vcf:
    tabix -p vcf ~{comp_vcf}
    tabix -p vcf ~{base_vcf}

    #extract vcf headers
    bcftools view -h ~{comp_vcf} > comp.header
    bcftools view -h ~{base_vcf} > base.header

    #generate headers of query and ref:
    echo -e "#chrom\tstart\tend\tname\tSVTYPE\tSVLEN" > query.header
    echo -e "#chrom\tstart\tend\tVID\tsvtype\tlength\tAF\tsamples" > ref.header

    #generate bed files:
    #svtk vcf2bed -i SVTYPE -i SVLEN ~{comp_vcf} comp.bed
    #svtk vcf2bed -i SVTYPE -i SVLEN ~{base_vcf} base.bed
    bcftools view -H ~{comp_vcf} \
    | awk '{ split($3, a, "_"); if (a[5] < 0) a[5] = -a[5]; print a[1], a[2], a[3], $3, a[4], a[5] }' \
    | sed -e 's/ /\t/g' > comp.bed

    bcftools view -H ~{base_vcf} \
    | awk '{ split($3, a, "_"); if (a[5] < 0) a[5] = -a[5]; print a[1], a[2], a[3], $3, a[4], a[5] }' \
    | sed -e 's/ /\t/g' > base.bed


    #generate query files:
    cat query.header <(tail -n+2 comp.bed | awk '{if ($7!="NA") print}' ) | bgzip > comp.query.gz
    cat query.header <(tail -n+2 base.bed | awk '{if ($7!="NA") print}' ) | bgzip > base.query.gz

    #generate ref files:
    cat ref.header <(tail -n+2 comp.bed | awk '{if ($7!="NA") print}' | sed -e "s/$/\t0\tsample/" ) | bgzip > comp.ref.gz
    cat ref.header <(tail -n+2 base.bed | awk '{if ($7!="NA") print}'  | sed -e "s/$/\t0\tsample/" ) | bgzip > base.ref.gz

    #run comparison:
    bash ~{benchmark_bash} -O comp_vs_base.bed -p comp_vs_base comp.query.gz base.ref.gz
    bash ~{benchmark_bash} -O base_vs_comp.bed -p base_vs_comp base.query.gz comp.ref.gz  


    # use R script to add GC to the vcf
    Rscript -e '

    comp_dat=read.table("base_vs_comp.bed", header=T, sep="\t", comment.char="")
    base_dat=read.table("comp_vs_base.bed", header=T, sep="\t", comment.char="")

    fp_comp = comp_dat[comp_dat[,8]=="NO_OVR",]
    tp_comp = comp_dat[comp_dat[,8]!="NO_OVR",]
    fn_base = base_dat[base_dat[,8]=="NO_OVR",]
    tp_base = base_dat[base_dat[,8]!="NO_OVR",]

    comp_vcf = read.table("~{comp_vcf}")
    base_vcf = read.table("~{base_vcf}")

    fp_comp_vcf = comp_vcf[comp_vcf[,3]%in%fp_comp[,4],]
    tp_comp_vcf = comp_vcf[comp_vcf[,3]%in%tp_comp[,4],]
    fn_base_vcf = base_vcf[base_vcf[,3]%in%fn_base[,4],]
    tp_base_vcf = base_vcf[base_vcf[,3]%in%tp_base[,4],]

    write.table(fp_comp_vcf, "fp_comp.vcf", quote=F, sep="\t", col.names=F, row.names=F)
    write.table(tp_comp_vcf, "tp_comp.vcf", quote=F, sep="\t", col.names=F, row.names=F)
    write.table(fn_base_vcf, "fn_base.vcf", quote=F, sep="\t", col.names=F, row.names=F)
    write.table(tp_base_vcf, "tp_base.vcf", quote=F, sep="\t", col.names=F, row.names=F)

    '


    #generate output vcf

    cat comp.header fp_comp.vcf | bgzip > "~{prefix}.fp_comp.vcf.gz"
    cat comp.header tp_comp.vcf | bgzip > "~{prefix}.tp_comp.vcf.gz"
    cat base.header fn_base.vcf | bgzip > "~{prefix}.fn_base.vcf.gz"
    cat base.header tp_base.vcf | bgzip > "~{prefix}.tp_base.vcf.gz"

    >>>


  output {
    File fp_vcf = "~{prefix}.fp_comp.vcf.gz"
    File fn_vcf = "~{prefix}.fn_base.vcf.gz"
    File tp_comp_vcf = "~{prefix}.tp_comp.vcf.gz"
    File tp_base_vcf = "~{prefix}.tp_base.vcf.gz"
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

task BenchmarkSVsV2{
  #this function would benchmark SVs using the benchmark_bash script, the SV into will be parsed from SVID column

  input{
    File comp_vcf
    File base_vcf
    File benchmark_bash
    File banchmark_helper_R
    String prefix
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(comp_vcf,"GiB") + size(base_vcf, "GiB"))*2,
    disk_gb: 15 + ceil(size(comp_vcf,"GiB") + size(base_vcf, "GiB"))*2,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euxo pipefail

    #index the input vcf:
    tabix -p vcf ~{comp_vcf}
    tabix -p vcf ~{base_vcf}

    #extract vcf headers
    bcftools view -h ~{comp_vcf} > comp.header
    bcftools view -h ~{base_vcf} > base.header

    #generate headers of query and ref:
    echo -e "#chrom\tstart\tend\tname\tSVTYPE\tSVLEN" > query.header
    echo -e "#chrom\tstart\tend\tVID\tsvtype\tlength\tAF\tsamples" > ref.header

    #generate query files:
    paste <(bcftools view -H ~{comp_vcf} | cut -f1,2) <(bcftools view -H ~{comp_vcf} | cut -f3 | sed -e 's/_/\t/g' | cut -f3) <(bcftools view -H ~{comp_vcf} | cut -f3) <(bcftools view -H ~{comp_vcf} | cut -f3 | sed -e 's/_/\t/g' | cut -f4,5) | cat query.header - | bgzip > comp.query.gz

    paste <(bcftools view -H ~{base_vcf} | cut -f1,2) <(bcftools view -H ~{base_vcf} | cut -f3 | sed -e 's/_/\t/g' | cut -f3) <(bcftools view -H ~{base_vcf} | cut -f3) <(bcftools view -H ~{base_vcf} | cut -f3 | sed -e 's/_/\t/g' | cut -f4,5) | cat query.header - | bgzip > base.query.gz

    #generate ref files:
    paste <(bcftools view -H ~{comp_vcf} | cut -f1,2) <(bcftools view -H ~{comp_vcf} | cut -f3 | sed -e 's/_/\t/g' | cut -f3) <(bcftools view -H ~{comp_vcf} | cut -f3) <(bcftools view -H ~{comp_vcf} | cut -f3 | sed -e 's/_/\t/g' | cut -f4,5) | sed -e "s/$/\t0\tsample/" | cat ref.header - | bgzip > comp.ref.gz

    paste <(bcftools view -H ~{base_vcf} | cut -f1,2) <(bcftools view -H ~{base_vcf} | cut -f3 | sed -e 's/_/\t/g' | cut -f3) <(bcftools view -H ~{base_vcf} | cut -f3) <(bcftools view -H ~{base_vcf} | cut -f3 | sed -e 's/_/\t/g' | cut -f4,5) | sed -e "s/$/\t0\tsample/" | cat ref.header - | bgzip > base.ref.gz


    #run comparison:
    bash ~{benchmark_bash} -O comp_vs_base.bed -p comp_vs_base comp.query.gz base.ref.gz
    bash ~{benchmark_bash} -O base_vs_comp.bed -p base_vs_comp base.query.gz comp.ref.gz  


    # use R script to add GC to the vcf
    Rscript -e '

    comp_dat=read.table("base_vs_comp.bed", header=T, sep="\t", comment.char="")
    base_dat=read.table("comp_vs_base.bed", header=T, sep="\t", comment.char="")

    fp_comp = comp_dat[comp_dat[,8]=="NO_OVR",]
    tp_comp = comp_dat[comp_dat[,8]!="NO_OVR",]
    fn_base = base_dat[base_dat[,8]=="NO_OVR",]
    tp_base = base_dat[base_dat[,8]!="NO_OVR",]

    comp_vcf = read.table("~{comp_vcf}")
    base_vcf = read.table("~{base_vcf}")

    fp_comp_vcf = comp_vcf[comp_vcf[,3]%in%fp_comp[,4],]
    tp_comp_vcf = comp_vcf[comp_vcf[,3]%in%tp_comp[,4],]
    fn_base_vcf = base_vcf[base_vcf[,3]%in%fn_base[,4],]
    tp_base_vcf = base_vcf[base_vcf[,3]%in%tp_base[,4],]

    write.table(fp_comp_vcf, "fp_comp.vcf", quote=F, sep="\t", col.names=F, row.names=F)
    write.table(tp_comp_vcf, "tp_comp.vcf", quote=F, sep="\t", col.names=F, row.names=F)
    write.table(fn_base_vcf, "fn_base.vcf", quote=F, sep="\t", col.names=F, row.names=F)
    write.table(tp_base_vcf, "tp_base.vcf", quote=F, sep="\t", col.names=F, row.names=F)

    '


    #generate output vcf

    cat comp.header fp_comp.vcf | bgzip > "~{prefix}.fp_comp.vcf.gz"
    cat comp.header tp_comp.vcf | bgzip > "~{prefix}.tp_comp.vcf.gz"
    cat base.header fn_base.vcf | bgzip > "~{prefix}.fn_base.vcf.gz"
    cat base.header tp_base.vcf | bgzip > "~{prefix}.tp_base.vcf.gz"

    >>>


  output {
    File fp_vcf = "~{prefix}.fp_comp.vcf.gz"
    File fn_vcf = "~{prefix}.fn_base.vcf.gz"
    File tp_comp_vcf = "~{prefix}.tp_comp.vcf.gz"
    File tp_base_vcf = "~{prefix}.tp_base.vcf.gz"
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

task CalculateInheritanceTable {
  input {
    File input_vcf
    File input_vcf_idx
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  Int disk_size = 30 + ceil(size(input_vcf,"GB") * 3)
  Int mem_size =  20 + ceil(size(input_vcf,"GB") * 2)

  String prefix = basename(input_vcf, ".vcf.gz")

  command <<<
    set -euxo pipefail

    bcftools view -H ~{input_vcf} | cut -f10- > gt_table.tsv

    Rscript -e '

    dat=read.table("gt_table.tsv")
    for(i in 1:ncol(dat)){
      dat[,i] = sapply(dat[,i], function(x){gsub("/","|", strsplit(as.character(x),":")[[1]][1])})
    }

    out = data.frame(table(dat))
    out = out[out[,ncol(out)]>0, ]
    write.table(out[,c(ncol(out), 1:(ncol(out)-1))],  "~{prefix}.inheri.stat", quote=F, sep="\t", col.names=F, row.names=F)

    '

  >>>

  output {
    File inheri_stat = "~{prefix}.inheri.stat"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: mem_size,
    disk_gb: disk_size,
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

task ConcatVcfs {
  input {
    Array[File] vcfs
    Array[File]? vcfs_idx
    Boolean merge_sort = false
    Boolean remove_dup = true
    String? outfile_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String merge_flag = if merge_sort then "--allow-overlaps" else ""

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(vcfs, "GB")
  Float compression_factor = 10.0
  Float base_disk_gb = 20.0
  Float base_mem_gb = 10.0

  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + compression_factor)),
    cpu_cores: 1,
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
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    VCFS="~{write_lines(vcfs)}"
    if ~{!defined(vcfs_idx)}; then
      cat ${VCFS} | xargs -n1 tabix
    fi

    # Merge the input VCFs
    bcftools concat -a ~{merge_flag} --output-type z --file-list ${VCFS} --output merged.tmp.vcf.gz

    # Index the temp merged VCF
    tabix -p vcf merged.tmp.vcf.gz


    if [[ ~{remove_dup} == "true" ]]; then
      # Remove duplicates
      bcftools norm -d exact -Oz -o ~{outfile_prefix}.vcf.gz merged.tmp.vcf.gz
    else
      mv merged.tmp.vcf.gz ~{outfile_prefix}.vcf.gz
    fi


    # Index the final deduplicated VCF
    tabix -p vcf ~{outfile_prefix}.vcf.gz

 >>>

  output {
    File concat_vcf = "~{outfile_prefix}.vcf.gz"
    File concat_vcf_idx =  "~{outfile_prefix}.vcf.gz.tbi"
  }
}

task MergeVcfs {
  input {
    Array[File] vcfs
    Array[File]? vcfs_idx
    Boolean merge_sort = false
    String output_name
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override

  }


  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(vcfs, "GB")
  Float compression_factor = 10.0
  Float base_disk_gb = 20.0
  Float base_mem_gb = 10.0

  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + compression_factor)),
    cpu_cores: 1,
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
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    VCFS="~{write_lines(vcfs)}"
    if ~{!defined(vcfs_idx)}; then
      cat ${VCFS} | xargs -n1 tabix
    fi

    bcftools merge --merge none ~{sep=' ' vcfs} -Oz -o ~{output_name}
    tabix -p vcf ~{output_name}

 >>>

  output {
    File output_merged_vcf = "~{output_name}"
    File output_merged_vcf_idx =  "~{output_name}.tbi"
  }
}

task ConvertBubblesToBiallelic{
  input {
    File input_vcf
    File input_vcf_idx

    File panel_biallelic_vcf
    File panel_biallelic_vcf_idx
    File convert_to_biallelic_script

    Boolean sort = false
    String docker_image
    RuntimeAttr? runtime_attr_override
  }


  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(input_vcf, "GB")
  Float compression_factor = 10.0
  Float base_disk_gb = 20.0
  Float base_mem_gb = 10.0

  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + compression_factor)),
    cpu_cores: 1,
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
    docker: docker_image
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  String prefix = basename(input_vcf, ".vcf.gz")

  command <<<
    set -euo pipefail

    zcat ~{input_vcf} | python3 ~{convert_to_biallelic_script} "~{panel_biallelic_vcf}" > "temp.vcf"
    bgzip  "temp.vcf"

    if [[ ~{sort} == "true" ]]; then
      # sort vcf
      bcftools sort "temp.vcf.gz" -O z -o "~{prefix}_biallelic.vcf.gz"
    else
      mv temp.vcf.gz "~{prefix}_biallelic.vcf.gz"
    fi


    tabix -p vcf  "~{prefix}_biallelic.vcf.gz"

 >>>

  output {
    File biallelic_vcf = "~{prefix}_biallelic.vcf.gz"
    File biallelic_vcf_idx =  "~{prefix}_biallelic.vcf.gz.tbi"
  }
}

task EvaluateInheriByGQ {
  input {
    File vcf_file         # Input VCF file (bgzipped or plain)
    File vcf_idx_file
    File inheri_stat  # 2-column TSV: SVID <tab> annotation
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(vcf_file,'.vcf.gz')

  command <<<
    set -euxo pipefail

    # use R script to add GC to the vcf
    Rscript -e '

    init_dataframe <- function(colnames, nrows = 0) {
      if (!is.character(colnames)) {
        stop("Column names must be a character vector.")
      }
      if (!is.numeric(nrows) || nrows < 0 || nrows %% 1 != 0) {
        stop("Number of rows must be a non-negative integer.")
      }

      # Create a list of empty character vectors of length nrows
      data_list <- setNames(
        replicate(length(colnames), rep(NA_character_, nrows), simplify = FALSE),
        colnames
      )
      
      # Create data frame
      df <- as.data.frame(data_list, stringsAsFactors = FALSE)
      return(df)
    }

    read_in_vcf_gt_gq<-function(vcf_file){
      # Read in all lines, excluding meta-information lines (those starting with "##")
      vcf_lines <- readLines(vcf_file)
      vcf_lines <- vcf_lines[!startsWith(vcf_lines, "##")]

      # Split the last header line to get column names (starts with "#CHROM")
      header <- strsplit(vcf_lines[1], "\t")[[1]]
      sample_names <- header[10:length(header)]

      # Read the VCF body (data part)
      vcf_data <- read.table(text = vcf_lines[-1], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      colnames(vcf_data) <- header

      #generate the table of GT and GQ of each sample
      sample = sample_names[1]
      print(sample)
      result = data.frame(sapply(vcf_data[,colnames(vcf_data) == sample], function(x){strsplit(as.character(x),":")[[1]][1]}))
      colnames(result)[ncol(result)] = paste(sample,"GT", sep="_")
      result[,2]  = data.frame(sapply(vcf_data[,colnames(vcf_data) == sample], function(x){strsplit(as.character(x),":")[[1]][2]}))
      colnames(result)[ncol(result)] = paste(sample,"GQ", sep="_")
      for(sample in sample_names[2:length(sample_names)]){
        print(sample)
        result[,ncol(result)+1]  = data.frame(sapply(vcf_data[,colnames(vcf_data) == sample], function(x){strsplit(as.character(x),":")[[1]][1]}))
        colnames(result)[ncol(result)] = paste(sample,"GT", sep="_")
        result[,ncol(result)+1]  = data.frame(sapply(vcf_data[,colnames(vcf_data) == sample], function(x){strsplit(as.character(x),":")[[1]][2]}))
        colnames(result)[ncol(result)] = paste(sample,"GQ", sep="_")
      }

      # Optional: Add CHROM and POS if you want to track variant position
      result$CHROM <- vcf_data$`#CHROM`
      result$POS <- vcf_data$POS

      # Reorder columns to put CHROM and POS first
      result <- result[, c("CHROM", "POS", setdiff(names(result), c("CHROM", "POS")))]
      return(result)
    }

    generate_gt_stat<-function(dat,inheri_table){
      #assuming 6 columns in dat: fa_GT, fa_GQ, mo_GT, mo_GQ, pb_GT, pb_GQ
      dat=dat[dat[,2]!="." & dat[,4]!="." & dat[,6]!=".",]
      dat = dat[dat[,1]!="0/0" | dat[,3]!="0/0" | dat[,5]!="0/0", ]
      dat[,2] = as.integer(dat[,2])
      dat[,4] = as.integer(dat[,4])
      dat[,6] = as.integer(dat[,6])
      dat=dat[order(dat[,6]),]

      #add decile to the GQ of proband
      dat[,7] = as.integer((c(1:nrow(dat))-1) /nrow(dat) *10)

      stat = data.frame(table(dat[,c(1,3,5,7)]))
      stat[,1] = gsub("/","|", stat[,1])
      stat[,2] = gsub("/","|", stat[,2])
      stat[,3] = gsub("/","|", stat[,3])
      colnames(stat)=c("fa","mo","pb","gt_cate","SV_count")
      stat = merge(stat, inheri_table, by=c("fa","mo","pb"))
      stat$SV_count = as.integer(stat$SV_count)
      return(stat)
    }

    categorize_gt_stat<-function(stat, dat){
      #assuming 6 columns in dat: fa_GT, fa_GQ, mo_GT, mo_GQ, pb_GT, pb_GQ
      dat=dat[dat[,2]!="." & dat[,4]!="." & dat[,6]!=".",]
      dat = dat[dat[,1]!="0/0" | dat[,3]!="0/0" | dat[,5]!="0/0", ]
      dat[,2] = as.integer(dat[,2])
      dat[,4] = as.integer(dat[,4])
      dat[,6] = as.integer(dat[,6])
      dat=dat[order(dat[,6]),]

      #add decile to the GQ of proband
      dat[,7] = as.integer((c(1:nrow(dat))-1) /nrow(dat) *10)
      
      #define the categories that infer de novo or mendelian errors:
      dnv = c("de_novo")
      me = c("Mendelian_Error_Fa_Pb" , "Mendelian_Error_Mo_Pb", "Mendelian_Error_Fa_Mo_Pb", "Mendelian_Error")

      error_rate_table = data.frame(table(stat$gt_cate))

      #calculate dnv counts by GQ decile:
      error_rate_table[,2] = sapply(error_rate_table[,1], function(x){sum(stat[stat$gt_cate == x & stat$category%in%dnv,]$SV_count)})
      #calculate mendelian error counts by GQ decile:
      error_rate_table[,3] = sapply(error_rate_table[,1], function(x){sum(stat[stat$gt_cate == x & stat$category%in%me & stat$pb!="0|0",]$SV_count)})
      #calculate SV counts by GQ decile:
      error_rate_table[,4] = sapply(error_rate_table[,1], function(x){sum(stat[stat$gt_cate == x & stat$pb!="0|0",]$SV_count)})

      #calculate mendelian error counts by GQ decile:
      error_rate_table[,5] = sapply(error_rate_table[,1], function(x){sum(stat[stat$gt_cate == x & stat$category%in%me,]$SV_count)})
      #calculate SV counts by GQ decile:
      error_rate_table[,6] = sapply(error_rate_table[,1], function(x){sum(stat[stat$gt_cate == x & stat$category!="REF",]$SV_count)})
      
      colnames(error_rate_table) = c("gt_cate","dnv","Mendelian_Error_child","in_child","Mendelian_Error_fam","in_fam")
      
      #calculate accumulative dnv rate
      #calculate accumulative dnv rate
      error_rate_table[,7] = sapply(c(1:nrow(error_rate_table)), function(x){sum(error_rate_table[c(x:nrow(error_rate_table)),]$dnv) / sum(error_rate_table[c(x:nrow(error_rate_table)),]$in_child)})
      #calculate accumulative mendelian error rate in child
      error_rate_table[,8] = sapply(c(1:nrow(error_rate_table)), function(x){sum(error_rate_table[c(x:nrow(error_rate_table)),]$Mendelian_Error_child) / sum(error_rate_table[c(x:nrow(error_rate_table)),]$in_child)})
      #calculate accumulative SV count in child
      error_rate_table[,9] = sapply(c(1:nrow(error_rate_table)), function(x){sum(error_rate_table[c(x:nrow(error_rate_table)),]$in_child)})
      #calculate accumulative mendelian error rate in family
      error_rate_table[,10] = sapply(c(1:nrow(error_rate_table)), function(x){sum(error_rate_table[c(x:nrow(error_rate_table)),]$Mendelian_Error_fam) / sum(error_rate_table[c(x:nrow(error_rate_table)),]$in_fam)})
      #calculate accumulative SV count in family
      error_rate_table[,11] = sapply(c(1:nrow(error_rate_table)), function(x){sum(error_rate_table[c(x:nrow(error_rate_table)),]$in_fam)})
      
      #add the GQ range for each decile
      error_rate_table[,12] = sapply(error_rate_table$gt_cate, function(x){min(dat[dat[,7]==x,][,6])})
      error_rate_table[,13] = sapply(error_rate_table$gt_cate, function(x){max(dat[dat[,7]==x,][,6])})
      colnames(error_rate_table)[c(5:ncol(error_rate_table))] = c("dnv_rate","Mendelian_Error_rate", "total_sv_count", "min_GQ","max_GQ")
      return(error_rate_table)
    }

    # Function to print help message
    print_help <- function() {
      cat("
      Usage: Rscript run_gt_categorization.R -i <input.vcf> -t <inheri_table.tsv> -o <output.tsv>

      Options:
        -i    Path to input VCF file
        -t    Path to inheritance table (tab-delimited, with header)
        -o    Path to output table (tab-delimited)
        -h, --help    Show this help message and exit

      Description:
        This script reads a VCF file, processes genotypes (GT) and genotype quality (GQ),
        integrates with an inheritance table, generates GT statistics,
        categorizes them, and writes the results to an output file.
      ")
    }

    # Assign variables
    inheri_table <- read.table("~{inheri_stat}", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    gt_table <- read_in_vcf_gt_gq("~{vcf_file}")
    gt_stat <- generate_gt_stat(gt_table[, 3:ncol(gt_table)], inheri_table)
    gt_cate <- categorize_gt_stat(gt_stat, gt_table[, 3:ncol(gt_table)])

    # Write output
    write.table(gt_cate, "~{prefix}.inheri_by_GQ.stat", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

    '

  >>>

  output {
    File inheri_by_GQ_stat = "~{prefix}.inheri_by_GQ.stat"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(vcf_file,"GiB")*3),
    disk_gb: 15 + ceil(size(vcf_file,"GiB")*3),
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

task ExtractVariantSites {
  input {
    File input_vcf
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(input_vcf, '.vcf.gz')
  command <<<
    set -euxo pipefail

    python3 <<CODE
    import pysam

    def determine_svtype(ref, alt, rec):
        if "<" in alt and ">" in alt:
            if "SVTYPE" in rec.info.keys():
              return rec.info["SVTYPE"]
            else:
              return gsub(">", "", gsub("<", "", alt))
        else:
          if len(ref) == 1 and len(alt) == 1:
              return "SNV"
          elif len(ref) < len(alt):
              return "INS"
          elif len(ref) > len(alt):
              return "DEL"
          else:
              return "OTH"

    def determine_svlen(ref, alt, rec):
        if "<" in alt and ">" in alt:
          if "SVLEN" in rec.info.keys():
            return str(rec.info["SVLEN"])
          else:
            return 50
        else:
          return str(abs(len(alt) - len(ref)))


    vcf_in = pysam.VariantFile("~{input_vcf}")
    #rewrite the vcfs in the a new file with the SVID updated to be consistent with the bed file
    vcf_out= pysam.VariantFile("~{prefix}.SVID_updated.vcf.gz", 'w', header = vcf_in.header)
    
    with open("~{prefix}.variant_sites.bed", "w") as out:
        for rec in vcf_in.fetch():
            chrom = rec.contig
            pos = str(rec.pos)
            ref = rec.ref
            end = str(rec.pos + len(ref) -1)  # .stop is preferred over parsing INFO['END']
            alt = rec.alts[0] if rec.alts else "."
            if alt!="<INV>":
              svtype = determine_svtype(ref, alt, rec)
              svlen  = determine_svlen(ref, alt, rec)
              ID = f"{chrom}_{pos}_{end}_{svtype}_{svlen}"
              rec.id = ID
              out.write('\t'.join([chrom, pos, end, ID, svtype, svlen]) + "\n")
              vcf_out.write(rec)
    CODE

    tabix -p vcf "~{prefix}.SVID_updated.vcf.gz"

  >>>

  output {
    File variant_sites = "~{prefix}.variant_sites.bed"
    File updated_vcf = "~{prefix}.SVID_updated.vcf.gz"
    File updated_vcf_idx = "~{prefix}.SVID_updated.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 20,
    disk_gb: 20 + ceil(size(input_vcf,"GiB")*2),
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

task ExtractTrioVCF {
  input {
    File input_vcf
    File sample_file
    String family_id
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euxo pipefail

    bcftools view -S ~{sample_file} -Oz -o ~{family_id}.trio.vcf.gz ~{input_vcf}
    tabix -p vcf "~{family_id}.trio.vcf.gz"

  >>>

  output {
    File output_vcf = "~{family_id}.trio.vcf.gz"
    File output_vcf_idx = "~{family_id}.trio.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: 20,
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

task ExtractChromosomeVcf {
  input {
    File input_vcf
    File input_vcf_idx
    String chromosome
    String docker_image
    File? monitoring_script
    RuntimeAttr? runtime_attr_override
  }

  
  command <<<
    set -euxo pipefail

    # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
    touch monitoring.log
    if [ -s ~{monitoring_script} ]; then
        bash ~{monitoring_script} > monitoring.log &
    fi


    bcftools view -r ~{chromosome} ~{input_vcf} -Oz -o ~{chromosome}.vcf.gz
    tabix -p vcf "~{chromosome}.vcf.gz"

  >>>

  output {
    File output_vcf = "~{chromosome}.vcf.gz"
    File output_vcf_idx = "~{chromosome}.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: ceil(size(input_vcf,"GiB")*2)+10,
    disk_gb: ceil(size(input_vcf,"GiB")*2)+20,
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

task ExtractVariantIndividualGenome {
  input {
    File vcf_file
    String sample_id
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euxo pipefail

    # Extract variants for the sample
    bcftools view -s ~{sample_id} ~{vcf_file} -Oz -o ~{sample_id}.vcf.gz
    tabix -p vcf ~{sample_id}.vcf.gz

    # Filter informative, alternative genotypes
    bcftools view -c 1 ~{sample_id}.vcf.gz -Oz -o ~{sample_id}.non_ref.vcf.gz
    tabix -p vcf ~{sample_id}.non_ref.vcf.gz

  >>>

  output {
    File non_ref_vcf = "~{sample_id}.non_ref.vcf.gz"
    File non_ref_vcf_idx = "~{sample_id}.non_ref.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: 20,
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

task ExtractChromosomeVariants {
  input {
    File? input_vcf       # .vcf.gz
    File? input_vcf_index # .vcf.gz.tbi
    String chromosome    # e.g. "chr1" or "1"
    String output_name   # e.g. "chr1.vcf.gz"
    String docker_image = "biocontainers/bcftools:v1.17-1-deb-py3"
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euxo pipefail

    # Extract chromosome-specific variants
    bcftools view -r ~{chromosome} ~{input_vcf} -Oz -o ~{output_name}

    # Index the output VCF
    tabix -p vcf ~{output_name}
  >>>

  output {
    File chr_vcf = output_name
    File chr_vcf_idx = output_name + ".tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb:  10,
    disk_gb: 10 + ceil(size(input_vcf,"GiB")*1.5),
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

task FilterVcfByAnotherVcf {
  #this task will filter vcf_file to the variants in filter_vcf_b by matching chrom, pos, svid, ref, and alts
  input {
    File vcf_file
    File vcf_idx
    File? vcf_file_b
    File? vcf_file_b_idx
    String docker_image
    RuntimeAttr? runtime_attr_override
  }


  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb:  10 + ceil(size(vcf_file,"GiB") + size(vcf_file_b, "GiB"))*4,
    disk_gb: 20 + ceil(size(vcf_file,"GiB") + size(vcf_file_b, "GiB"))*4,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String prefix = basename(vcf_file, "vcf.gz")

  command <<<
    set -euo pipefail

    # Write Python script
    python3 <<CODE
    import sys
    import pysam

    vcf_a_path = "~{vcf_file}"
    vcf_b_path = "~{vcf_file_b}"
    output_path = "~{prefix}.filtered.vcf.gz"

    def variant_key(record):
        return (record.contig, record.pos, record.id, record.ref, tuple(record.alts))

    # Load B into set
    vcf_b = pysam.VariantFile(vcf_b_path)
    b_variants = set(variant_key(rec) for rec in vcf_b)
    vcf_b.close()

    # Filter A
    vcf_a = pysam.VariantFile(vcf_a_path)
    vcf_out = pysam.VariantFile(output_path, "w", header=vcf_a.header)

    for rec in vcf_a:
        if variant_key(rec) in b_variants:
            vcf_out.write(rec)

    vcf_a.close()
    vcf_out.close()
    
    CODE

    tabix -p vcf "~{prefix}.filtered.vcf.gz"
  >>>

  output {
    File filtered_vcf = "~{prefix}.filtered.vcf.gz"
    File filtered_vcf_idx = "~{prefix}.filtered.vcf.gz.tbi"
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

task CalcuCompStat{
  input {
    File tp_query 
    File tp_ref
    File fp_query 
    File fp_ref

    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(tp_query, ".tp_query.vcf.gz")

  command <<<
    set -euxo pipefail

    # Index VCF if needed
    if [[ ! -f "~{tp_query}.tbi" && ! -f "~{tp_query}.csi" ]]; then
      bcftools index ~{tp_query}
    fi

    if [[ ! -f "~{fp_query}.tbi" && ! -f "~{fp_query}.csi" ]]; then
      bcftools index ~{fp_query}
    fi

    if [[ ! -f "~{tp_ref}.tbi" && ! -f "~{tp_ref}.csi" ]]; then
      bcftools index ~{tp_ref}
    fi

    if [[ ! -f "~{fp_ref}.tbi" && ! -f "~{fp_ref}.csi" ]]; then
      bcftools index ~{fp_ref}
    fi


    echo -e "Genomic_Context\tSVTYPE\tSizeRange\tfp_query\ttp_query\tfp_ref\ttp_ref" > "~{prefix}.stat"

    bcftools  view -i 'INFO/GC="US"'  ~{fp_query} | cut -f1-3 > fp_query.US
    bcftools  view -i 'INFO/GC="US"'  ~{tp_query} | cut -f1-3 > tp_query.US
    bcftools  view -i 'INFO/GC="US"'  ~{fp_ref} | cut -f1-3 > fp_ref.US
    bcftools  view -i 'INFO/GC="US"'  ~{tp_ref} | cut -f1-3 > tp_ref.US

    paste <(grep -v "#" fp_query.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="SNV") print}'  | wc -l) <(grep -v "#" tp_query.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="SNV") print}'  | wc -l) <(grep -v "#" fp_ref.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="SNV") print}'  | wc -l) <(grep -v "#" tp_ref.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="SNV") print}'  | wc -l) | sed -e 's/^/US\tSNV\tSNV\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5<30) print}'  | wc -l) <(grep -v "#" tp_query.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5<30) print}'  | wc -l) <(grep -v "#" fp_ref.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5<30) print}'  | wc -l) <(grep -v "#" tp_ref.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5<30) print}'  | wc -l) | sed -e 's/^/US\tDEL\tIndel_sm\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" tp_query.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" fp_ref.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" tp_ref.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>29 && $5<50) print}'  | wc -l) | sed -e 's/^/US\tDEL\tIndel_lg\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>49 ) print}'  | wc -l) <(grep -v "#" tp_query.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>49 ) print}'  | wc -l) <(grep -v "#" fp_ref.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>49 ) print}'  | wc -l) <(grep -v "#" tp_ref.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>49 ) print}'  | wc -l) | sed -e 's/^/US\tDEL\tSV\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5<30) print}'  | wc -l) <(grep -v "#" tp_query.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5<30) print}'  | wc -l) <(grep -v "#" fp_ref.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5<30) print}'  | wc -l) <(grep -v "#" tp_ref.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5<30) print}'  | wc -l) | sed -e 's/^/US\tINS\tIndel_sm\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" tp_query.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" fp_ref.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" tp_ref.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>29 && $5<50) print}'  | wc -l) | sed -e 's/^/US\tINS\tIndel_lg\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>49 ) print}'  | wc -l) <(grep -v "#" tp_query.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>49 ) print}'  | wc -l) <(grep -v "#" fp_ref.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>49 ) print}'  | wc -l) <(grep -v "#" tp_ref.US | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>49 ) print}'  | wc -l) | sed -e 's/^/US\tINS\tSV\t/'  >> "~{prefix}.stat"


    bcftools  view -i 'INFO/GC="RM"'  ~{fp_query} | cut -f1-3 > fp_query.RM
    bcftools  view -i 'INFO/GC="RM"'  ~{tp_query} | cut -f1-3 > tp_query.RM
    bcftools  view -i 'INFO/GC="RM"'  ~{fp_ref} | cut -f1-3 > fp_ref.RM
    bcftools  view -i 'INFO/GC="RM"'  ~{tp_ref} | cut -f1-3 > tp_ref.RM

    paste <(grep -v "#" fp_query.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="SNV") print}'  | wc -l) <(grep -v "#" tp_query.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="SNV") print}'  | wc -l) <(grep -v "#" fp_ref.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="SNV") print}'  | wc -l) <(grep -v "#" tp_ref.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="SNV") print}'  | wc -l) | sed -e 's/^/RM\tSNV\tSNV\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5<30) print}'  | wc -l) <(grep -v "#" tp_query.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5<30) print}'  | wc -l) <(grep -v "#" fp_ref.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5<30) print}'  | wc -l) <(grep -v "#" tp_ref.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5<30) print}'  | wc -l) | sed -e 's/^/RM\tDEL\tIndel_sm\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" tp_query.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" fp_ref.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" tp_ref.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>29 && $5<50) print}'  | wc -l) | sed -e 's/^/RM\tDEL\tIndel_lg\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>49 ) print}'  | wc -l) <(grep -v "#" tp_query.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>49 ) print}'  | wc -l) <(grep -v "#" fp_ref.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>49 ) print}'  | wc -l) <(grep -v "#" tp_ref.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>49 ) print}'  | wc -l) | sed -e 's/^/RM\tDEL\tSV\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5<30) print}'  | wc -l) <(grep -v "#" tp_query.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5<30) print}'  | wc -l) <(grep -v "#" fp_ref.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5<30) print}'  | wc -l) <(grep -v "#" tp_ref.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5<30) print}'  | wc -l) | sed -e 's/^/RM\tINS\tIndel_sm\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" tp_query.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" fp_ref.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" tp_ref.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>29 && $5<50) print}'  | wc -l) | sed -e 's/^/RM\tINS\tIndel_lg\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>49 ) print}'  | wc -l) <(grep -v "#" tp_query.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>49 ) print}'  | wc -l) <(grep -v "#" fp_ref.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>49 ) print}'  | wc -l) <(grep -v "#" tp_ref.RM | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>49 ) print}'  | wc -l) | sed -e 's/^/RM\tINS\tSV\t/'  >> "~{prefix}.stat"


    bcftools  view -i 'INFO/GC="SD"'  ~{fp_query} | cut -f1-3 > fp_query.SD
    bcftools  view -i 'INFO/GC="SD"'  ~{tp_query} | cut -f1-3 > tp_query.SD
    bcftools  view -i 'INFO/GC="SD"'  ~{fp_ref} | cut -f1-3 > fp_ref.SD
    bcftools  view -i 'INFO/GC="SD"'  ~{tp_ref} | cut -f1-3 > tp_ref.SD

    paste <(grep -v "#" fp_query.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="SNV") print}'  | wc -l) <(grep -v "#" tp_query.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="SNV") print}'  | wc -l) <(grep -v "#" fp_ref.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="SNV") print}'  | wc -l) <(grep -v "#" tp_ref.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="SNV") print}'  | wc -l) | sed -e 's/^/SD\tSNV\tSNV\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5<30) print}'  | wc -l) <(grep -v "#" tp_query.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5<30) print}'  | wc -l) <(grep -v "#" fp_ref.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5<30) print}'  | wc -l) <(grep -v "#" tp_ref.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5<30) print}'  | wc -l) | sed -e 's/^/SD\tDEL\tIndel_sm\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" tp_query.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" fp_ref.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" tp_ref.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>29 && $5<50) print}'  | wc -l) | sed -e 's/^/SD\tDEL\tIndel_lg\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>49 ) print}'  | wc -l) <(grep -v "#" tp_query.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>49 ) print}'  | wc -l) <(grep -v "#" fp_ref.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>49 ) print}'  | wc -l) <(grep -v "#" tp_ref.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>49 ) print}'  | wc -l) | sed -e 's/^/SD\tDEL\tSV\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5<30) print}'  | wc -l) <(grep -v "#" tp_query.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5<30) print}'  | wc -l) <(grep -v "#" fp_ref.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5<30) print}'  | wc -l) <(grep -v "#" tp_ref.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5<30) print}'  | wc -l) | sed -e 's/^/SD\tINS\tIndel_sm\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" tp_query.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" fp_ref.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" tp_ref.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>29 && $5<50) print}'  | wc -l) | sed -e 's/^/SD\tINS\tIndel_lg\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>49 ) print}'  | wc -l) <(grep -v "#" tp_query.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>49 ) print}'  | wc -l) <(grep -v "#" fp_ref.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>49 ) print}'  | wc -l) <(grep -v "#" tp_ref.SD | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>49 ) print}'  | wc -l) | sed -e 's/^/SD\tINS\tSV\t/'  >> "~{prefix}.stat"


    bcftools  view -i 'INFO/GC="SR"'  ~{fp_query} | cut -f1-3 > fp_query.SR
    bcftools  view -i 'INFO/GC="SR"'  ~{tp_query} | cut -f1-3 > tp_query.SR
    bcftools  view -i 'INFO/GC="SR"'  ~{fp_ref} | cut -f1-3 > fp_ref.SR
    bcftools  view -i 'INFO/GC="SR"'  ~{tp_ref} | cut -f1-3 > tp_ref.SR

    paste <(grep -v "#" fp_query.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="SNV") print}'  | wc -l) <(grep -v "#" tp_query.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="SNV") print}'  | wc -l) <(grep -v "#" fp_ref.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="SNV") print}'  | wc -l) <(grep -v "#" tp_ref.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="SNV") print}'  | wc -l) | sed -e 's/^/SR\tSNV\tSNV\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5<30) print}'  | wc -l) <(grep -v "#" tp_query.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5<30) print}'  | wc -l) <(grep -v "#" fp_ref.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5<30) print}'  | wc -l) <(grep -v "#" tp_ref.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5<30) print}'  | wc -l) | sed -e 's/^/SR\tDEL\tIndel_sm\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" tp_query.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" fp_ref.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" tp_ref.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>29 && $5<50) print}'  | wc -l) | sed -e 's/^/SR\tDEL\tIndel_lg\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>49 ) print}'  | wc -l) <(grep -v "#" tp_query.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>49 ) print}'  | wc -l) <(grep -v "#" fp_ref.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>49 ) print}'  | wc -l) <(grep -v "#" tp_ref.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="DEL" && $5>49 ) print}'  | wc -l) | sed -e 's/^/SR\tDEL\tSV\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5<30) print}'  | wc -l) <(grep -v "#" tp_query.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5<30) print}'  | wc -l) <(grep -v "#" fp_ref.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5<30) print}'  | wc -l) <(grep -v "#" tp_ref.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5<30) print}'  | wc -l) | sed -e 's/^/SR\tINS\tIndel_sm\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" tp_query.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" fp_ref.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>29 && $5<50) print}'  | wc -l) <(grep -v "#" tp_ref.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>29 && $5<50) print}'  | wc -l) | sed -e 's/^/SR\tINS\tIndel_lg\t/'  >> "~{prefix}.stat"
    paste <(grep -v "#" fp_query.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>49 ) print}'  | wc -l) <(grep -v "#" tp_query.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>49 ) print}'  | wc -l) <(grep -v "#" fp_ref.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>49 ) print}'  | wc -l) <(grep -v "#" tp_ref.SR | cut -f3 | sed -e 's/_/\t/g' | awk '{if ($4=="INS" && $5>49 ) print}'  | wc -l) | sed -e 's/^/SR\tINS\tSV\t/'  >> "~{prefix}.stat"

  >>>

  output {
    File comp_stat = "~{prefix}.stat"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: 20 + ceil(size(tp_query, "GiB") + size(fp_query, "GiB") + size(tp_ref, "GiB") + size(fp_ref, "GiB"))*4,
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

task CalcuPlotCompResults{
  input {
    File tp_query 
    File tp_ref
    File fp_query 
    File fp_ref

    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(fp_query, ".fp_query.vcf.gz")

  command <<<
    set -euxo pipefail

    Rscript -e '

    init_empty_sv_table <- function() {
      df <- data.frame(
        Genomic_Context = character(),
        SVTYPE = character(),
        Freq = numeric(),
        Freq = numeric(),
        V4 = character(),
        stringsAsFactors = FALSE
      )
      return(df)
    }

    calcu_counts_by_variant_type_and_GC<-function(dat){
      snv_stat = data.frame(table(dat[dat$SVTYPE=="SNV",c("Genomic_Context", "SVTYPE")]))
      indel_sm_stat = data.frame(table(dat[dat$SVLEN>0 & dat$SVLEN<30,c("Genomic_Context", "SVTYPE")]))
      indel_lg_stat = data.frame(table(dat[dat$SVLEN>29 & dat$SVLEN<50,c("Genomic_Context", "SVTYPE")]))
      sv_stat = data.frame(table(dat[dat$SVLEN>49,c("Genomic_Context", "SVTYPE")]))
      if(nrow(snv_stat)>0){
        snv_stat[,ncol(snv_stat)+1] = "SNV"
      } else{ snv_stat = init_empty_sv_table()  }

      if(nrow(indel_sm_stat)>0){
        indel_sm_stat[,ncol(indel_sm_stat)+1] = "Indel_sm"
      } else{ indel_sm_stat = init_empty_sv_table()  }

      if(nrow(indel_lg_stat)>0){
        indel_lg_stat[,ncol(indel_lg_stat)+1] = "Indel_lg"
      } else{ indel_lg_stat = init_empty_sv_table()  }

      if(nrow(sv_stat)>0){
        sv_stat[,ncol(sv_stat)+1] = "SV"
      } else{ sv_stat = init_empty_sv_table()  }
      
      out = rbind(snv_stat, indel_sm_stat, indel_lg_stat, sv_stat)
      return(out)      
    }

    readin_benchmark_results<-function(filename){
      print("read in data ... ")
      dat=read.table(filename)
      print("extract svtype information ... ")
      dat[,ncol(dat)+1] = sapply(dat[,3], function(x){strsplit(as.character(x),"_")[[1]][4]})
      colnames(dat)[ncol(dat)] = "SVTYPE"
      print("extract svlen information ... ")
      dat[,ncol(dat)+1] = sapply(dat[,3], function(x){strsplit(as.character(x),"_")[[1]][5]})
      colnames(dat)[ncol(dat)] = "SVLEN"
      dat$SVLEN = as.integer(dat$SVLEN)

      if(nrow(dat[grepl("<",dat[,5]),])>0){
        tmp1 = dat[!grepl("<",dat[,5]),]
        tmp = dat[grepl("<",dat[,5]),]
        tmp$SVTYPE = sapply(tmp[,8], function(x){extract_info_col(x,"SVTYPE")})
        tmp$SVLEN = sapply(tmp[,8], function(x){extract_info_col(x,"SVLEN")})
        tmp$SVLEN = abs(as.integer(tmp$SVLEN))
        dat = rbind(tmp, tmp1)
      }


      print("extract genomic context information ... ")
      dat[,ncol(dat)+1] = sapply(dat[,8], function(x){extract_info_col(x,"GC")})
      colnames(dat)[ncol(dat)] = "Genomic_Context"
      return(dat)
    }

    extract_info_col<-function(info, col){
      #col represents columns in info, eg. GC, SVTYPE, SVLEN
      tmp = strsplit(as.character(info),";")[[1]]
      out = tmp[grepl(paste(col, "=", sep=""), tmp)]
      return(strsplit(as.character(out),"=")[[1]][2])
    }

    extract_GC<-function(info){
      tmp = strsplit(as.character(info),";")[[1]]
      out = tmp[grepl("GC=", tmp)]
      return(strsplit(as.character(out),"=")[[1]][2])
    }

    extract_SVTYPE_from_info<-function(info){
      tmp = strsplit(as.character(info),";")[[1]]
      out = tmp[grepl("SVTYPE=", tmp)]
      return(strsplit(as.character(out),"=")[[1]][2])
    }

    extract_SVLEN_from_info<-function(info){
      tmp = strsplit(as.character(info),";")[[1]]
      out = tmp[grepl("SVLEN=", tmp)]
      return(strsplit(as.character(out),"=")[[1]][2])
    }

    add_ovr_bar.US_RM<-function(query_vs_ref.stat,size_range, svtype, y_axis){
      rec_width = .3
      rect(-1,y_axis - rec_width, 1, y_axis+rec_width, col=colorBlindGrey8[1])
      snv_usrm = colSums(query_vs_ref.stat[query_vs_ref.stat$SizeRange==size_range & query_vs_ref.stat$SVTYPE==svtype & query_vs_ref.stat$Genomic_Context%in%c("US","RM"),c(4:7)])
      rect(-snv_usrm[2]/sum(snv_usrm[c(1,2)]),y_axis - rec_width, snv_usrm[4]/sum(snv_usrm[c(3,4)]), y_axis+rec_width, col=colorBlindGrey8[2])
      text(-1.2, y_axis+rec_width, format(sum(snv_usrm[c(1,2)]), big.mark = ","), adj=c(0,0))
      text(1.2, y_axis+rec_width, format(sum(snv_usrm[c(3,4)]), big.mark = ","), adj=c(1,0))
    }

    add_ovr_bar.SD_SR<-function(query_vs_ref.stat,size_range, svtype, y_axis){
      rec_width = .3
      rect(-1,y_axis - rec_width, 1, y_axis+rec_width, col=colorBlindGrey8[1])
      snv_usrm = colSums(query_vs_ref.stat[query_vs_ref.stat$SizeRange==size_range & query_vs_ref.stat$SVTYPE==svtype & query_vs_ref.stat$Genomic_Context%in%c("SD","SR"),c(4:7)])
      rect(-snv_usrm[2]/sum(snv_usrm[c(1,2)]),y_axis - rec_width, snv_usrm[4]/sum(snv_usrm[c(3,4)]), y_axis+rec_width, col=colorBlindGrey8[2])
      text(-1.2, y_axis+rec_width, format(sum(snv_usrm[c(1,2)]), big.mark = ","), adj=c(0,0))
      text(1.2, y_axis+rec_width, format(sum(snv_usrm[c(3,4)]), big.mark = ","), adj=c(1,0))
    }

    colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    plot_ovr_bar<-function(fp_query_file,  tp_query_file, fp_ref_file, tp_ref_file, prefix){
      fp_query=readin_benchmark_results(fp_query_file)
      tp_query=readin_benchmark_results(tp_query_file)
      fp_query = fp_query[!fp_query$V3%in%tp_query$V3,]
      fp_ref=readin_benchmark_results(fp_ref_file)
      tp_ref=readin_benchmark_results(tp_ref_file)
      fp_ref = fp_ref[!fp_ref$V3%in%tp_ref$V3,]
      
      #generate statistics
      fp_query.stat = calcu_counts_by_variant_type_and_GC(fp_query)
      tp_query.stat = calcu_counts_by_variant_type_and_GC(tp_query)
      fp_ref.stat = calcu_counts_by_variant_type_and_GC(fp_ref)
      tp_ref.stat = calcu_counts_by_variant_type_and_GC(tp_ref)
      #integrate statistics
      query.stat = merge(fp_query.stat, tp_query.stat, by=c("Genomic_Context","SVTYPE","V4"))
      ref.stat = merge(fp_ref.stat, tp_ref.stat, by=c("Genomic_Context","SVTYPE","V4"))
      stat = merge(query.stat, ref.stat, by=c("Genomic_Context","SVTYPE","V4"))
      colnames(stat) = c("Genomic_Context","SVTYPE","SizeRange","fp_query","tp_query","fp_ref","tp_ref")

      write.table(stat, paste(prefix, "stat", sep="."), quote=F, sep="\t", col.names=T, row.names=F)

      pdf(paste(prefix, "pdf", sep="."), height = 4, width = 8)
      par(mfrow=c(1,2))
      par(mar=c(2,3,4,2))
      plot(c(-1.2,1.2),c(0.2,8),frame.plot = F, type = "n", xlab = "", ylab = "", xaxt="n", yaxt="n", main = "(US/RM)")
      axis(2,c(7.5,6:4,3:1-.5), labels = c("SNV","1-30","30-50",">50","1-30","30-50",">50"), las=2,mgp = c(1,.5,0))
      axis(1,c(-5:5)/5, c(5:0,1:5)/5,mgp = c(1,.5,0))
      abline(v =c(-5:5)/5, col="grey")
      abline(v =0, col="black",lwd=2)
      
      add_ovr_bar.US_RM(stat, "SNV","SNV", 7.5)
      add_ovr_bar.US_RM(stat, "Indel_sm","DEL", 6)
      add_ovr_bar.US_RM(stat, "Indel_lg","DEL", 5)
      add_ovr_bar.US_RM(stat, "SV","DEL", 4)
      add_ovr_bar.US_RM(stat, "Indel_sm","INS", 2.5)
      add_ovr_bar.US_RM(stat, "Indel_lg","INS", 1.5)
      add_ovr_bar.US_RM(stat, "SV","INS", 0.5)
      
      plot(c(-1.2,1.2),c(0.2,8),frame.plot = F, type = "n", xlab = "", ylab = "", xaxt="n", yaxt="n", main = "(SD/SR)")
      axis(2,c(7.5,6:4,3:1-.5), labels = c("SNV","1-30","30-50",">50","1-30","30-50",">50"), las=2,mgp = c(1,.5,0))
      axis(1,c(-5:5)/5, c(5:0,1:5)/5,mgp = c(1,.5,0))
      abline(v =c(-5:5)/5, col="grey")
      abline(v =0, col="black",lwd=2)
      rec_width = .3
      add_ovr_bar.SD_SR(stat, "SNV","SNV", 7.5)
      add_ovr_bar.SD_SR(stat, "Indel_sm","DEL", 6)
      add_ovr_bar.SD_SR(stat, "Indel_lg","DEL", 5)
      add_ovr_bar.SD_SR(stat, "SV","DEL", 4)
      add_ovr_bar.SD_SR(stat, "Indel_sm","INS", 2.5)
      add_ovr_bar.SD_SR(stat, "Indel_lg","INS", 1.5)
      add_ovr_bar.SD_SR(stat, "SV","INS", 0.5)
     
      dev.off() 
    }

    plot_ovr_bar("~{fp_query}", "~{tp_query}", "~{fp_ref}", "~{tp_ref}", "~{prefix}")

    '

  >>>


  output {
    File figure = "~{prefix}.pdf"
    File stat = "~{prefix}.stat"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 20 + ceil(size(tp_query,"GiB") + size(tp_ref,"GiB") + size(fp_query,"GiB") + size(fp_ref,"GiB"))*2,
    disk_gb: 25 + ceil(size(tp_query,"GiB") + size(tp_ref,"GiB") + size(fp_query,"GiB") + size(fp_ref,"GiB"))*2,
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

task PlotCompResults{
  input {
    File comp_stat

    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(comp_stat, ".stat")

  command <<<
    set -euxo pipefail

    Rscript -e '

    add_ovr_bar.US_RM<-function(query_vs_ref.stat,size_range, svtype, y_axis){
      rec_width = .3
      rect(-1,y_axis - rec_width, 1, y_axis+rec_width, col=colorBlindGrey8[1])
      snv_usrm = colSums(query_vs_ref.stat[query_vs_ref.stat$SizeRange==size_range & query_vs_ref.stat$SVTYPE==svtype & query_vs_ref.stat$Genomic_Context%in%c("US","RM"),c(4:7)])
      rect(-snv_usrm[2]/sum(snv_usrm[c(1,2)]),y_axis - rec_width, snv_usrm[4]/sum(snv_usrm[c(3,4)]), y_axis+rec_width, col=colorBlindGrey8[2])
      text(-1.2, y_axis+rec_width, format(sum(snv_usrm[c(1,2)]), big.mark = ","), adj=c(0,0))
      text(1.2, y_axis+rec_width, format(sum(snv_usrm[c(3,4)]), big.mark = ","), adj=c(1,0))
    }

    add_ovr_bar.SD_SR<-function(query_vs_ref.stat,size_range, svtype, y_axis){
      rec_width = .3
      rect(-1,y_axis - rec_width, 1, y_axis+rec_width, col=colorBlindGrey8[1])
      snv_usrm = colSums(query_vs_ref.stat[query_vs_ref.stat$SizeRange==size_range & query_vs_ref.stat$SVTYPE==svtype & query_vs_ref.stat$Genomic_Context%in%c("SD","SR"),c(4:7)])
      rect(-snv_usrm[2]/sum(snv_usrm[c(1,2)]),y_axis - rec_width, snv_usrm[4]/sum(snv_usrm[c(3,4)]), y_axis+rec_width, col=colorBlindGrey8[2])
      text(-1.2, y_axis+rec_width, format(sum(snv_usrm[c(1,2)]), big.mark = ","), adj=c(0,0))
      text(1.2, y_axis+rec_width, format(sum(snv_usrm[c(3,4)]), big.mark = ","), adj=c(1,0))
    }

    plot_ovr_bar<-function(stat, prefix){

      pdf(paste(prefix, "pdf", sep="."), height = 4, width = 8)
      par(mfrow=c(1,2))
      par(mar=c(2,3,4,2))
      plot(c(-1.2,1.2),c(0.2,8),frame.plot = F, type = "n", xlab = "", ylab = "", xaxt="n", yaxt="n", main = "(US/RM)")
      axis(2,c(7.5,6:4,3:1-.5), labels = c("SNV","1-30","30-50",">50","1-30","30-50",">50"), las=2,mgp = c(1,.5,0))
      axis(1,c(-5:5)/5, c(5:0,1:5)/5,mgp = c(1,.5,0))
      abline(v =c(-5:5)/5, col="grey")
      abline(v =0, col="black",lwd=2)
      
      add_ovr_bar.US_RM(stat, "SNV","SNV", 7.5)
      add_ovr_bar.US_RM(stat, "Indel_sm","DEL", 6)
      add_ovr_bar.US_RM(stat, "Indel_lg","DEL", 5)
      add_ovr_bar.US_RM(stat, "SV","DEL", 4)
      add_ovr_bar.US_RM(stat, "Indel_sm","INS", 2.5)
      add_ovr_bar.US_RM(stat, "Indel_lg","INS", 1.5)
      add_ovr_bar.US_RM(stat, "SV","INS", 0.5)
      
      plot(c(-1.2,1.2),c(0.2,8),frame.plot = F, type = "n", xlab = "", ylab = "", xaxt="n", yaxt="n", main = "(SD/SR)")
      axis(2,c(7.5,6:4,3:1-.5), labels = c("SNV","1-30","30-50",">50","1-30","30-50",">50"), las=2,mgp = c(1,.5,0))
      axis(1,c(-5:5)/5, c(5:0,1:5)/5,mgp = c(1,.5,0))
      abline(v =c(-5:5)/5, col="grey")
      abline(v =0, col="black",lwd=2)
      rec_width = .3
      add_ovr_bar.SD_SR(stat, "SNV","SNV", 7.5)
      add_ovr_bar.SD_SR(stat, "Indel_sm","DEL", 6)
      add_ovr_bar.SD_SR(stat, "Indel_lg","DEL", 5)
      add_ovr_bar.SD_SR(stat, "SV","DEL", 4)
      add_ovr_bar.SD_SR(stat, "Indel_sm","INS", 2.5)
      add_ovr_bar.SD_SR(stat, "Indel_lg","INS", 1.5)
      add_ovr_bar.SD_SR(stat, "SV","INS", 0.5)
     
      dev.off() 
    }

    colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    stat = read.table("~{comp_stat}", header=T)
    plot_ovr_bar(stat, "~{prefix}")

    '

  >>>


  output {
    File figure = "~{prefix}.pdf"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 5,
    disk_gb: 10,
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

task SplitVcfIntoShards {
  input {
    File input_vcf       # .vcf.gz
    File input_vcf_index # .vcf.gz.tbi
    Int variants_per_shard    # number of varians in each splitted variants
    String output_prefix   # e.g. "chr1.vcf.gz"
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euxo pipefail

    mkdir chunks

    # Extract header
    bcftools view -h ~{input_vcf} > chunks/header.vcf

    # Output body lines (no header), then split
    bcftools view -H ~{input_vcf} | split -l ~{variants_per_shard} - chunks/body_

    # Reconstruct chunked VCFs with header
    for body in chunks/body_*; do
      chunk_name=chunks/~{output_prefix}_$(basename "$body")
      cat chunks/header.vcf "$body" | bgzip -c > "${chunk_name}.vcf.gz"
      tabix -p vcf "${chunk_name}.vcf.gz"
    done
 
   >>>

  output {
    Array[File] split_vcfs = glob("chunks/*.vcf.gz")
    Array[File] split_vcf_indexes = glob("chunks/*.vcf.gz.tbi")
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb:  20,
    disk_gb: 20 + ceil(size(input_vcf,"GiB")*5),
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

task SplitVcfByAnnotationPython {
  input {
    File vcf_file         # Input VCF file (bgzipped or plain)
    File vcf_idx
    File svid_annotation  # 2-column TSV: SVID <tab> annotation
    String docker_image
    String anno_list
    String midfix
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(vcf_file,'.vcf.gz')

  command <<<

    python3 <<CODE
    import pysam
    import os
    from collections import defaultdict

    def determine_svtype(ref, alt):
        if len(ref) == 1 and len(alt) == 1:
            return "SNV"
        elif len(ref) < len(alt):
            return "INS"
        elif len(ref) > len(alt):
            return "DEL"
        else:
            return "OTH"

    def determine_svlen(ref, alt):
        return str(abs(len(alt) - len(ref)))

    # Load SVID -> annotation mapping
    annotation_list = "~{anno_list}".split(',')
    svid_list = []

    print('reading in SVIDs ... ')
    with open("~{svid_annotation}", "r") as f:
        for line in f:
            if line.strip():
                svid, annotation = line.strip().split('\t')
                if annotation in annotation_list:
                  svid_list.append(svid)
    print('complete !')

    # Prepare output VCF writers per annotation
    vcf_in = pysam.VariantFile("~{vcf_file}", "r")
    vcf_out = pysam.VariantFile("~{prefix}.~{midfix}.vcf.gz", 'w', header = vcf_in.header)

    for rec in vcf_in.fetch():
        if rec.id in svid_list:
            print(svid)
            vcf_out.write(rec)

    vcf_out.close()

    CODE

    tabix -p vcf  "~{prefix}.~{midfix}.vcf.gz"
  >>>

  output {
    File split_vcf = "~{prefix}.~{midfix}.vcf.gz"
    File split_vcf_idx = "~{prefix}.~{midfix}.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: 20,
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

task SplitVcfByAnnotationR {
  input {
    File vcf_file         # Input VCF file (bgzipped or plain)
    File vcf_idx
    File svid_annotation  # 2-column TSV: SVID <tab> annotation
    String docker_image
    String anno_list
    String midfix
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(vcf_file,'.vcf.gz')

  command <<<
    set -euxo pipefail

    R --vanilla <<EOF

    annotation_list <- strsplit("~{anno_list}", ',')[[1]]
    svid_gc <- read.table("~{svid_annotation}", header = TRUE)
    svid_list <- svid_gc[svid_gc[,2] %in% annotation_list, 1]
    vcf_in <- read.table("~{vcf_file}", header = FALSE)
    vcf_out <- vcf_in[vcf_in[,3] %in% svid_list, ]
    write.table(vcf_out, file = "~{prefix}.~{midfix}.vcf", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
    
    EOF

    # Re-attach VCF header and compress + index
    
    bcftools view -h ~{vcf_file} | cat - ~{prefix}.~{midfix}.vcf | bgzip > ~{prefix}.~{midfix}.vcf.gz
    tabix -p vcf ~{prefix}.~{midfix}.vcf.gz

  >>>

  output {
    File split_vcf = "~{prefix}.~{midfix}.vcf.gz"
    File split_vcf_idx = "~{prefix}.~{midfix}.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(vcf_file, "GiB")*5),
    disk_gb: 15 + ceil(size(vcf_file, "GiB")*5),
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

task SplitVcfToSites {
  input {
    File vcf_file         # Input VCF file (bgzipped or plain)
    File vcf_idx
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(vcf_file,'.vcf.gz')

  command <<<
    set -euxo pipefail

    zcat ~{vcf_file} | cut -f1-8 | bgzip > "~{prefix}.sites.gz"
    tabix -p vcf "~{prefix}.sites.gz"

  >>>

  output {
    File vcf_sites = "~{prefix}.sites.gz"
    File vcf_sites_idx =  "~{prefix}.sites.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(vcf_file, "GiB")*2),
    disk_gb: 15 + ceil(size(vcf_file, "GiB")*2),
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

task SplitVariantsBySize {
  input {
    File input_vcf
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(input_vcf, ".vcf.gz")
  command <<<
    set -euxo pipefail

    # Index VCF if needed
    if [[ ! -f "~{input_vcf}.tbi" && ! -f "~{input_vcf}.csi" ]]; then
      bcftools index ~{input_vcf}
    fi

    python3 <<CODE
  
    import pysam
    import os

    def get_variant_type_size(rec):
        ref = rec.ref
        alt = rec.alts[0]
        len_ref = len(ref)
        len_alt = len(alt)
        svlen = abs(len_ref - len_alt)
        if "<" in alt and ">" in alt:
          if "SVLEN" in rec.info.keys():
            return "SV_GT_50", rec.info["SVLEN"]
          else:
            return "SV_GT_50", 50
        elif len_ref == 1 and len_alt == 1:
            return "SNV", 0
        elif svlen <= 30:
            return "INDEL_1_30", svlen
        elif 30 < svlen <= 50:
            return "INDEL_30_50", svlen
        else:
            return "SV_GT_50", svlen

    def main(input_vcf):
        # Open input VCF
        vcf_in = pysam.VariantFile(input_vcf)

        # Prepare output files
        outputs = {
            "SNV": pysam.VariantFile(f"snvs.vcf", 'w', header=vcf_in.header),
            "INDEL_1_30": pysam.VariantFile(f"indels_1_30.vcf", 'w', header=vcf_in.header),
            "INDEL_30_50": pysam.VariantFile(f"indels_30_50.vcf", 'w', header=vcf_in.header),
            "SV_GT_50": pysam.VariantFile(f"svs_gt_50.vcf", 'w', header=vcf_in.header)
        }

        for rec in vcf_in.fetch():
            if len(rec.alts) != 1:
                continue  # skip multiallelics for now
            vtype, svlen = get_variant_type_size(rec)
            outputs[vtype].write(rec)

        # Close all files
        for f in outputs.values():
            f.close()
        vcf_in.close()

    main("~{input_vcf}")
    
    CODE

    mv snvs.vcf "~{prefix}.snv.vcf"
    mv indels_1_30.vcf "~{prefix}.indel_1_30.vcf"
    mv indels_30_50.vcf "~{prefix}.indel_31_50.vcf"
    mv svs_gt_50.vcf "~{prefix}.sv_over50.vcf"

    bgzip "~{prefix}.snv.vcf"
    bgzip "~{prefix}.indel_1_30.vcf"
    bgzip "~{prefix}.indel_31_50.vcf"
    bgzip "~{prefix}.sv_over50.vcf"

    tabix -p vcf "~{prefix}.snv.vcf.gz"
    tabix -p vcf "~{prefix}.indel_1_30.vcf.gz"
    tabix -p vcf "~{prefix}.indel_31_50.vcf.gz"
    tabix -p vcf "~{prefix}.sv_over50.vcf.gz"

  >>>

  output {
    File snv_vcf = "~{prefix}.snv.vcf.gz"
    File indel_1_30_vcf = "~{prefix}.indel_1_30.vcf.gz"
    File indel_31_50_vcf = "~{prefix}.indel_31_50.vcf.gz"
    File sv_vcf = "~{prefix}.sv_over50.vcf.gz"

    File snv_vcf_idx = "~{prefix}.snv.vcf.gz.tbi"
    File indel_1_30_vcf_idx = "~{prefix}.indel_1_30.vcf.gz.tbi"
    File indel_31_50_vcf_idx = "~{prefix}.indel_31_50.vcf.gz.tbi"
    File sv_vcf_idx = "~{prefix}.sv_over50.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: 20 + ceil(size(input_vcf, "GiB")*5),
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

task SplitVariantsByType {
  input {
    File input_vcf
    File input_vcf_idx
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(input_vcf, ".vcf.gz")
  command <<<
    set -euxo pipefail

    bcftools view -i 'ID ~ "DEL"' ~{input_vcf} -Oz -o ~{prefix}.DEL.vcf.gz
    bcftools view -i 'ID ~ "INS"' ~{input_vcf} -Oz -o ~{prefix}.INS.vcf.gz
    tabix -p vcf ~{prefix}.DEL.vcf.gz
    tabix -p vcf ~{prefix}.INS.vcf.gz

  >>>

  output {
    File del_vcf = "~{prefix}.DEL.vcf.gz"
    File ins_vcf = "~{prefix}.INS.vcf.gz"
    File del_idx = "~{prefix}.DEL.vcf.gz.tbi"
    File ins_idx = "~{prefix}.INS.vcf.gz.tbi"

  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: 20 + ceil(size(input_vcf, "GiB")*3),
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

task SplitVariantsBySizeAt20bp {
  input {
    File input_vcf
    File input_vcf_idx
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(input_vcf, ".vcf.gz")
  command <<<
    set -euxo pipefail

    # Index VCF if needed
    if [[ ! -f "~{input_vcf}.tbi" && ! -f "~{input_vcf}.csi" ]]; then
      bcftools index ~{input_vcf}
    fi

    python3 <<CODE
  
    import pysam
    import os

    def get_variant_type_size(ref, alt):
        len_ref = len(ref)
        len_alt = len(alt)
        svlen = abs(len_ref - len_alt)
        if "<" in alt and ">" in alt:
          return "SV_GT_20", 0
        elif len_ref == 1 and len_alt == 1:
            return "INDEL_LT_20", 0
        elif svlen < 20:
            return "INDEL_LT_20", svlen
        else:
            return "SV_GT_20", svlen

    def main(input_vcf):
        # Open input VCF
        vcf_in = pysam.VariantFile(input_vcf)

        # Prepare output files
        outputs = {
            "INDEL_LT_20": pysam.VariantFile(f"indels_lt_20.vcf", 'w', header=vcf_in.header),
            "SV_GT_20": pysam.VariantFile(f"svs_gt_20.vcf", 'w', header=vcf_in.header)
        }

        for rec in vcf_in.fetch():
            if len(rec.alts) != 1:
                continue  # skip multiallelics for now
            ref = rec.ref
            alt = rec.alts[0]
            vtype, svlen = get_variant_type_size(ref, alt)
            outputs[vtype].write(rec)

        # Close all files
        for f in outputs.values():
            f.close()
        vcf_in.close()

    main("~{input_vcf}")
    
    CODE

    mv indels_lt_20.vcf "~{prefix}.indels_lt_20.vcf"
    mv svs_gt_20.vcf "~{prefix}.svs_gt_20.vcf"

    bgzip "~{prefix}.indels_lt_20.vcf"
    bgzip "~{prefix}.svs_gt_20.vcf"

    tabix -p vcf "~{prefix}.indels_lt_20.vcf.gz"
    tabix -p vcf "~{prefix}.svs_gt_20.vcf.gz"

  >>>

  output {
    File indels_lt_20_vcf = "~{prefix}.indels_lt_20.vcf.gz"
    File svs_gt_20_vcf = "~{prefix}.svs_gt_20.vcf.gz"

    File indels_lt_20_vcf_idx = "~{prefix}.indels_lt_20.vcf.gz.tbi"
    File svs_gt_20_vcf_idx = "~{prefix}.svs_gt_20.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: 20 + ceil(size(input_vcf, "GiB")*5),
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

task IntegrateInheritanceTable {
  input {
    File inheri_table_snv
    File inheri_table_indel_sm
    File inheri_table_indel_lg
    File inheri_table_sv
    File inheri_table_ref
    String family_id
    String prefix

    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  Int disk_size = 10 
  Int mem_size =  15

  command <<<
    set -euxo pipefail

    Rscript -e '

    process_and_collapse <- function(df) {
      if(nrow(df)>0){
        # Step 1: Replace "." with "0" in columns 2 to 4
        for(i in c(2:4)){
          df[,i] = sapply(df[,i], function(x){strsplit(as.character(x),":")[[1]][1]})
        }
        df = df[df[,2]!="./." | df[,3]!="./." | df[,4]!="./.",]

        df[, 2:4] <- lapply(df[, 2:4], function(col) gsub("/", "|", col))
        df[, 2:4] <- lapply(df[, 2:4], function(col) gsub("\\.", "0", col))
        for(i in c(2:4)){
          if(nrow(df[df[,i]==0,])>0){
            df[df[,i]==0,][,i] = "0|0"  
          }
          if(nrow(df[df[,i]==1,])>0){
            df[df[,i]==1,][,i] = "0|1"  
          }
          if(nrow(df[df[,i]==2,])>0){
            df[df[,i]==2,][,i] = "1|1"  
          }
        }


        # Step 2: Collapse rows by columns 2–4, summing column 1
        # Convert column 1 to numeric (if not already)
        df[, 1] <- as.numeric(df[, 1])

        # Use aggregate to group and sum
        collapsed_df <- aggregate(df[, 1], by = df[, 2:4], FUN = sum)

        # Rename the columns appropriately
        colnames(collapsed_df) <- c("fa", "mo", "pb", "count_variants")

        return(collapsed_df)
      } else {
          df <- data.frame(
            fa = character(),
            mo = character(),
            pb = character(),
            count_variants = character(),
            stringsAsFactors = FALSE
          )
      }
    }


    read_or_empty <- function(file1) {
      if (!file.exists(file1)) {
        stop("File does not exist.")
      }
      
      # Check if file has any lines
      line_count <- length(readLines(file1, warn = FALSE))
      
      if (line_count > 0) {
        df <- read.table(file1, header = FALSE, stringsAsFactors = FALSE)
      } else {
        df <- data.frame(
          V1 = character(),
          V2 = character(),
          V3 = character(),
          V4 = character(),
          stringsAsFactors = FALSE
        )
      }
      
      return(df)
    }

    read_and_merge_five_tables <- function(file1, file2, file3, file4, file5) {
      # Read tables 1-4 without headers
      t1 <- read_or_empty(file1)
      t2 <- read_or_empty(file2)
      t3 <- read_or_empty(file3)
      t4 <- read_or_empty(file4)

      # Process tables 1–4
      p1 <- process_and_collapse(t1)
      p2 <- process_and_collapse(t2)
      p3 <- process_and_collapse(t3)
      p4 <- process_and_collapse(t4)

      # Rename count columns for clarity
      colnames(p1)[4] <- "count_snv"
      colnames(p2)[4] <- "count_indel_sm"
      colnames(p3)[4] <- "count_indel_lg"
      colnames(p4)[4] <- "count_sv"

      # Merge processed tables by col2, col3, col4
      merged <- Reduce(function(x, y) merge(x, y, by = c("fa", "mo", "pb"), all = TRUE),
                       list(p1, p2, p3, p4))

      # Read table 5 with header
      t5 <- read.table(file5, header = TRUE, stringsAsFactors = FALSE)

      # Final merge with table 5
      final <- merge(merged, t5, by.x = c("fa", "mo", "pb"), by.y = c("fa", "mo", "pb"), all = TRUE)

      return(final)
    }

    file1 = "~{inheri_table_snv}"
    file2 = "~{inheri_table_indel_sm}"
    file3 = "~{inheri_table_indel_lg}"
    file4 = "~{inheri_table_sv}"
    file5 = "~{inheri_table_ref}"

    integrated_inheri_table = read_and_merge_five_tables(file1, file2, file3, file4, file5)
    write.table(integrated_inheri_table, "~{family_id}.~{prefix}.integrated_inheritance_table.tsv", quote=F, sep="\t", col.names=T, row.names=F)
    
    '

  >>>

  output {
    File integrated_inheri_stat = "~{family_id}.~{prefix}.integrated_inheritance_table.tsv"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: mem_size,
    disk_gb: disk_size,
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

task IntegrateInheriByGQTable {
  input {
    File inheri_gq_table_snv
    File inheri_gq_table_indel_sm
    File inheri_gq_table_indel_lg
    File inheri_gq_table_sv
    String family_id
    String prefix

    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  Int disk_size = 10 
  Int mem_size =  15

  command <<<
    set -euxo pipefail

    sed -e 's/$/\tSNV/' ~{inheri_gq_table_snv} > "~{family_id}.~{prefix}.inheri_by_gq.stat"
    sed -e 's/$/\tIndel_sm/' ~{inheri_gq_table_indel_sm} >> "~{family_id}.~{prefix}.inheri_by_gq.stat"
    sed -e 's/$/\tIndel_lg/' ~{inheri_gq_table_indel_lg} >> "~{family_id}.~{prefix}.inheri_by_gq.stat"
    sed -e 's/$/\tSV/' ~{inheri_gq_table_sv} >> "~{family_id}.~{prefix}.inheri_by_gq.stat"

  >>>

  output {
    File integrated_inheri_by_gq_stat = "~{family_id}.~{prefix}.inheri_by_gq.stat"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: mem_size,
    disk_gb: disk_size,
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

task WriteTrioSampleFile {
  input {
    Array[String] family  # [father, mother, child]
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command {
    echo -e "${family[0]}\n${family[1]}\n${family[2]}" > trio_${family[0]}_${family[2]}.samples.txt
  }

  output {
    File sample_file = "trio_${family[0]}_${family[2]}.samples.txt"
    String family_id = "${family[0]}_${family[2]}"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: 5,
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

task IndexPanGenieRefPanel {
    input {
        File panel_vcf_gz
        File panel_vcf_gz_tbi
        File reference_fasta
        Array[String] chromosomes
        String output_prefix

        String docker_image
        File? monitoring_script
        Int? kmer_length = 31
        String? extra_args

        RuntimeAttr? runtime_attr_override
    }

    Int num_chromosomes = length(chromosomes)

    command {
        set -euxo pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # subset reference and panel VCF to chromosomes
        samtools faidx -r <(echo -e "~{sep="\n" chromosomes}") ~{reference_fasta} > reference.subset.fa
        bcftools view ~{panel_vcf_gz} -r ~{sep="," chromosomes} > panel.subset.vcf

        NPROC=$(nproc)
        # NUM_CHROMOSOMES=~{num_chromosomes}
        # NUM_THREADS=$(( NPROC < NUM_CHROMOSOMES ? NPROC : NUM_CHROMOSOMES ))
        NUM_THREADS=$NPROC

        /pangenie/build/src/PanGenie-index \
            -o ~{output_prefix} \
            -r reference.subset.fa \
            -t $NUM_THREADS \
            -v panel.subset.vcf \
            -k ~{kmer_length} \
            ~{extra_args}
    }

    RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 100,
      disk_gb: 100 + ceil(size(panel_vcf_gz, "GiB")*3),
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

    output {
        File monitoring_log = "monitoring.log"
        Array[File] pangenie_index_chromosome_graphs = glob("~{output_prefix}_*_Graph.cereal")
        Array[File] pangenie_index_chromosome_kmers = glob("~{output_prefix}_*_kmers.tsv.gz")
        File pangenie_index_unique_kmers_map = "~{output_prefix}_UniqueKmersMap.cereal"
        File pangenie_index_path_segments_fasta = "~{output_prefix}_path_segments.fasta"
    }
}

task IndexPanGenieCaseReads {
    input {
        File input_cram

        String docker_image
        File? monitoring_script

        RuntimeAttr? runtime_attr_override
    }

    # if we instead simply use File cram_idx = "~{input_cram}.crai" in the output block,
    # Terra tries to localize an index adjacent to the CRAM in PreprocessCaseReads???
    String output_prefix = basename(input_cram, ".cram")

    command {
        set -euxo pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        samtools index -@ $(nproc) ~{input_cram} ~{output_prefix}.cram.crai
    }

    RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 2,
      disk_gb: 5,
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

    output {
        File monitoring_log = "monitoring.log"
        File cram_idx = "~{output_prefix}.cram.crai"
    }
}

task PreprocessBiallelicRefPanelVcf{
    input{
      File input_vcf
      File? input_vcf_idx

      String docker_image
      File? monitoring_script

      RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(input_vcf, ".vcf.gz")

    command <<<
      set -euxo pipefail

      # Index VCF if needed
      if [[ ! -f "~{input_vcf}.tbi" && ! -f "~{input_vcf}.csi" ]]; then
        bcftools index ~{input_vcf}
      fi

      python3 <<CODE

      import pysam
      import sys

      def extract_info_id_to_record_id(input_vcf, output_vcf):
          vcf_in = pysam.VariantFile(input_vcf, "r")

          vcf_out = pysam.VariantFile(output_vcf, "w", header=vcf_in.header)

          for rec in vcf_in:
              # Get ID from INFO column
              info_id = rec.info.get("ID", None)
              if info_id:
                  # Set record.id to the info value
                  rec.id = str(info_id[0])
              vcf_out.write(rec)

          vcf_in.close()
          vcf_out.close()
          print(f"Updated VCF written to: {output_vcf}")

      extract_info_id_to_record_id("~{input_vcf}", "~{prefix}.preprocessed.vcf.gz")    

      CODE

      tabix -p vcf "~{prefix}.preprocessed.vcf.gz"

    >>>

    RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 5,
      disk_gb: 10 + ceil(size(input_vcf, "GiB")*2),
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

    output {
      File preprocessed_vcf = "~{prefix}.preprocessed.vcf.gz"
      File preprocessed_vcf_idx = "~{prefix}.preprocessed.vcf.gz.tbi"
    }
}

task PreprocessPanGenieCaseReads {
    input {
        File input_cram
        File input_cram_idx
        File reference_fasta
        File reference_fasta_fai
        Array[String] chromosomes
        String output_prefix

        String docker_image
        File? monitoring_script

        RuntimeAttr? runtime_attr_override
    }

    String filter_N_regex = "/^>/{N;/^>.*\\n.*N.*/d}"

    command {
        set -euxo pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # hacky way to get chromosomes into bed file
        grep -P '~{sep="\\t|" chromosomes}\t' ~{reference_fasta_fai} | cut -f 1,2 | sed -e 's/\t/\t1\t/g' > chromosomes.bed

        # filter out read pairs containing N nucleotides
        # TODO move functionality into KAGE code
        samtools view --reference ~{reference_fasta} -@ $(nproc) -L chromosomes.bed -u ~{input_cram} | \
            samtools fasta --reference ~{reference_fasta} -@ $(nproc) | \
            sed -E '~{filter_N_regex}' > ~{output_prefix}.preprocessed.fasta
    }

    RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 30,
      disk_gb: 500,
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
    output {
        File monitoring_log = "monitoring.log"
        File preprocessed_fasta = "~{output_prefix}.preprocessed.fasta"
    }
}

task PreprocessPanGenieCaseReadsWithoutSubsetting {
    input {
        File input_cram
        File input_cram_idx
        File reference_fasta
        File reference_fasta_fai
        String output_prefix

        String docker_image
        File? monitoring_script

        RuntimeAttr? runtime_attr_override
    }

    String filter_N_regex = "/^>/{N;/^>.*\\n.*N.*/d}"

    command {
        set -euxo pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # filter out read pairs containing N nucleotides
        # TODO move functionality into KAGE code
        samtools fasta --reference ~{reference_fasta} -@ $(nproc) ~{input_cram} | sed -E '~{filter_N_regex}' > ~{output_prefix}.preprocessed.fasta
    }

    RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 30,
      disk_gb: 500,
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

    output {
        File monitoring_log = "monitoring.log"
        File preprocessed_fasta = "~{output_prefix}.preprocessed.fasta"
    }
}

task PanGenieGenotype {
    input {
        Array[File] pangenie_index_chromosome_graphs
        Array[File] pangenie_index_chromosome_kmers
        File pangenie_index_unique_kmers_map
        File pangenie_index_path_segments_fasta
        String index_prefix
        File input_fasta
        String sample_name
        String output_prefix

        String docker_image
        File? monitoring_script
        String? extra_args

        RuntimeAttr? runtime_attr_override
    }

    command {
        set -euxo pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        NPROC=$(nproc)
        NUM_THREADS=$NPROC
        NUM_THREADS=32

        mv ~{sep=" " pangenie_index_chromosome_graphs} .
        mv ~{sep=" " pangenie_index_chromosome_kmers} .
        mv ~{pangenie_index_unique_kmers_map} .
        mv ~{pangenie_index_path_segments_fasta} .

        /pangenie/build/src/PanGenie \
            -i ~{input_fasta} \
            -f ~{index_prefix} \
            -o ~{output_prefix} \
            -s ~{sample_name} \
            -t $NUM_THREADS \
            -j $NUM_THREADS \
            ~{extra_args}

        bgzip -c ~{output_prefix}_genotyping.vcf > ~{output_prefix}_genotyping.vcf.gz
        bcftools index -t ~{output_prefix}_genotyping.vcf.gz

        # naively set all GTs to phased for Vcfdist evaluation
        bcftools +setGT ~{output_prefix}_genotyping.vcf.gz -Oz -o ~{output_prefix}_genotyping_naively_phased.vcf.gz -- -t a -n p
        bcftools index -t ~{output_prefix}_genotyping_naively_phased.vcf.gz
    }

    RuntimeAttr default_attr = object {
      cpu_cores: 1,
      mem_gb: 50,
      disk_gb: ceil(size(input_fasta, "GiB")*2),
      boot_disk_gb: 10,
      preemptible_tries: 1,
      max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
      cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
      memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
      disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
      bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
      docker: docker_image
      preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
      maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        File monitoring_log = "monitoring.log"
        File genotyping_vcf_gz = "~{output_prefix}_genotyping.vcf.gz"
        File genotyping_vcf_gz_tbi = "~{output_prefix}_genotyping.vcf.gz.tbi"
        File genotyping_naively_phased_vcf_gz = "~{output_prefix}_genotyping_naively_phased.vcf.gz"
        File genotyping_naively_phased_vcf_gz_tbi = "~{output_prefix}_genotyping_naively_phased.vcf.gz.tbi"
        File histogram = "~{output_prefix}_histogram.histo"
        #File path_segments_fasta = "~{output_prefix}_path_segments.fasta"
    }
}

struct RuntimeAttributes {
    Int? cpu
    Int? command_mem_gb
    Int? additional_mem_gb
    Int? disk_size_gb
    Int? boot_disk_size_gb
    Boolean? use_ssd
    Int? preemptible
    Int? max_retries
}

