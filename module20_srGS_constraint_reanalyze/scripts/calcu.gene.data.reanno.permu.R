#Process SNV gene data
read.snvdata <- function(SNVdata.in, gene.metadata.in, gene.utr.in, phaplo_ptriplo_file){
  #Read & clean
  snv.data <- read.table(SNVdata.in, header=T, sep=',') 
  #snv.data = snv.data[!is.na(snv.data$LOEUF) & !is.na(snv.data$pHaplo) & !is.na(snv.data$pTriplo),]
  
  colnames(snv.data) [1] = 'gene'
  metadata <- read.table(gene.metadata.in, header=T, comment.char="")  # comment.char=""
  merged <- merge(x=snv.data, y=metadata,  by="gene", sort=F)
  merged <- merged[which(!(merged$X.chr %in% c("chrX", "chrY"))), ]
  
  utr = read.table(gene.utr.in, header = T)
  utr[,1] = gsub(';','', utr[,1])
  colnames(utr)[5] = 'gene'
  merged <-merge(merged, utr, by='gene')
  
  
  #Assign oe deciles
  #   merged$mis_oe_dec <- ceiling(10*rank(merged$oe_mis_upper)/(nrow(merged)+1))  # mis_oe --> oe_mis_upper
  merged$ptv_oe_dec <- ceiling(10*rank(merged$LOEUF)/(nrow(merged)+1))  # ptv_oe --> LOEUF
  #Assign oe to 40 bins
  #   merged$mis_oe_binrank <- ceiling(40*rank(merged$oe_mis_upper)/(nrow(merged)+1))
  merged$ptv_oe_binrank <- ceiling(40*rank(merged$LOEUF)/(nrow(merged)+1))
  #Assign oe percentiles
  #   merged$mis_oe_cent <- ceiling(100*rank(merged$oe_mis_upper)/(nrow(merged)+1))
  merged$ptv_oe_cent <- ceiling(100*rank(merged$LOEUF)/(nrow(merged)+1))
  #Return formatted data
  
  
  phaplo_ptriplo = read.table(phaplo_ptriplo_file, header = T, comment.char = "")
  colnames(phaplo_ptriplo)[1] = 'gene'
  snv.data.2 = merge_loeuf_phaplo_ptriplo(merged[,-6], phaplo_ptriplo)
  return(snv.data.2)

 }

readin.re_anno_svtype<-function(re_anno_SV_file, sv_info){
  re_anno = read.table(re_anno_SV_file, header = T)
  re_anno = re_anno[re_anno[,1]!="SVID",]
  colnames(re_anno)[1] = 'name'
  re_anno[,ncol(re_anno)+1] = paste(re_anno$X5_prime_utr,re_anno$X3_prime_utr, sep = ',')
  colnames(re_anno)[ncol(re_anno)] = 'utr'
  re_anno$utr = gsub('NA,','', re_anno$utr)
  re_anno$utr = gsub(',NA','', re_anno$utr)
  re_anno[re_anno$utr=="NA",]$utr = NA
  re_anno=merge(re_anno, sv_info[,c("name","SVTYPE")], by='name')
  return(re_anno)
}

permutate.snvdata<-function(snv.data){
  snv.data.out = snv.data
  for(i in c(1:100)){
    print(i)
    tmp = snv.data
    tmp$gene = paste(tmp$gene, '.permi_', i, sep='')
    snv.data.out = rbind(snv.data.out, tmp)
  }
  return(snv.data.out)
}

getSVdat.reannotated.permutate <-function(dat,re_anno_svtype, genes, prefix=NULL){
  re_anno_SVs = re_anno_svtype[re_anno_svtype$name%in%dat$name,]
  #whole transcript overlap
  whole_transcript_overlap.any <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$whole_transcript_overlap),]$whole_transcript_overlap), split = ',')))
  whole_transcript_overlap.del <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$whole_transcript_overlap) & re_anno_SVs$SVTYPE=="DEL",]$whole_transcript_overlap), split = ',')))
  whole_transcript_overlap.dup <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$whole_transcript_overlap) & re_anno_SVs$SVTYPE=="DUP",]$whole_transcript_overlap), split = ',')))
  whole_transcript_overlap.ins <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$whole_transcript_overlap) & re_anno_SVs$SVTYPE=="INS",]$whole_transcript_overlap), split = ',')))
  whole_transcript_overlap.inv <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$whole_transcript_overlap) & re_anno_SVs$SVTYPE=="INV",]$whole_transcript_overlap), split = ',')))
  
  #intact exon overlap
  intact_exon_overlap.any <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$intact_exon_overlap),]$intact_exon_overlap), split = ',')))
  intact_exon_overlap.del <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$intact_exon_overlap) & re_anno_SVs$SVTYPE=="DEL",]$intact_exon_overlap), split = ',')))
  intact_exon_overlap.dup <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$intact_exon_overlap) & re_anno_SVs$SVTYPE=="DUP",]$intact_exon_overlap), split = ',')))
  intact_exon_overlap.ins <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$intact_exon_overlap) & re_anno_SVs$SVTYPE=="INS",]$intact_exon_overlap), split = ',')))
  intact_exon_overlap.inv <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$intact_exon_overlap) & re_anno_SVs$SVTYPE=="INV",]$intact_exon_overlap), split = ',')))
  
  #partial exon overlap
  partial_exon_overlap.any <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$partial_exon_overlap),]$partial_exon_overlap), split = ',')))
  partial_exon_overlap.del <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$partial_exon_overlap) & re_anno_SVs$SVTYPE=="DEL",]$partial_exon_overlap), split = ',')))
  partial_exon_overlap.dup <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$partial_exon_overlap) & re_anno_SVs$SVTYPE=="DUP",]$partial_exon_overlap), split = ',')))
  partial_exon_overlap.ins <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$partial_exon_overlap) & re_anno_SVs$SVTYPE=="INS",]$partial_exon_overlap), split = ',')))
  partial_exon_overlap.inv <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$partial_exon_overlap) & re_anno_SVs$SVTYPE=="INV",]$partial_exon_overlap), split = ',')))
  
  #tss transcript overlap
  tss_transcripts_overlap.any <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$tss_transcripts_overlap),]$tss_transcripts_overlap), split = ',')))
  tss_transcripts_overlap.del <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$tss_transcripts_overlap) & re_anno_SVs$SVTYPE=="DEL",]$tss_transcripts_overlap), split = ',')))
  tss_transcripts_overlap.dup <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$tss_transcripts_overlap) & re_anno_SVs$SVTYPE=="DUP",]$tss_transcripts_overlap), split = ',')))
  tss_transcripts_overlap.ins <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$tss_transcripts_overlap) & re_anno_SVs$SVTYPE=="INS",]$tss_transcripts_overlap), split = ',')))
  tss_transcripts_overlap.inv <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$tss_transcripts_overlap) & re_anno_SVs$SVTYPE=="INV",]$tss_transcripts_overlap), split = ',')))
  
  #partial transcript overlap
  partial_transcripts_overlap.any <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$partial_transcripts_overlap),]$partial_transcripts_overlap), split = ',')))
  partial_transcripts_overlap.del <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$partial_transcripts_overlap) & re_anno_SVs$SVTYPE=="DEL",]$partial_transcripts_overlap), split = ',')))
  partial_transcripts_overlap.dup <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$partial_transcripts_overlap) & re_anno_SVs$SVTYPE=="DUP",]$partial_transcripts_overlap), split = ',')))
  partial_transcripts_overlap.ins <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$partial_transcripts_overlap) & re_anno_SVs$SVTYPE=="INS",]$partial_transcripts_overlap), split = ',')))
  partial_transcripts_overlap.inv <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$partial_transcripts_overlap) & re_anno_SVs$SVTYPE=="INV",]$partial_transcripts_overlap), split = ',')))

  #tss transcript overlap
  tss_coding_transcripts_overlap.any <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$tss_coding_transcripts_overlap),]$tss_coding_transcripts_overlap), split = ',')))
  tss_coding_transcripts_overlap.del <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$tss_coding_transcripts_overlap) & re_anno_SVs$SVTYPE=="DEL",]$tss_coding_transcripts_overlap), split = ',')))
  tss_coding_transcripts_overlap.dup <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$tss_coding_transcripts_overlap) & re_anno_SVs$SVTYPE=="DUP",]$tss_coding_transcripts_overlap), split = ',')))
  tss_coding_transcripts_overlap.ins <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$tss_coding_transcripts_overlap) & re_anno_SVs$SVTYPE=="INS",]$tss_coding_transcripts_overlap), split = ',')))
  tss_coding_transcripts_overlap.inv <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$tss_coding_transcripts_overlap) & re_anno_SVs$SVTYPE=="INV",]$tss_coding_transcripts_overlap), split = ',')))
  
  #partial transcript overlap
  partial_coding_transcripts_overlap.any <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$partial_coding_transcripts_overlap),]$partial_coding_transcripts_overlap), split = ',')))
  partial_coding_transcripts_overlap.del <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$partial_coding_transcripts_overlap) & re_anno_SVs$SVTYPE=="DEL",]$partial_coding_transcripts_overlap), split = ',')))
  partial_coding_transcripts_overlap.dup <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$partial_coding_transcripts_overlap) & re_anno_SVs$SVTYPE=="DUP",]$partial_coding_transcripts_overlap), split = ',')))
  partial_coding_transcripts_overlap.ins <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$partial_coding_transcripts_overlap) & re_anno_SVs$SVTYPE=="INS",]$partial_coding_transcripts_overlap), split = ',')))
  partial_coding_transcripts_overlap.inv <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$partial_coding_transcripts_overlap) & re_anno_SVs$SVTYPE=="INV",]$partial_coding_transcripts_overlap), split = ',')))
  
  #5' overlap
  X5_prime_utr.any <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$X5_prime_utr),]$X5_prime_utr), split = ',')))
  X5_prime_utr.del <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$X5_prime_utr) & re_anno_SVs$SVTYPE=="DEL",]$X5_prime_utr), split = ',')))
  X5_prime_utr.dup <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$X5_prime_utr) & re_anno_SVs$SVTYPE=="DUP",]$X5_prime_utr), split = ',')))
  X5_prime_utr.ins <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$X5_prime_utr) & re_anno_SVs$SVTYPE=="INS",]$X5_prime_utr), split = ',')))
  X5_prime_utr.inv <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$X5_prime_utr) & re_anno_SVs$SVTYPE=="INV",]$X5_prime_utr), split = ',')))
  
  #3' overlap
  X3_prime_utr.any <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$X3_prime_utr),]$X3_prime_utr), split = ',')))
  X3_prime_utr.del <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$X3_prime_utr) & re_anno_SVs$SVTYPE=="DEL",]$X3_prime_utr), split = ',')))
  X3_prime_utr.dup <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$X3_prime_utr) & re_anno_SVs$SVTYPE=="DUP",]$X3_prime_utr), split = ',')))
  X3_prime_utr.ins <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$X3_prime_utr) & re_anno_SVs$SVTYPE=="INS",]$X3_prime_utr), split = ',')))
  X3_prime_utr.inv <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$X3_prime_utr) & re_anno_SVs$SVTYPE=="INV",]$X3_prime_utr), split = ',')))
  
  #utr overlap regardless of direction
  utr.any <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$utr),]$utr), split = ',')))
  utr.del <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$utr) & re_anno_SVs$SVTYPE=="DEL",]$utr), split = ',')))
  utr.dup <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$utr) & re_anno_SVs$SVTYPE=="DUP",]$utr), split = ',')))
  utr.ins <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$utr) & re_anno_SVs$SVTYPE=="INS",]$utr), split = ',')))
  utr.inv <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$utr) & re_anno_SVs$SVTYPE=="INV",]$utr), split = ',')))
  
  #SVs within exon
  inside_exons.any <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$inside_exons),]$inside_exons), split = ',')))
  inside_exons.del <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$inside_exons) & re_anno_SVs$SVTYPE=="DEL",]$inside_exons), split = ',')))
  inside_exons.dup <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$inside_exons) & re_anno_SVs$SVTYPE=="DUP",]$inside_exons), split = ',')))
  inside_exons.ins <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$inside_exons) & re_anno_SVs$SVTYPE=="INS",]$inside_exons), split = ',')))
  inside_exons.inv <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$inside_exons) & re_anno_SVs$SVTYPE=="INV",]$inside_exons), split = ',')))

  #SVs within introns
  inside_introns.any <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$inside_introns),]$inside_introns), split = ',')))
  inside_introns.del <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$inside_introns) & re_anno_SVs$SVTYPE=="DEL",]$inside_introns), split = ',')))
  inside_introns.dup <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$inside_introns) & re_anno_SVs$SVTYPE=="DUP",]$inside_introns), split = ',')))
  inside_introns.ins <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$inside_introns) & re_anno_SVs$SVTYPE=="INS",]$inside_introns), split = ',')))
  inside_introns.inv <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$inside_introns) & re_anno_SVs$SVTYPE=="INV",]$inside_introns), split = ',')))

  #SVs within exon
  promoter_overlap.any <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$promoter_overlap),]$promoter_overlap), split = ',')))
  promoter_overlap.del <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$promoter_overlap) & re_anno_SVs$SVTYPE=="DEL",]$promoter_overlap), split = ',')))
  promoter_overlap.dup <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$promoter_overlap) & re_anno_SVs$SVTYPE=="DUP",]$promoter_overlap), split = ',')))
  promoter_overlap.ins <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$promoter_overlap) & re_anno_SVs$SVTYPE=="INS",]$promoter_overlap), split = ',')))
  promoter_overlap.inv <- as.character(unlist(strsplit(as.character(re_anno_SVs[!is.na(re_anno_SVs$promoter_overlap) & re_anno_SVs$SVTYPE=="INV",]$promoter_overlap), split = ',')))
  
  #Collect vector per gene
  res <- as.data.frame(t(sapply(genes, function(gene){
    g.whole_transcript_overlap.any <- length(which(whole_transcript_overlap.any==gene))
    g.whole_transcript_overlap.del <- length(which(whole_transcript_overlap.del==gene))
    g.whole_transcript_overlap.dup <- length(which(whole_transcript_overlap.dup==gene))
    g.whole_transcript_overlap.ins <- length(which(whole_transcript_overlap.ins==gene))
    g.whole_transcript_overlap.inv <- length(which(whole_transcript_overlap.inv==gene))
    
    g.tss_transcripts_overlap.any <- length(which(tss_transcripts_overlap.any==gene))
    g.tss_transcripts_overlap.del <- length(which(tss_transcripts_overlap.del==gene))
    g.tss_transcripts_overlap.dup <- length(which(tss_transcripts_overlap.dup==gene))
    g.tss_transcripts_overlap.ins <- length(which(tss_transcripts_overlap.ins==gene))
    g.tss_transcripts_overlap.inv <- length(which(tss_transcripts_overlap.inv==gene))
    
    g.partial_transcripts_overlap.any <- length(which(partial_transcripts_overlap.any==gene))
    g.partial_transcripts_overlap.del <- length(which(partial_transcripts_overlap.del==gene))
    g.partial_transcripts_overlap.dup <- length(which(partial_transcripts_overlap.dup==gene))
    g.partial_transcripts_overlap.ins <- length(which(partial_transcripts_overlap.ins==gene))
    g.partial_transcripts_overlap.inv <- length(which(partial_transcripts_overlap.inv==gene))
    
    g.tss_coding_transcripts_overlap.any <- length(which(tss_coding_transcripts_overlap.any==gene))
    g.tss_coding_transcripts_overlap.del <- length(which(tss_coding_transcripts_overlap.del==gene))
    g.tss_coding_transcripts_overlap.dup <- length(which(tss_coding_transcripts_overlap.dup==gene))
    g.tss_coding_transcripts_overlap.ins <- length(which(tss_coding_transcripts_overlap.ins==gene))
    g.tss_coding_transcripts_overlap.inv <- length(which(tss_coding_transcripts_overlap.inv==gene))
    
    g.partial_coding_transcripts_overlap.any <- length(which(partial_coding_transcripts_overlap.any==gene))
    g.partial_coding_transcripts_overlap.del <- length(which(partial_coding_transcripts_overlap.del==gene))
    g.partial_coding_transcripts_overlap.dup <- length(which(partial_coding_transcripts_overlap.dup==gene))
    g.partial_coding_transcripts_overlap.ins <- length(which(partial_coding_transcripts_overlap.ins==gene))
    g.partial_coding_transcripts_overlap.inv <- length(which(partial_coding_transcripts_overlap.inv==gene))
    
    g.intact_exon_overlap.any <- length(which(intact_exon_overlap.any==gene))
    g.intact_exon_overlap.del <- length(which(intact_exon_overlap.del==gene))
    g.intact_exon_overlap.dup <- length(which(intact_exon_overlap.dup==gene))
    g.intact_exon_overlap.ins <- length(which(intact_exon_overlap.ins==gene))
    g.intact_exon_overlap.inv <- length(which(intact_exon_overlap.inv==gene))
    
    g.partial_exon_overlap.any <- length(which(partial_exon_overlap.any==gene))
    g.partial_exon_overlap.del <- length(which(partial_exon_overlap.del==gene))
    g.partial_exon_overlap.dup <- length(which(partial_exon_overlap.dup==gene))
    g.partial_exon_overlap.ins <- length(which(partial_exon_overlap.ins==gene))
    g.partial_exon_overlap.inv <- length(which(partial_exon_overlap.inv==gene))
    
    g.X5_prime_utr.any <- length(which(X5_prime_utr.any==gene))
    g.X5_prime_utr.del <- length(which(X5_prime_utr.del==gene))
    g.X5_prime_utr.dup <- length(which(X5_prime_utr.dup==gene))
    g.X5_prime_utr.ins <- length(which(X5_prime_utr.ins==gene))
    g.X5_prime_utr.inv <- length(which(X5_prime_utr.inv==gene))
    
    g.X3_prime_utr.any <- length(which(X3_prime_utr.any==gene))
    g.X3_prime_utr.del <- length(which(X3_prime_utr.del==gene))
    g.X3_prime_utr.dup <- length(which(X3_prime_utr.dup==gene))
    g.X3_prime_utr.ins <- length(which(X3_prime_utr.ins==gene))
    g.X3_prime_utr.inv <- length(which(X3_prime_utr.inv==gene))
    
    g.utr.any <- length(which(utr.any==gene))
    g.utr.del <- length(which(utr.del==gene))
    g.utr.dup <- length(which(utr.dup==gene))
    g.utr.ins <- length(which(utr.ins==gene))
    g.utr.inv <- length(which(utr.inv==gene))
    
    g.inside_exons.any <- length(which(inside_exons.any==gene))
    g.inside_exons.del <- length(which(inside_exons.del==gene))
    g.inside_exons.dup <- length(which(inside_exons.dup==gene))
    g.inside_exons.ins <- length(which(inside_exons.ins==gene))
    g.inside_exons.inv <- length(which(inside_exons.inv==gene))

    g.inside_introns.any <- length(which(inside_introns.any==gene))
    g.inside_introns.del <- length(which(inside_introns.del==gene))
    g.inside_introns.dup <- length(which(inside_introns.dup==gene))
    g.inside_introns.ins <- length(which(inside_introns.ins==gene))
    g.inside_introns.inv <- length(which(inside_introns.inv==gene))

    g.promoter_overlap.any <- length(which(promoter_overlap.any==gene))
    g.promoter_overlap.del <- length(which(promoter_overlap.del==gene))
    g.promoter_overlap.dup <- length(which(promoter_overlap.dup==gene))
    g.promoter_overlap.ins <- length(which(promoter_overlap.ins==gene))
    g.promoter_overlap.inv <- length(which(promoter_overlap.inv==gene))
    
    g.out <- as.integer(c(
      g.whole_transcript_overlap.any, g.whole_transcript_overlap.del, g.whole_transcript_overlap.dup, g.whole_transcript_overlap.ins, g.whole_transcript_overlap.inv,
      g.tss_transcripts_overlap.any,  g.tss_transcripts_overlap.del,  g.tss_transcripts_overlap.dup,  g.tss_transcripts_overlap.ins,  g.tss_transcripts_overlap.inv, 
      g.partial_transcripts_overlap.any,g.partial_transcripts_overlap.del,g.partial_transcripts_overlap.dup,g.partial_transcripts_overlap.ins,g.partial_transcripts_overlap.inv,     
      g.tss_coding_transcripts_overlap.any,  g.tss_coding_transcripts_overlap.del,  g.tss_coding_transcripts_overlap.dup,  g.tss_coding_transcripts_overlap.ins,  g.tss_coding_transcripts_overlap.inv, 
      g.partial_coding_transcripts_overlap.any,g.partial_coding_transcripts_overlap.del,g.partial_coding_transcripts_overlap.dup,g.partial_coding_transcripts_overlap.ins,g.partial_coding_transcripts_overlap.inv,     
      g.intact_exon_overlap.any,      g.intact_exon_overlap.del,      g.intact_exon_overlap.dup,      g.intact_exon_overlap.ins,      g.intact_exon_overlap.inv,
      g.partial_exon_overlap.any,     g.partial_exon_overlap.del,     g.partial_exon_overlap.dup,     g.partial_exon_overlap.ins,     g.partial_exon_overlap.inv,     
      g.X5_prime_utr.any,             g.X5_prime_utr.del,             g.X5_prime_utr.dup,             g.X5_prime_utr.ins,             g.X5_prime_utr.inv,             
      g.X3_prime_utr.any,             g.X3_prime_utr.del,             g.X3_prime_utr.dup,             g.X3_prime_utr.ins,             g.X3_prime_utr.inv,             
      g.inside_exons.any,             g.inside_exons.del,             g.inside_exons.dup,             g.inside_exons.ins,             g.inside_exons.inv,  
      g.inside_introns.any,           g.inside_introns.del,           g.inside_introns.dup,           g.inside_introns.ins,           g.inside_introns.inv,  
      g.promoter_overlap.any,         g.promoter_overlap.del,         g.promoter_overlap.dup,         g.promoter_overlap.ins,         g.promoter_overlap.inv,  
      g.utr.any,             g.utr.del,             g.utr.dup,             g.utr.ins,             g.utr.inv))
    g.out[which(is.na(g.out))] <- 0
    g.out <- c(gene, g.out)
    return(g.out)
    
  })))  
  
  colnames(res) <- c("gene", 
                     'whole_transcript_overlap.any', 'whole_transcript_overlap.del','whole_transcript_overlap.dup','whole_transcript_overlap.ins','whole_transcript_overlap.inv',
                     'tss_transcripts_overlap.any',  'tss_transcripts_overlap.del',  'tss_transcripts_overlap.dup',  'tss_transcripts_overlap.ins',  'tss_transcripts_overlap.inv', 
                     'partial_transcripts_overlap.any','partial_transcripts_overlap.del','partial_transcripts_overlap.dup','partial_transcripts_overlap.ins','partial_transcripts_overlap.inv',     
                     'tss_coding_transcripts_overlap.any',  'tss_coding_transcripts_overlap.del',  'tss_coding_transcripts_overlap.dup',  'tss_coding_transcripts_overlap.ins',  'tss_coding_transcripts_overlap.inv', 
                     'partial_coding_transcripts_overlap.any','partial_coding_transcripts_overlap.del','partial_coding_transcripts_overlap.dup','partial_coding_transcripts_overlap.ins','partial_coding_transcripts_overlap.inv',     
                     'intact_exon_overlap.any', 'intact_exon_overlap.del','intact_exon_overlap.dup','intact_exon_overlap.ins','intact_exon_overlap.inv',
                     'partial_exon_overlap.any', 'partial_exon_overlap.del','partial_exon_overlap.dup','partial_exon_overlap.ins','partial_exon_overlap.inv',
                     'X5_prime_utr.any', 'X5_prime_utr.del','X5_prime_utr.dup','X5_prime_utr.ins','X5_prime_utr.inv',
                     'X3_prime_utr.any', 'X3_prime_utr.del','X3_prime_utr.dup','X3_prime_utr.ins','X3_prime_utr.inv',
                     'inside_exons.any', 'inside_exons.del','inside_exons.dup','inside_exons.ins','inside_exons.inv',
                     'inside_introns.any', 'inside_introns.del','inside_introns.dup','inside_introns.ins','inside_introns.inv',
                     'promoter_overlap.any', 'promoter_overlap.del','promoter_overlap.dup','promoter_overlap.ins','promoter_overlap.inv',
                      'utr.any', 'utr.del','utr.dup','utr.ins','utr.inv')
  if(!is.null(prefix)){    colnames(res)[-1] <- paste(prefix, colnames(res)[-1], sep=".")   }
  rownames(res) <- 1:nrow(res)
  return(res)
}


getSVdat.all.reannotated.permutate <- function(dat, re_anno_svtype, snv.data, include.cpx=T, require.SR=F){
  #Gather data
  genes <- sort(unique(as.character(snv.data$gene)))
  if(include.cpx==F){
    dat <- dat[which(dat$SVTYPE != "CPX"), ]
  }
  if(require.SR==T){
    dat <- dat[grep("SR", dat$EVIDENCE), ]
  }
  
  print("processing all SVs ...")
  all.counts <- getSVdat.reannotated.permutate(dat, re_anno_svtype, genes, prefix="all")
  print("processing common SVs ...")
  common.counts <- getSVdat.reannotated.permutate(dat[which(dat$AF>=0.01), ], re_anno_svtype, genes, prefix="common")
  print("processing rare SVs ...")
  rare.counts <- getSVdat.reannotated.permutate(dat[which(dat$AF<0.01), ], re_anno_svtype, genes, prefix="rare")
  print("processing ultra_rare SVs ...")
  ultra_rare.counts <- getSVdat.reannotated.permutate(dat[which(dat$AF<0.001), ], re_anno_svtype, genes, prefix="under01perc")
  print("processing singleton SVs ...")
  singleton.counts <- getSVdat.reannotated.permutate(dat[which(dat$AC==1), ], re_anno_svtype, genes, prefix="singleton")
  print("processing large rare SVs ...")
  rare.large.counts <- getSVdat.reannotated.permutate(dat[which(dat$AF<0.01 & dat$SVLEN>100000), ], re_anno_svtype, genes, prefix="rare_large")
  #Merge data
  merged <- merge(x=snv.data, y=all.counts,  by="gene", all.x=T, sort=F)
  merged <- merge(x=merged, y=common.counts, by="gene", all.x=T, sort=F)
  merged <- merge(x=merged, y=rare.counts, by="gene", all.x=T, sort=F)
  merged <- merge(x=merged, y=ultra_rare.counts, by="gene", all.x=T, sort=F)
  merged <- merge(x=merged, y=singleton.counts, by="gene", all.x=T, sort=F)
  merged <- merge(x=merged, y=rare.large.counts, by="gene", all.x=T, sort=F)
  
  return(merged)
}

merge_loeuf_phaplo_ptriplo<-function(snv.data, phaplo_ptriplo){
  snv.data.2 = merge(snv.data, phaplo_ptriplo, by='gene')
  snv.data.2 = snv.data.2[!is.na(snv.data.2$pTriplo),]
  snv.data.2 = snv.data.2[order(snv.data.2$LOEUF),]
  snv.data.2$LOEUF_tile =as.integer(c(1:nrow(snv.data.2) - 1)/nrow(snv.data.2)*10)+1
  snv.data.2$LOEUF_percentile =as.integer(c(1:nrow(snv.data.2) - 1)/nrow(snv.data.2)*100)+1
  snv.data.2 = snv.data.2[order(snv.data.2$pHaplo),]
  snv.data.2$pHaplo_tile = 10 - as.integer(c(1:nrow(snv.data.2) - 1)/nrow(snv.data.2)*10)
  snv.data.2$pHaplo_percentile = 100 - as.integer(c(1:nrow(snv.data.2) - 1)/nrow(snv.data.2)*100)
  snv.data.2 = snv.data.2[order(snv.data.2$pTriplo),]
  snv.data.2$pTriplo_tile = 10 - as.integer(c(1:nrow(snv.data.2) - 1)/nrow(snv.data.2)*10)
  snv.data.2$pTriplo_percentile =100 - as.integer(c(1:nrow(snv.data.2) - 1)/nrow(snv.data.2)*100)
  return(snv.data.2)
}

# ---- Main ----

library("optparse")

option_list <- list(
  make_option(c("-p", "--permutate"),    type = "character", default = NULL,
              help = "Permutation round label, e.g. permi_1", metavar = "character"),
  make_option(c("-s", "--sv"),      type = "character", default = NULL,
              help = "SV info table (tsv/tsv.gz)", metavar = "character"),
  make_option(c("-g", "--gene"),    type = "character", default = NULL,
              help = "Gene info table (tsv/tsv.gz)", metavar = "character"),
  make_option(c("-r", "--reanno"), type = "character", default = NULL,
              help = "Re-annotated SV-vs-gene file (e.g. gnomAD_SV_v3.vs.hg38_gencode_v39.permuted_seed1.integrated)", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output .RData path", metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

required_args <- c("sv", "gene", "reanno", "output", "permutate")
missing_args  <- required_args[sapply(required_args, function(x) is.null(opt[[x]]) || opt[[x]] == "")]
if (length(missing_args) > 0) {
  stop(paste("Missing required arguments:", paste(paste0("-", substr(missing_args, 1, 1)), collapse = ", ")))
}

permu <- opt$permutate

print('reading in sv info ...')
sv_info <- read.table(opt$sv, header = TRUE, comment.char = "", sep = '\t')
sm_depth_only_dup <- sv_info[sv_info$SVTYPE == "DUP" & sv_info$ALGORITHMS == "depth" & sv_info$SVLEN < 20000, ]
sv_info <- sv_info[!sv_info$name %in% sm_depth_only_dup$name, ]
sv_info <- sv_info[sv_info$SVLEN < 1000000, ]

print('reading in gene info ...')
snv.data.permu <- read.table(opt$gene, header = TRUE, sep = '\t', comment.char = "")
snv.data.permu$gene <- paste(snv.data.permu$gene_name,  permu, sep = '.')

print('reading in sv vs. genes ...')
re_anno_svtype.permutate <- readin.re_anno_svtype(opt$reanno, sv_info)
re_anno_svtype.permutate <- re_anno_svtype.permutate[!re_anno_svtype.permutate$name %in% sm_depth_only_dup$name, ]
re_anno_svtype.permutate <- merge(re_anno_svtype.permutate, sv_info[, c("name", "SVLEN", "AF", "AC")])

print('organize sv vs. genes ...')
gene.data.reanno.permu <- getSVdat.all.reannotated.permutate(
  dat            = sv_info,
  re_anno_svtype = re_anno_svtype.permutate,
  snv.data       = snv.data.permu
)

save(gene.data.reanno.permu, file = opt$output)

print('writing tsv.gz output ...')
tsv_gz_output <- sub("(\\.rData|\\.RData|\\.rdata)$", ".tsv.gz", opt$output)
if (tsv_gz_output == opt$output) tsv_gz_output <- paste0(opt$output, ".tsv.gz")
con <- gzcon(file(tsv_gz_output, open = "wb"))
write.table(gene.data.reanno.permu, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
close(con)
cat(sprintf("Written %d rows x %d cols -> %s\n", nrow(gene.data.reanno.permu), ncol(gene.data.reanno.permu), tsv_gz_output))

