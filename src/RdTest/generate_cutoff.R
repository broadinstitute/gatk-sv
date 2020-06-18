#! /usr/bin/env Rscript

##meant to generate cutoffs for different copy states##
##first input is a Rd-test generated median_geno file from genotyping command on a set of know sites##
## second input is how many copy states to assess##
## third input is output filename for final matrix of cutoffs##

args <- commandArgs(TRUE)

file=args[1]
max_copystate=as.numeric(args[2])
outputfile=args[3]

median_matrix=read.table(file,header=TRUE, na.string='.')

##remove any sites without bincov coverage##
median_matrix=median_matrix[!is.na(median_matrix[,5]),]

##use expected cutoffs for genotypeing (0.25,.75,1.25,1.75)
cutoffs <- seq(0.25, max_copystate*0.5+0.25 , by = .5)

cnv_median=unlist(median_matrix[,5:ncol(median_matrix)])

prev_cutoff=0
copystate<-cnv_median
output_matrix=matrix(,ncol=3,nrow=max_copystate+1)
for (i in 0:max_copystate) {
  copystate.mean<-mean(cnv_median[which(cnv_median <= cutoffs[i+1] & cnv_median > prev_cutoff) ])
  copystate.sd<-sd(cnv_median[which(cnv_median <= cutoffs[i+1] & cnv_median > prev_cutoff) ])
  ##Assign 0.5*copystate  if nothing to train on for given copy state##
  if (is.na(copystate.mean)==TRUE){
    copystate.mean=(i+1)*0.5
  }
  ##Assign previous sd if nothing to train on for given copy state##
  if (is.na(copystate.sd)==TRUE){
    copystate.sd=output_matrix[i,3]
  }  
  output_matrix[i+1,]<-c(i,copystate.mean,copystate.sd)
  prev_cutoff=cutoffs[i+1]
} 

zcalc <- function(mean1,sigma1,mean2,sigma2)  { return((mean2/sigma2+mean1/sigma1)/(1/sigma1+1/sigma2))}

transistion_cutoffs<-c()
for (i in 0:max_copystate-1) {
   mean1=output_matrix[i+1,2]
   sigma1=output_matrix[i+1,3]
   mean2=output_matrix[i+2,2]
   sigma2=output_matrix[i+2,3]
   transistion_cutoffs[i+1]<-zcalc(mean1,sigma1,mean2,sigma2)
}

transistion_cutoffs<-c(transistion_cutoffs,max_copystate*0.5+0.25)

final_matrix<-cbind(output_matrix,transistion_cutoffs)

colnames(final_matrix)<-c("copy_state","mean","sd","cutoffs")

write.table(final_matrix,outputfile,row.names = FALSE,quote = FALSE,sep= "\t")
 