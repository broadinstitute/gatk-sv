#! /usr/bin/env Rscript

##meant to generate cutoffs for different copy states##
##first input is a PE file##
## second input is output filename for final matrix of cutoffs##

args <- commandArgs(TRUE)

file=args[1]
outputfile=args[2]

median_matrix=read.table(file)

output_matrix=matrix(,ncol=3,nrow=5)
for (i in 0:3) {
  copystate.mean<-mean(median_matrix[which(median_matrix[,3] == i ) ,2])
  copystate.sd<-sd(median_matrix[which(median_matrix[,3] == i ) ,2])    
  output_matrix[i+1,]<-c(i,copystate.mean,copystate.sd)
  } 

zcalc <- function(mean1,sigma1,mean2,sigma2)  { return((mean2/sigma2+mean1/sigma1)/(1/sigma1+1/sigma2))}


transistion_cutoffs<-c()

##pick between 1 and 3 which has lower sd##
copy_state0<-output_matrix[3,]
copy_state1=output_matrix[which(output_matrix[,3]==min(output_matrix[c(2,4),3])),]
copy_state1.mean=copy_state1[2]
copy_state1.sd=copy_state1[3]

transistion_cutoffs.copystate1<-c("transition1",zcalc(copy_state1.mean,copy_state1.sd,output_matrix[3,2],output_matrix[3,3]))

copy_state2=output_matrix[which(output_matrix[,3]==min(output_matrix[c(1),3])),]
copy_state2.mean=copy_state2[2]
copy_state2.sd=copy_state2[3]

transistion_cutoffs.copystate2<-c("transition2",zcalc(copy_state2.mean,copy_state2.sd,copy_state1.mean,copy_state1.sd))

metric_matrix<-rbind(copy_state0,copy_state1,copy_state2)

colnames(metric_matrix)<-c("copy_state","mean","sd")

final_matrix<-rbind(transistion_cutoffs.copystate1,transistion_cutoffs.copystate2)

colnames(final_matrix)<-c("trainsition_state","cutoffs")

write.table(final_matrix,outputfile,row.names = FALSE,quote = FALSE,sep= "\t")

write.table(metric_matrix,paste(outputfile,".metric",sep=""),row.names = FALSE,quote = FALSE,sep= "\t") 
