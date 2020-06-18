#! /usr/bin/env Rscript

library("caret")

args <- commandArgs(TRUE)

##inputs##
##Note file header should have with SVID,Status,Metric1...MetricX###
##Metric name is unique to the test and can describe the test statistic##
##Status in train file should indicate whether variant is passing or failing##

trainfile=args[1]
allfile=args[2]
seed=args[3]
output=args[4]

##read file and assign SVIDs as row names since they are not a metric for SV trainging##
metricsTrain<-read.table(trainfile,header=TRUE,row.names=1)
metricsAll<-read.table(allfile,header=TRUE,row.names=1)

##Parameter Tune##
fitControl <- trainControl(method = "oob")

##Use seed to ensure reproducibility## 
set.seed(seed)

##Run Model##       
rfresult<- train(as.factor(Status) ~ ., data = metricsTrain, 
                 method = 'rf', 
                 trControl = fitControl, metric = "Accuracy")

##Probablity Precdiction across all samples##
finalpred<-predict(rfresult,metricsAll,type='prob')


##output a probability of being valid for each variant##
output_matrix<-cbind(row.names(metricsAll), finalpred[,2])
colnames(output_matrix)<-c("name","prob")
write.table(output_matrix, paste(output,".pred",sep=""), 
            sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)

##save model for future use if necessary##
saveRDS(rfresult, paste(output,".rf",sep=""))
