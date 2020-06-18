#!/usr/bin/env Rscript

#---------------------------
# R package check for SV-pipeline
# Talkowski Laboratory
#
# Harrison Brand and Joseph Glessner
# Update June 2017
#---------------------------

#Loads required packages; installs if necessary
RPackages <- c("caret","optparse", "plyr", "MASS","ROCR", "zoo","data.table","methods","metap", "e1071", "fpc", "BSDA", "DAAG", "pwr", "reshape", "perm", "hash")
for (i in RPackages)
{
  if (i %in% rownames(installed.packages()) == FALSE) {
      print(paste("installing package:",i,sep=""))
      install.packages((i), repos = "http://cran.rstudio.com")
      library(i, character.only = TRUE)
  } else {
    library(i, character.only = TRUE)
  }
}
