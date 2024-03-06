# This script creates a ped file with all families >3 into trios
# to match desired input for de novo SV filtering

# Define args and load libraries
args = commandArgs(trailingOnly=TRUE)
require(plyr)
library(data.table)
library(tidyverse)

input_ped <- args[1]

# Read input ped
ped <- fread(input_ped,header = TRUE)

# Add a column with size of family
ped %>%
  group_by(FamID) %>%
  mutate(FamSize = n()) -> ped2

# Steph's code
distfam<-ped2%>%distinct(FamID,FamSize)
fam_breakdown<-distfam%>%group_by(FamSize)%>%tally()

# subset the families > 3
large_families<-ped2%>%filter(FamSize>=4)
large_family_ids <- large_families$FamID # FamIDs of families >3

large_families_prob<-large_families%>%filter(MotherID!=0|FatherID!=0) # probands from families > 3
large_families_prob$newFamID <- NA
families<-unique(large_families_prob$FamID)

# function that splits by probands, adds _# to the end of the FamID and then finds the parents and adds _# to them as well
split_fam<-function(probands,pedigree){
  trio <- data.frame(matrix(ncol = 8, nrow = 0))
  families<-unique(probands$FamID)
  for (i in 1:length(families)){
    # reassign family IDs for >=quads (fam1 => fam1_1, fam1_2)
    perfam<-subset(probands,FamID==families[i])
    for (x in 1:nrow(perfam)){
      probands$newFamID[probands$IndividualID==as.character(perfam[x,2])]<-paste0(probands$FamID[probands$IndividualID==as.character(perfam[x,2])],"_",x)
    }
  }
  # each proband row should now have a new fam_id
  # pedigree needs a new column
  pedigree$newFamID<-NA
  # now go back and pick up parents info and rename famID of parents
  for (y in 1:nrow(probands)){
    paternal<-subset(pedigree, IndividualID == as.character(probands[y, 3]))
    paternal$newFamID <- as.character(probands[y,8])
    maternal<-subset(pedigree, IndividualID == as.character(probands[y, 4]))
    maternal$newFamID<-as.character(probands[y,8])
    trio<-rbind(data.frame(trio), data.frame(probands[y,]),data.frame(paternal),data.frame(maternal))
  }
  return(trio)
}

output<-split_fam(large_families_prob, large_families)

# subset out newFamID column and rename FamID
new_output <- subset(output, select=c("newFamID", "IndividualID", "MotherID", "FatherID", "Gender", "Affected"))
colnames(new_output)[1] <- "FamID"

# subset to just families with all three in the new ped
new_output %>%
  group_by(FamID) %>%
  mutate(FamSize = n()) %>%
  filter(FamSize==3)-> trios

new_trios <- subset(trios, select=c("FamID", "IndividualID", "MotherID", "FatherID", "Gender", "Affected"))

# If mom and dad are not in original ped file as an individual, remove them so that they do not cause an error in the de novo script
new_trios$FatherID[!(new_trios$FatherID %in% ped$IndividualID)] <- 0
new_trios$MotherID[!(new_trios$MotherID %in% ped$IndividualID)] <- 0

# subset ped to only include families that are not > 3
ped_without_large_families <- subset(ped, !(FamID %in% large_family_ids))

# If mom and dad are not in original ped file as an individual, remove them so that they do not cause an error in the de novo script
ped_without_large_families$FatherID[!(ped_without_large_families$FatherID %in% ped$IndividualID)] <- 0
ped_without_large_families$MotherID[!(ped_without_large_families$MotherID %in% ped$IndividualID)] <- 0

# concat new trio pedigree and old pedigree without families > 3
ped_without_large_families_df <- data.frame(ped_without_large_families)
new_trios_df <- data.frame(new_trios)

final_ped<-rbind(new_trios_df, ped_without_large_families_df)

write.table(final_ped,file="cleaned_ped.txt",quote=F,row.names=F,col.names=T,sep="\t")
