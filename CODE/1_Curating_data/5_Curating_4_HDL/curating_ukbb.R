##############
#INTRODUCTION#
##############

#This is a check: let's go!!

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#LOADING FUNCTIONS#
###################

parse_chr_pos <- function(id){
  
  chr_ <- strsplit(id, ":")[[1]][1] #this works for those that do not present the weird formatting :)
  pos_ <- strsplit(id, ":")[[1]][2] #this works for those that do not present the weird formatting :)
  
  chr_pos <- paste("chr", chr_, ":", pos_, sep = "")
  
  return(chr_pos)
  
}

##################
#LOADING BMI DATA#
##################

load("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/UKBB_ref_panel/snp.dictionary.imputed.rda")

snp.dictionary$chr_pos <- as.character(unlist(sapply(snp.dictionary$variant, parse_chr_pos)))

######################
#SAVE THE CURATED BMI#
######################

fwrite(snp.dictionary, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/UKBB_ref_panel/ukbb_ref_panel_curated.txt", sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)
