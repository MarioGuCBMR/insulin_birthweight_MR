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

parse_rsid <- function(rsid){
  
  rsid_ <- strsplit(rsid, ":")[[1]][1] #this works for those that do not present the weird formatting :)
  
  return(rsid_)
  
}

##################
#LOADING BMI DATA#
##################

bmi <- fread("C:/Users/zlc436/Desktop/PULIT_data_reservoir/Combined/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")

#We need to be careful, some SNPs from Pulit et al do not present correct chr_pos.
#The mistake is in those without INFO data, thus, the ones in GIANT data.
#We are gonna remove those.

summary(bmi$INFO) 
summary(bmi$Freq_Tested_Allele) 
summary(bmi$CHR) 
summary(bmi$N) 

bmi <- bmi[which(is.na(bmi$INFO) == FALSE),]

summary(bmi$INFO) #all good
summary(bmi$Freq_Tested_Allele) #all good
summary(bmi$CHR) #all good
summary(bmi$N) #all good

#It seems all is OK.
#Now we need to clean the RSIDs, because the have the effect allele and the other allele.
#Which is actually fine.

bmi$RSID <- as.character(unlist(sapply(bmi$SNP, parse_rsid)))

bmi$chr_pos <- paste("chr", bmi$CHR, ":", bmi$POS, sep = "")

######################
#SAVE THE CURATED BMI#
######################

fwrite(bmi, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/BMI/bmi_curated.txt", sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)
