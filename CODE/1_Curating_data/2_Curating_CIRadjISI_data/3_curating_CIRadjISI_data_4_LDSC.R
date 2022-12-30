##############
#INTRODUCTION#
##############

#This is code to clean FI data for GC (HDL).

#I checked beforehand, and there is a perfect overlap with the HapMap2 imputed data.
#There is no need to do some checkity checks.

#################
#Loading libries#
#################

library(data.table)
library(cause)
library(tidyverse)

#########################
#Making space for memory#
#########################

memory.limit(size=80000000)

###################
#Loading functions#
###################

chr_parser <- function(chr_pos){
  #Simple function: it introduces a SNP in chr1:1000 format and returns the chromosomes as a number.
  
  options(scipen=999)
  
  chr_ <- strsplit(chr_pos, ":")[[1]][1]
  chr_end <- as.numeric(as.character(strsplit(chr_, "chr")[[1]][2]))
  
  return(chr_end)
  
}

pos_parser <- function(chr_pos){
  #Simple function: it introduces a SNP in chr1:1000 format and returns the positions as a number.
  
  options(scipen=999) #to avoid any weird transformations.
  
  pos_ <- as.numeric(as.character(strsplit(chr_pos, ":")[[1]][2]))
  
  return(pos_)
  
}

##############
#Loading data#
##############

CIRadjISI <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/CIRadjISI/CIRadjISI_combined_Curated.txt")

##############################
#CLEANING MHC REGION FOR LDSC#
##############################

CIRadjISI$chr <- as.numeric(as.character(unlist(sapply(CIRadjISI$chr_pos_37, chr_parser))))
CIRadjISI$pos <- as.numeric(as.character(unlist(sapply(CIRadjISI$chr_pos_37, pos_parser))))

CIRadjISI_mhc <- CIRadjISI[which(CIRadjISI$chr == 6 & CIRadjISI$pos > 26000000 & CIRadjISI$pos < 34000000),]

summary(CIRadjISI_mhc$chr) #perfect
summary(CIRadjISI_mhc$pos) #perfect

CIRadjISI <- CIRadjISI[which(!(CIRadjISI$chr_pos_37%in%CIRadjISI_mhc$chr_pos_37)),]

######################################################
#PERFECT: LET'S GO WITH SELECTING THE COLUMNS FOR HDL#
######################################################

CIRadjISI_4_ldsc <- CIRadjISI%>%
  select(snp, effect, effect_allele, other_allele, maf, pvalue)

colnames(CIRadjISI_4_ldsc) <- c("SNP", "BETA", "A1", "A2", "MAF", "P")

fwrite(CIRadjISI_4_ldsc, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/DATA_4_MUNGING/CIRadjISI_4_LDSC.txt", sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

python munge_sumstats.py --sumstats "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/DATA_4_MUNGING/CIRadjISI_4_LDSC.txt" --N 4789 --out "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/MUNGED_DATA/CIRadjISI_4_LDSC.txt" --merge-alleles w_hm3.snplist
