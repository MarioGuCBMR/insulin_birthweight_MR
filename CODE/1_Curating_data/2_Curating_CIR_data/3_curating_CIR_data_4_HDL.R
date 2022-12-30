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

CIR <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/CIR/CIR_combined_Curated.txt")

ukbb <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/UKBB_ref_panel/ukbb_ref_panel_curated.txt")

length(which(CIR$snp%in%ukbb$rsid)) #not enough for this reference panel, but enough for the other one.

######################################################
#PERFECT: LET'S GO WITH SELECTING THE COLUMNS FOR HDL#
######################################################

CIR$N <- 5318

CIR_4_hdl <- CIR %>%
  select(snp, effect_allele, other_allele, N, effect, stderr)

colnames(CIR_4_hdl) <- c("SNP", "A1", "A2", "N", "b", "se")

fwrite(CIR_4_hdl, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/MUNGED_DATA/CIR_4_HDL.txt", sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)
