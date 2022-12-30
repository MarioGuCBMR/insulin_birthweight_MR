##############
#INTRODUCTION#
##############

#This is code to clean MBW data for GC (HDL).

#We cannot map the data perfectly, because we are missing some SNPs.
#The reason behind this is that many SNPs are missing the rsID.
#But I have a way of solving that ;) Using the imputed reference panel and loading them.

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


se_calculator <- function(beta_, pval_){
  
  z = (as.numeric(sign(beta_)) * abs( qnorm(as.numeric(pval_)/2)))
  
  se_ <- (as.numeric(beta_)/as.numeric(z)) #this works really well.
  
  return(se_)
  
}

new_rsid_parser <- function(rsid){
  #This allele takes the rs,rs format and retrieves the first rs.
  
  rsid <-strsplit(rsid, ",")[[1]][1]
  
  return(rsid)
  
}

merged_rsid_parser <- function(rsid){
  #This allele takes the rs,rs format and retrieves the first rs.
  
  rsid <-strsplit(rsid, ",")[[1]][2]
  
  return(rsid)
  
}


##############
#Loading data#
##############

MBW <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/OUTPUT/MBW/Birthweight_offspring_mothers2021_curated_FULL.txt")

ukbb <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/UKBB_ref_panel/ukbb_ref_panel_curated.txt")

##########################
#Now, starts the fun baby#
##########################

MBW_rsid <- MBW[which(str_detect(MBW$rsID, "rs") == TRUE),]
MBW_rsid_fail <- MBW[which(str_detect(MBW$rsID, "rs") == FALSE),]

#Now let's do some checks...

length(which(MBW_rsid$rsID%in%ukbb$rsid)) #(1028666)/(1029876) Most of them!! but we are missing some fellas.
length(which(MBW_rsid$merged_rsID%in%ukbb$rsid)) #enough to be OK.

#Nothing we can do...
#We will just have to work with this.

######################################################
#PERFECT: LET'S GO WITH SELECTING THE COLUMNS FOR HDL#
######################################################

MBW$N <- 270002

MBW_4_hdl <- MBW %>%
  select("rsID", "A1", "A0", "N", "Beta-A1", "SE")

colnames(MBW_4_hdl) <- c("SNP", "A1", "A2", "N", "b", "se")

MBW_4_hdl$SNP <- ifelse(str_detect(MBW_4_hdl$SNP, "rs") == FALSE, "-", MBW_4_hdl$SNP) #usually we expect them to have mean 0 and error of 1 so, we are gonna set it up like that.
MBW_4_hdl$se <- ifelse(is.na(MBW_4_hdl$se) == TRUE, 1, MBW_4_hdl$se) #usually we expect them to have mean 0 and error of 1 so, we are gonna set it up like that.

fwrite(MBW_4_hdl, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/MUNGED_DATA/MBW_4_HDL.txt", sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)
