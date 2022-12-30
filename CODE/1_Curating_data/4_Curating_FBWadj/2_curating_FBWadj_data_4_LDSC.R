##############
#INTRODUCTION#
##############

#This is code to clean fbwadjmbw data for GC (HDL).

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

fbwadjmbw <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/FBWadj/FBWAdj_curated.txt")

######################################################
#PERFECT: LET'S GO WITH SELECTING THE COLUMNS FOR HDL#
######################################################

fbwadjmbw$maf <- ifelse(fbwadjmbw$eaf > 0.5, 1-fbwadjmbw$eaf, fbwadjmbw$eaf)

fbwadjmbw_4_ldsc <- fbwadjmbw %>%
  select("RSID", "beta", "ea", "nea", "maf", "p", "n_ownBW")

colnames(fbwadjmbw_4_ldsc) <- c("SNP", "BETA", "A1", "A2", "MAF", "P", "N")

fwrite(fbwadjmbw_4_ldsc, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/DATA_4_MUNGING/FBWadj_4_LDSC.txt", sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "-")

##################################
#LET'S MUNGE THE HELL OUT OF THIS#
##################################

#cd J:\CBMR\SUN-CBMR-Kilpelainen-Group\Team projects\Hermina_and_Mario\IR_BW_AIM4\CODE\2_Genetic_Correlations\LDSC\ldsc-master\ldsc-master

python munge_sumstats.py --sumstats "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/DATA_4_MUNGING/FBWadj_4_LDSC.txt" --out "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/MUNGED_DATA/FBWadj_4_LDSC.txt" --merge-alleles w_hm3.snplist
