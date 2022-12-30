##############
#INTRODUCTION#
##############

#This is code to clean FBW data for GC (HDL).

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

FBW <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/OUTPUT/FBW/Birthweight2021_Curated.txt")

######################################################
#PERFECT: LET'S GO WITH SELECTING THE COLUMNS FOR HDL#
######################################################

FBW$eaf <- as.numeric(FBW$`EGG-frq`)/100

FBW$maf <- ifelse(FBW$eaf > 0.5, 1-FBW$eaf, FBW$eaf)

FBW_4_ldsc <- FBW %>%
  select("rsID", "Beta-A1", "A1", "A0", "maf", "P", "IS-info")

colnames(FBW_4_ldsc) <- c("SNP", "BETA", "A1", "A2", "MAF", "P", "INFO")

fwrite(FBW_4_ldsc, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/DATA_4_MUNGING/FBW_4_LDSC.txt", sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "-")

#########################
#Let's do some checks...#
#########################

check <- FBW_4_ldsc[which(FBW_4_ldsc$SNP != ""),]
check <- check[which(is.na(check$MAF) == FALSE),]
check <- check[which(is.na(check$BETA) == FALSE),]

fwrite(check, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/DATA_4_MUNGING/FBW_4_LDSC.txt", sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "-")

##################################
#LET'S MUNGE THE HELL OUT OF THIS#
##################################

#cd J:\CBMR\SUN-CBMR-Kilpelainen-Group\Team projects\Hermina_and_Mario\IR_BW_AIM4\CODE\2_Genetic_Correlations\LDSC\ldsc-master\ldsc-master

python munge_sumstats.py --sumstats "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/DATA_4_MUNGING/FBW_4_LDSC.txt" --N 423683 --out "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/MUNGED_DATA/FBW_4_LDSC.txt" --merge-alleles w_hm3.snplist
