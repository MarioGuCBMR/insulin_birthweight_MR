##############
#INTRODUCTION#
##############

#This code is to curate the data from the Icelandic+EGG+UKBB.

#The main issue is that we do not have the SE, so let's generate it!

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)

memory.limit(size=8000000)

###################
#Loading functions#
###################

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

mbw <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Icelandic_BW_data/Birthweight_offspring_mothers2021.gz")

###############################
#REMOVING INSERTIONS/DELETIONS#
###############################

yes_vect <- c("A", "C", "G", "T")

mbw_clean <- mbw[which(mbw$A0%in%yes_vect),]
mbw_clean <- mbw_clean[which(mbw_clean$A1%in%yes_vect),]

###############
#Generating SE#
###############

z <- (as.numeric(sign(mbw_clean$`Beta-A1`)) * abs( qnorm(as.numeric(mbw_clean$P)/2)))

se_ <- (as.numeric(mbw_clean$`Beta-A1`)/as.numeric(z)) #this works really well.

mbw_clean$SE <- se_

head(mbw_clean) #PERFECT

#Small test with 1st SNP in the supplementary table:

mbw_clean[which(mbw_clean$rsID == "rs17367504"),]

#Chr      Pos       rsID A0 A1 IS-frq IS-info EGG-frq Beta-A1
#1: chr1 11802721 rs17367504  A  G   15.5       1   16.67   0.032
#P     I2 P-het          SE
#1: 1.9e-17 63.365 0.099 0.003764817

#PERFECT.

rm(mbw)
mbw <- mbw_clean
rm(mbw_clean)

#####################################
#QC ON DATA TO MAKE THINGS GO FASTER#
#####################################

#liftover dies every time we want to change from build 38 to 37.
#I am going to remove some of the variants before.

#INFO < 0.7
#MAF < 0.01.

summary(as.numeric(mbw$`IS-frq`)) #they are in percentage...
summary(as.numeric(mbw$`EGG-frq`)) #16049576 NAs. Mother of God.
summary(as.numeric(mbw$`IS-info`)) #INFO is all good, incredibly.

#That means that we are gonna remove those with MAF > 0.01 in the Icelandic data that has data on AF on all
#and then check on EGG data whether they have it or not.

mbw_NA_egg <- mbw[which(is.na(as.numeric(mbw$`EGG-frq`)) == TRUE),]
mbw_ok_egg <- mbw[which(is.na(as.numeric(mbw$`EGG-frq`)) == FALSE),]

summary(as.numeric(mbw_ok_egg$`IS-frq`)) #PERF
summary(as.numeric(mbw_ok_egg$`EGG-frq`)) #PERF

mbw_ok_egg <- mbw_ok_egg[which(as.numeric(mbw_ok_egg$`EGG-frq`) > 1),] #we go from 33M to 13M 
mbw_ok_egg <- mbw_ok_egg[which(as.numeric(mbw_ok_egg$`EGG-frq`) < 99),] #we stay in the 13;

summary(as.numeric(mbw_ok_egg$`IS-frq`)) #PERF
summary(as.numeric(mbw_ok_egg$`EGG-frq`)) #PERF. There are some mismatches with IS data, but we trust EGG data, not IS.

mbw_NA_egg <- mbw_NA_egg[which(as.numeric(mbw_NA_egg$`IS-frq`) > 1),] #we go from 33M to 13M 
mbw_NA_egg <- mbw_NA_egg[which(as.numeric(mbw_NA_egg$`IS-frq`) < 99),] #we stay in the 13;

summary(as.numeric(mbw_NA_egg$`IS-frq`)) #PERF
summary(as.numeric(mbw_NA_egg$`EGG-frq`)) #PERF. There are some mismatches with IS data, but we trust EGG data, not IS.

mbw <- rbind(mbw_NA_egg, mbw_ok_egg)

rm(mbw_NA_egg)
rm(mbw_ok_egg)

######################
#Taking care of rsIDs#
######################

#I know that some rsIDs have some special characters.
#We are going to investigate those.

mbw <- mbw[order(mbw$rsID),]

head(mbw$rsID) #these are normal.
tail(mbw$rsID) #we have NAs, but that is to be expected

#Let's do some checks:

mbw_rsid_na <- mbw[which(is.na(mbw$rsID) == TRUE),] #2M approximately.
mbw_rsid_merged <- mbw[which(str_detect(mbw$rsID, ",") == TRUE),]
mbw_rsid_normal <- mbw[which(!(mbw$rsID%in%mbw_rsid_merged$rsID | mbw$rsID%in%mbw_rsid_na$rsID)),]

mbw_rsid_normal <- mbw_rsid_normal[order(mbw_rsid_normal$rsID),]

head(mbw_rsid_normal$rsID) #these are normal.
tail(mbw_rsid_normal$rsID) #awesome!! This data seems cleaned enough.

#Let's clean the rsID in the merged variants section, just in case.
#It seems that the first rsID is the newest, and then they add the merged one.
#Not that I 100% care because we are going to use the data the chr_pos, but in any case.

mbw_rsid_merged$merged_rsID <- as.character(unlist(sapply(mbw_rsid_merged$rsID, merged_rsid_parser)))
mbw_rsid_merged$rsID <- as.character(unlist(sapply(mbw_rsid_merged$rsID, new_rsid_parser)))

#Now let's get the rest of the dataframes to have the same format:

mbw_rsid_normal$merged_rsID <- "-"
mbw_rsid_na$merged_rsID <- "-"

mbw_rsid_end <- rbind(mbw_rsid_merged, mbw_rsid_na, mbw_rsid_normal)

##############################################################
#Finally, let's just take care of the chromosome and position#
##############################################################

#I am going to save it as RAW since I have not removed MHC, or MAF<0.01 and such:

mbw_rsid_end$chr_pos_38 <- paste(mbw_rsid_end$Chr, mbw_rsid_end$Pos, sep = ":") 

#The build is 38!!! Don´t be mislead by the variants that have the same 
#position in chromosome 1.
  
fwrite(mbw_rsid_end, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Icelandic_BW_data/Birthweight_offspring_mothers2021_with_SE.txt")
