##############
#INTRODUCTION#
##############

#This code is to curate the data from the Icelandic+EGG+UKBB.
#The main issue is that we do not have the SE, so let's generate it!
#We are also are going to remove the insertions and deletions today.

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

fbw <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Icelandic_BW_data/Birthweight2021.gz")

###############################
#REMOVING INSERTIONS/DELETIONS#
###############################

yes_vect <- c("A", "C", "G", "T")

fbw_clean <- fbw[which(fbw$A0%in%yes_vect),]
fbw_clean <- fbw_clean[which(fbw_clean$A1%in%yes_vect),]

###############
#Generating SE#
###############

z <- (as.numeric(sign(fbw_clean$`Beta-A1`)) * abs( qnorm(as.numeric(fbw_clean$P)/2)))

se_ <- (as.numeric(fbw_clean$`Beta-A1`)/as.numeric(z)) #this works really well.

fbw_clean$SE <- se_

head(fbw_clean) #PERFECT

#Small test with 1st SNP in the supplementary table:

fbw_clean[which(fbw_clean$rsID == "rs12401656"),]

#Chr      Pos       rsID A0 A1 IS-frq IS-info EGG-frq Beta-A1       P     I2 P-het          SE
#1: chr1 42991096 rs12401656  G  A   13.2       1   13.49  -0.028 2.4e-14 78.960 0.029 0.003671096

#PERFECT: it worked like a charm.

rm(fbw)
fbw <- fbw_clean
rm(fbw_clean)

#####################################
#QC ON DATA TO MAKE THINGS GO FASTER#
#####################################

#Let's check the AF.

summary(as.numeric(fbw$`IS-frq`)) #they are in percentage...
summary(as.numeric(fbw$`EGG-frq`)) #15463564 NAs. Mother of God.
summary(as.numeric(fbw$`IS-info`)) #INFO is all good, incredibly.

fbw_NA_egg <- fbw[which(is.na(as.numeric(fbw$`EGG-frq`)) == TRUE),]
fbw_ok_egg <- fbw[which(is.na(as.numeric(fbw$`EGG-frq`)) == FALSE),]

summary(as.numeric(fbw_ok_egg$`IS-frq`)) #PERF
summary(as.numeric(fbw_ok_egg$`EGG-frq`)) #PERF

fbw_ok_egg <- fbw_ok_egg[which(as.numeric(fbw_ok_egg$`EGG-frq`) > 1),] #we go from 33M to 13M 
fbw_ok_egg <- fbw_ok_egg[which(as.numeric(fbw_ok_egg$`EGG-frq`) < 99),] #we stay in the 13;

summary(as.numeric(fbw_ok_egg$`IS-frq`)) #PERF
summary(as.numeric(fbw_ok_egg$`EGG-frq`)) #PERF. There are some mismatches with IS data, but we trust EGG data, not IS.

fbw_NA_egg <- fbw_NA_egg[which(as.numeric(fbw_NA_egg$`IS-frq`) > 1),] #we go from 33M to 13M 
fbw_NA_egg <- fbw_NA_egg[which(as.numeric(fbw_NA_egg$`IS-frq`) < 99),] #we stay in the 13;

summary(as.numeric(fbw_NA_egg$`IS-frq`)) #PERF
summary(as.numeric(fbw_NA_egg$`EGG-frq`)) #PERF. There are some mismatches with IS data, but we trust EGG data, not IS.

fbw <- rbind(fbw_NA_egg, fbw_ok_egg)

rm(fbw_NA_egg)
rm(fbw_ok_egg)

######################
#Taking care of rsIDs#
######################

#I know that some rsIDs have some special characters.
#We are going to investigate those.

fbw <- fbw[order(fbw$rsID),]

head(fbw$rsID) #these are normal.
tail(fbw$rsID) #we have NAs, but that is to be expected

#Let's do some checks:

fbw_rsid_na <- fbw[which(is.na(fbw$rsID) == TRUE),] #2M approximately.
fbw_rsid_merged <- fbw[which(str_detect(fbw$rsID, ",") == TRUE),]
fbw_rsid_normal <- fbw[which(!(fbw$rsID%in%fbw_rsid_merged$rsID | fbw$rsID%in%fbw_rsid_na$rsID)),]

fbw_rsid_normal <- fbw_rsid_normal[order(fbw_rsid_normal$rsID),]

head(fbw_rsid_normal$rsID) #these are normal.
tail(fbw_rsid_normal$rsID) #awesome!! This data seems cleaned enough.

#Let's clean the rsID in the merged variants section, just in case.
#It seems that the first rsID is the newest, and then they add the merged one.
#Not that I 100% care because we are going to use the data the chr_pos, but in any case.

fbw_rsid_merged$merged_rsID <- as.character(unlist(sapply(fbw_rsid_merged$rsID, merged_rsid_parser)))
fbw_rsid_merged$rsID <- as.character(unlist(sapply(fbw_rsid_merged$rsID, new_rsid_parser)))

#Now let's get the rest of the dataframes to have the same format:

fbw_rsid_normal$merged_rsID <- "-"
fbw_rsid_na$merged_rsID <- "-"

fbw_rsid_end <- rbind(fbw_rsid_merged, fbw_rsid_na, fbw_rsid_normal)

##############################################################
#Finally, let's just take care of the chromosome and position#
##############################################################

#I am going to save it as RAW since I have not removed MHC, or MAF<0.01 and such:

fbw_rsid_end$chr_pos_38 <- paste(fbw_rsid_end$Chr, fbw_rsid_end$Pos, sep = ":") 

#The build is 38!!! Don´t be mislead by the variants that have the same 
#position in chromosome 1.
  
fwrite(fbw_rsid_end, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Icelandic_BW_data/Birthweight2021_with_SE.txt")


