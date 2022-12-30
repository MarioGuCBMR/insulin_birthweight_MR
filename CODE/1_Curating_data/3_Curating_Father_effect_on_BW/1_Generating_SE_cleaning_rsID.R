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

father <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Icelandic_BW_data/Birthweight_offspring_fathers2021.gz")

###############################
#REMOVING INSERTIONS/DELETIONS#
###############################

yes_vect <- c("A", "C", "G", "T")

father_clean <- father[which(father$A0%in%yes_vect),]
father_clean <- father_clean[which(father_clean$A1%in%yes_vect),]

###############
#Generating SE#
###############

z <- (as.numeric(sign(father_clean$`Beta-A1`)) * abs( qnorm(as.numeric(father_clean$P)/2)))

se_ <- (as.numeric(father_clean$`Beta-A1`)/as.numeric(z)) #this works really well.

father_clean$SE <- se_

#Small test with 1st SNP in the supplementary table:

father_clean[which(father_clean$rsID == "rs1374204"),]

#Chr      Pos      rsID A0 A1 IS-frq  IS-info Beta-A1          P          SE
#1: chr2 46257066 rs1374204  T  C 27.917 0.998347  -0.062 9.9013e-16 0.007722896

#PERFECT: it worked like a charm.

rm(father)
father <- father_clean
rm(father_clean)

#####################################
#QC ON DATA TO MAKE THINGS GO FASTER#
#####################################

#liftover dies every time we want to change from build 38 to 37.
#I am going to remove some of the variants before.

#INFO < 0.7
#MAF < 0.01.

summary(as.numeric(father$`IS-frq`)) #they are in percentage...
summary(as.numeric(father$`IS-info`)) #INFO is all good, incredibly.

#That means that we are gonna remove those with MAF > 0.01 in the Icelandic data that has data on AF on all
#and then check on EGG data whether they have it or not.

father <- father[which(as.numeric(father$`IS-frq`) > 1),] #we go from 23M to 13M 
father <- father[which(as.numeric(father$`IS-frq`) < 99),] #we stay in the 13;

#We should be careful and read the paper before removing more
#But it is important to do some QC on the heterogeneity and the AF of the EGG that do not match with the icelandic.
#Let's do a small test here:

father_egg_icelandic <- father[which(father$I2 != "na"),]

#Let's worry about only those that are heterogenic:

father_egg_icelandic <- father_egg_icelandic[which(father_egg_icelandic$I2 > 25),] #2M. Not too bad

#plot(father_egg_icelandic$`EGG-frq`, father_egg_icelandic$`IS-frq`) #the plot takes a while to run, but it runs pretty well. It does what I want, but I am not entirely convinced.

#I think this should be enough for the time being, we can remove the rest in the follow-up cleaning.
#We will remove the MHC region in the build 37 using the 26M to 34M afterwords.

######################
#Taking care of rsIDs#
######################

#I know that some rsIDs have some special characters.
#We are going to investigate those.

father <- father[order(father$rsID),]

head(father$rsID) #these are normal.
tail(father$rsID) #we have NAs, but that is to be expected

#Let's do some checks:

father_rsid_na <- father[which(is.na(father$rsID) == TRUE),] #2M approximately.
father_rsid_merged <- father[which(str_detect(father$rsID, ",") == TRUE),]
father_rsid_normal <- father[which(!(father$rsID%in%father_rsid_merged$rsID | father$rsID%in%father_rsid_na$rsID)),]

father_rsid_normal <- father_rsid_normal[order(father_rsid_normal$rsID),]

head(father_rsid_normal$rsID) #these are normal.
tail(father_rsid_normal$rsID) #awesome!! This data seems cleaned enough.

#Let's clean the rsID in the merged variants section, just in case.
#It seems that the first rsID is the newest, and then they add the merged one.
#Not that I 100% care because we are going to use the data the chr_pos, but in any case.

father_rsid_merged$merged_rsID <- as.character(unlist(sapply(father_rsid_merged$rsID, merged_rsid_parser)))
father_rsid_merged$rsID <- as.character(unlist(sapply(father_rsid_merged$rsID, new_rsid_parser)))

#Now let's get the rest of the dataframes to have the same format:

father_rsid_normal$merged_rsID <- "-"
father_rsid_na$merged_rsID <- "-"

father_rsid_end <- rbind(father_rsid_merged, father_rsid_na, father_rsid_normal)

##############################################################
#Finally, let's just take care of the chromosome and position#
##############################################################

#I am going to save it as RAW since I have not removed MHC, or MAF<0.01 and such:

father_rsid_end$chr_pos_38 <- paste(father_rsid_end$Chr, father_rsid_end$Pos, sep = ":") 

#The build is 38!!! Don´t be mislead by the variants that have the same 
#position in chromosome 1.
  
fwrite(father_rsid_end, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Icelandic_BW_data/Birthweight_offspring_fathers2021_with_SE.txt")


