##############
#INTRODUCTION#
##############

#This code is to curate the data from the Icelandic+EGG+UKBB.
#Here we are going to try to get the build38 to build37.
#We are going to use the already curated MBW data which is gonna make life so much easier!!!

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)

memory.limit(size=800000000)

#Let's install liftover for BioConductor:

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("liftOver")
library(liftOver)

################
#Load functions#
################

parse_rsid <- function(rsid){
  
  rsid_ <- strsplit(rsid, ":")[[1]][1] #this works for those that do not present the weird formatting :)
  
  return(rsid_)
  
}


liftover_generator <- function(chr, pos){
  
  pos_1 <- as.numeric(pos)-1
  pos_2 <- as.numeric(pos)+1
  
  chr_pos_liftover <- paste(chr, ":", pos_1, "-", pos_2, sep = "")
  
  return(chr_pos_liftover)
  
}

clean_liftover_data <- function(chr_pos_liftover){
  #Takes liftover data (chr:pos1-pos2) in a build and transforms it to chr:pos data
  
  tmp_pos <- strsplit(chr_pos_liftover, "-")[[1]][2]
  chr_tmp <- strsplit(chr_pos_liftover, ":")[[1]][1]
  
  tmp_pos_clean <- as.numeric(tmp_pos)-1
  
  chr_pos <- paste(chr_tmp, ":", tmp_pos_clean, sep = "")
  
  return(chr_pos)
  
  
}

##############
#Loading data#
##############

father <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Icelandic_BW_data/Birthweight_offspring_fathers2021_with_SE.txt")

mbw_curated <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/OUTPUT/MBW/Birthweight_offspring_mothers2021_curated_FULL.txt")

###################
#Matching MBW data#
###################

mbw_match <- mbw_curated[which(mbw_curated$chr_pos_38%in%father$chr_pos_38),]
mbw_match_no_dupl <- mbw_match[which(duplicated(mbw_match$chr_pos_38) == FALSE),]

#We have duplicates, but we do not care about that!
#We care about the positions, not about the alleles.

father_match <- father[which(father$chr_pos_38%in%mbw_match_no_dupl$chr_pos_38),]

#It seems we also have duplicates in father.
#We will take care of those, of course.

father_match_non_dupl <- father_match[which(duplicated(father_match$chr_pos_38) == FALSE),] #that is a perfect match!!!
father_match_dupl <- father_match[which(duplicated(father_match$chr_pos_38) == TRUE),] #no DUPLthat is a perfect match!!!

###################################
#Let's match the data for non-dupl#
###################################

#We are going to first do this one.
#Then take the duplicates.
#And then takes those that mismatch.
#Easy peasy lemon squeazy.

father_match_non_dupl <- father_match_non_dupl[order(match(father_match_non_dupl$chr_pos_38, mbw_match_no_dupl$chr_pos_38)),]

length(which(father_match_non_dupl$chr_pos_38 == mbw_match_no_dupl$chr_pos_38)) #all of them!

father_match_non_dupl$chr_pos_37 <- mbw_match_no_dupl$chr_pos_37 #awesome.

###############################################################################
#Using Pulit data to get SNPs in build 37 on those that could not be recovered#
###############################################################################

father_mismatch <- father[which(!(father$chr_pos_38%in%mbw_match_no_dupl$chr_pos_38)),]

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

bmi_match <- bmi[which(bmi$RSID%in%father_mismatch$rsID),]
father_mismatch_match <- father_mismatch[which(father_mismatch$rsID%in%bmi$RSID),]

#We have duplicates, probably in both.
#Since we want all the data on father_data..., we are gonna divide the data into father_dupl and father_non_dupl and match accordingly.

#In order to this properly we need to divide this between duplicated or not duplicated:

dupl <- father_mismatch_match$rsID[which(duplicated(father_mismatch_match$rsID) == TRUE)]

#Now we can get the chromosome and position for those that are duplicated and for those that are not.
#First the non-duplicated, because it is gonna be easier.

father_mismatch_match_not_dupl <- father_mismatch_match[which(!(father_mismatch_match$rsID%in%dupl)),]
bmi_match_not_dupl <- bmi_match[which(bmi_match$RSID%in%father_mismatch_match_not_dupl$rsID),] #so close!! Some variants, while not being duplicates in father, they are for BMI.

#That is not problem here, though, we just want the chr_pos37.

bmi_match_not_dupl <- bmi_match_not_dupl[which(duplicated(bmi_match_not_dupl$RSID) == FALSE),] #this did the trick.

#Let's match the data...

bmi_match_not_dupl <- bmi_match_not_dupl[order(match(bmi_match_not_dupl$RSID, father_mismatch_match_not_dupl$rsID)),]

#Let's check that it has worked out...

length(which(bmi_match_not_dupl$RSID != father_mismatch_match_not_dupl$rsID)) #all of them.

father_mismatch_match_not_dupl$chr_pos_37 <- paste("chr", bmi_match_not_dupl$CHR, ":", bmi_match_not_dupl$POS, sep = "")

#There are no duplicates, so there is no issue here my friend.

#Finally, let's match all of these dataframes with the chr_pos_37

father_mismatch_match_end <- father_mismatch_match_not_dupl

######################################################################
#And now let's gonna run liftover for the data that we could not find#
######################################################################

father_mismatch_2 <- father_mismatch[which(!(father_mismatch$rsID%in%father_mismatch_match_end$rsID)),] #this takes also the ones that do not have RSID.

dim(father_mismatch_2) #11141

#With 5M variants..., I think we can run this directly on liftover if we are lucky!

options(scipen = 100) #important to not make liftover crazy!

pos1 <- as.numeric(father_mismatch_2$Pos)-1
pos2 <- as.numeric(father_mismatch_2$Pos)+1

father_mismatch_2$build38_ranges <- paste(father_mismatch_2$Chr, ":", pos1, "-", pos2, sep = "")

fwrite(as.data.frame(father_mismatch_2$build38_ranges), "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/Father/father_build38.txt", sep = " ", col.names = FALSE)

#########################################
#Liftover failed: let's check the dealio#
#########################################

failed <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/Father/Failed_Conversions.txt")

#CHECKING if all is good:

father_mismatch_2_failed <- father_mismatch_2[which(father_mismatch_2$build38_ranges%in%failed$`#Partially deleted in new`),]

dim(father_mismatch_2_failed) #323 -> no problem.

father_mismatch_2_clean <- father_mismatch_2[which(!(father_mismatch_2$build38_ranges%in%failed$`#Partially deleted in new`)),]
father_mismatch_2_lost <- father_mismatch_2[which(father_mismatch_2$build38_ranges%in%failed$`#Partially deleted in new`),]

fwrite(as.data.frame(father_mismatch_2_clean$build38_ranges), "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/Father/father_build38_clean.txt", sep = " ", col.names = FALSE)

#####################################
#Let's retrieve the data in build 37#
#####################################

data_37 <- as.data.frame(read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/father/father_build37_raw.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

father_mismatch_2_clean$chr_pos_37 <- as.character(unlist(sapply(data_37$V1, clean_liftover_data)))
father_mismatch_2_lost$chr_pos_37 <- "-"

####################################################
#NOW WE CAN MERGE ALL THE DATA POSSIBLE, THANKS GOD#
####################################################

father_mismatch_2_clean <- father_mismatch_2_clean %>%
  dplyr::select(-c(build38_ranges))

father_mismatch_2_lost <- father_mismatch_2_lost %>%
  dplyr::select(-c(build38_ranges))

father_mismatch_full <- rbind(father_mismatch_match_end, father_mismatch_2_clean, father_mismatch_2_lost)

#We can merge the final data now:

father_curated_end <- rbind(father_mismatch_full, father_match_non_dupl) 

fwrite(father_curated_end, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/father/Birthweight_offspring_fathers2021_curated_FULL.txt")

##########
#Checking#
##########

check <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/OUTPUT/father/Birthweight2021_Curated.txt")
