##############
#INTRODUCTION#
##############

#This code will serve to parse the data of FIadjBMI used in for the GRS and 2SMR.

###################
#LOADING FUNCTIONS#
###################

library(tidyverse)
library(xlsx)
library(data.table)

memory.size(8000000000)

################
#LOAD FUNCTIONS#
################

load_independent_snps <- function(path_){
  #This function takes a path and the initial part of the name of the chromosome data (the last thing that needs to be written is chr.)
  #And loads all the data in the folder, chromosome by chromosome and adds the SNPs in a vector.
  
  #path_ <- path_fi_combined_lipid_glgc
  
  independent_snp_lax <- c()
  
  for(i in seq(1,22)){
    
    path_end <- paste(path_, i, ".clumped", sep = "")
    
    tmp_ <- fread(path_end)
    
    independent_snp_lax <- c(independent_snp_lax, tmp_$SNP)
    
  }
  
  return(independent_snp_lax)
  
}

parse_rsid <- function(rsid){
  
  rsid_ <- strsplit(rsid, ":")[[1]][1] #this works for those that do not present the weird formatting :)
  
  return(rsid_)
  
}

##############
#LOADING DATA#
##############

fiadjbmi <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Glycemic_traits_raw_data/MAGIC1000G_FI_EUR.tsv.gz")
hb1ac <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Glycemic_traits_raw_data/MAGIC1000G_HbA1c_EUR.tsv.gz")

######################################################################
#STEP 1: OBTAIN RSIDS FOR ALL SNPS THAT ARE SHARED BETWEEN DATAFRAMES#
######################################################################

#Let's match them with chr:pos.

fiadjbmi$chr_pos <- paste(fiadjbmi$chromosome, ":", fiadjbmi$base_pair_location, sep = "")
hb1ac$chr_pos <- paste(hb1ac$chromosome, ":", hb1ac$base_pair_location, sep = "")

fiadjbmi_match <- fiadjbmi[which(fiadjbmi$chr_pos%in%hb1ac$chr_pos),]
#hb1ac_match <- hb1ac[which(hb1ac$chr_pos%in%fiadjbmi$chr_pos),]

#The duplicates battle is going to be fearsome in this one.

##########################
#FIadjBMI: non-duplicates#
##########################

fiadjbmi_match_non_dupl <- fiadjbmi_match[which(duplicated(fiadjbmi_match$chr_pos) == FALSE),]

hb1ac_match_non_dupl <- hb1ac[which(hb1ac$chr_pos%in%fiadjbmi_match_non_dupl$chr_pos),]
hb1ac_match_non_dupl <- hb1ac_match_non_dupl[which(duplicated(hb1ac_match_non_dupl$chr_pos) == FALSE),]

#Absolute perfect match!!

fiadjbmi_match_non_dupl <- fiadjbmi_match_non_dupl[order(match(fiadjbmi_match_non_dupl$chromosome, hb1ac_match_non_dupl$chr_pos)),]

length(which(fiadjbmi_match_non_dupl$chr_pos == hb1ac_match_non_dupl$chr_pos)) #all of them

fiadjbmi_match_non_dupl$RSID <- hb1ac_match_non_dupl$variant #we may not have all RSIDs, but we win most of them

#Let's check how many of these do not include RSID:

length(which(str_detect(fiadjbmi_match_non_dupl$RSID, "rs") == FALSE)) #there are some clearly missing.

#We will do what we can with these fellas.

#In any case, let's get the duplicated:

fiadjbmi_match_dupl <- fiadjbmi_match[which(duplicated(fiadjbmi_match$chr_pos) == TRUE),]

#We have double duplicates, of course:

fiadjbmi_match_dupl_1 <- fiadjbmi_match_dupl[which(duplicated(fiadjbmi_match_dupl$chr_pos) == FALSE),]
fiadjbmi_match_dupl_2 <- fiadjbmi_match_dupl[which(duplicated(fiadjbmi_match_dupl$chr_pos) == TRUE),]

fiadjbmi_match_ref_4_dupl_1 <- fiadjbmi_match_non_dupl[which(fiadjbmi_match_non_dupl$chr_pos%in%fiadjbmi_match_dupl_1$chr_pos),]
fiadjbmi_match_ref_4_dupl_2 <- fiadjbmi_match_non_dupl[which(fiadjbmi_match_non_dupl$chr_pos%in%fiadjbmi_match_dupl_2$chr_pos),]

fiadjbmi_match_dupl_1 <- fiadjbmi_match_dupl_1[order(match(fiadjbmi_match_dupl_1$chr_pos, fiadjbmi_match_ref_4_dupl_1$chr_pos)),]
fiadjbmi_match_dupl_2 <- fiadjbmi_match_dupl_1[order(match(fiadjbmi_match_dupl_2$chr_pos, fiadjbmi_match_ref_4_dupl_2$chr_pos)),]

length(which(fiadjbmi_match_dupl_1$chr_pos == fiadjbmi_match_dupl_1$chr_pos)) #perfect
length(which(fiadjbmi_match_dupl_2$chr_pos == fiadjbmi_match_dupl_2$chr_pos)) #perfect

fiadjbmi_match_dupl_1$RSID <- fiadjbmi_match_ref_4_dupl_1$RSID
fiadjbmi_match_dupl_2$RSID <- fiadjbmi_match_ref_4_dupl_2$RSID

fiadjbmi_match_end <- rbind(fiadjbmi_match_non_dupl, fiadjbmi_match_dupl_1, fiadjbmi_match_dupl_2)

#####################################################################################
#STEP 2: let's try to recover those that are not matching between them with BMI data#
#####################################################################################

fiadjbmi_mismatch <- fiadjbmi[which(!(fiadjbmi$chr_pos%in%hb1ac$chr_pos)),]

bmi <- fread("C:/Users/zlc436/Desktop/PULIT_data_reservoir/Combined/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")

#We need to be careful, some SNPs from Pulit et al do not present correct chr_pos.
#The mistake is in those without INFO data, thus, the ones in GIANT data.
#We are gonna remove those.

bmi <- bmi[which(is.na(bmi$INFO) == FALSE),]

summary(bmi$INFO) #all good
summary(bmi$Freq_Tested_Allele) #all good
summary(bmi$CHR) #all good
summary(bmi$N) #all good

bmi$chr_pos <- paste(bmi$CHR, ":", bmi$POS, sep = "")

#bmi$RSID <- as.character(unlist(sapply(bmi$SNP, parse_rsid)))

bmi_match <- bmi[which(bmi$chr_pos%in%fiadjbmi_mismatch$chr_pos),] #let's check if it is worthy.

#helluva worthy:

fiadjbmi_mismatch_in_bmi <- fiadjbmi_mismatch[which(fiadjbmi_mismatch$chr_pos%in%bmi_match$chr_pos),]

#Do we have duplicates?

fiadjbmi_mismatch_in_bmi_non_dupl <- fiadjbmi_mismatch_in_bmi[which(duplicated(fiadjbmi_mismatch_in_bmi$chr_pos) == FALSE),]

bmi_match_non_dupl <- bmi_match[which(bmi_match$chr_pos%in%fiadjbmi_mismatch_in_bmi_non_dupl$chr_pos),]

bmi_match_non_dupl <- bmi_match_non_dupl[which(duplicated(bmi_match_non_dupl$chr_pos) == FALSE),]

#Now let's order these fellas:

bmi_match_non_dupl <- bmi_match_non_dupl[order(match(bmi_match_non_dupl$chr_pos, fiadjbmi_mismatch_in_bmi_non_dupl$chr_pos)),]

length(which(bmi_match_non_dupl$chr_pos == fiadjbmi_mismatch_in_bmi_non_dupl$chr_pos)) #perfect

bmi_match_non_dupl$RSID <- as.character(unlist(sapply(bmi_match_non_dupl$SNP, parse_rsid)))

fiadjbmi_mismatch_in_bmi_non_dupl$RSID <- bmi_match_non_dupl$RSID

#Now let's go for the duplicates:

fiadjbmi_mismatch_in_bmi_dupl <- fiadjbmi_mismatch_in_bmi[which(duplicated(fiadjbmi_mismatch_in_bmi$chr_pos) == TRUE),]
fiadjbmi_mismatch_in_bmi_non_dupl_4_ref <- fiadjbmi_mismatch_in_bmi_non_dupl[which(fiadjbmi_mismatch_in_bmi_non_dupl$chr_pos%in%fiadjbmi_mismatch_in_bmi_dupl$chr_pos),]

fiadjbmi_mismatch_in_bmi_non_dupl_4_ref <- fiadjbmi_mismatch_in_bmi_non_dupl_4_ref[order(match(fiadjbmi_mismatch_in_bmi_non_dupl_4_ref$chr_pos, fiadjbmi_mismatch_in_bmi_dupl$chr_pos)),]

length(which(fiadjbmi_mismatch_in_bmi_non_dupl_4_ref$chr_pos == fiadjbmi_mismatch_in_bmi_dupl$chr_pos))

fiadjbmi_mismatch_in_bmi_dupl$RSID <- fiadjbmi_mismatch_in_bmi_non_dupl_4_ref$RSID

fiadjbmi_mismatch_in_bmi_end <- rbind(fiadjbmi_mismatch_in_bmi_non_dupl, fiadjbmi_mismatch_in_bmi_dupl)

#########
#AWESOME#
#########

#Just to check we are gonna see how many are we missing, but we are not going to do 
#anything.
#If we solve them..., we are solving them when doing specific merges.

fiadjbmi_mismatch_end <- fiadjbmi_mismatch[which(!(fiadjbmi_mismatch$chr_pos%in%bmi_match$chr_pos)),] #around a million!!

#Still, not too much we can do there.
#Let's save the data and do what we can!!

fiadjbmi_mismatch_end$RSID <- "-"

####################################################################
#Let's combine these fellas that we have obtained and save the data#
####################################################################

fiadjbmi_curated_end <- rbind(fiadjbmi_mismatch_end, fiadjbmi_mismatch_in_bmi_end, fiadjbmi_match_end)

fwrite(fiadjbmi_curated_end, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/FIadjBMI/FIadjBMI_curated.txt")
