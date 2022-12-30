##############
#INTRODUCTION#
##############

#This code is to curate the data from the Icelandic+EGG+UKBB.
#Here we are going to try to get the build38 to build37.
#Here goes nothing...

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

mbw <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Icelandic_BW_data/Birthweight_offspring_mothers2021_with_SE.txt")

##########################################
#Using Pulit data to get SNPs in build 37#
##########################################

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

bmi_match <- bmi[which(bmi$RSID%in%mbw$rsID),]
mbw_match <- mbw[which(mbw$rsID%in%bmi$RSID),]

#We have duplicates, probably in both.
#Since we want all the data on mbw_data..., we are gonna divide the data into mbw_dupl and mbw_non_dupl and match accordingly.

#In order to this properly we need to divide this between duplicated or not duplicated:

dupl <- mbw_match$rsID[which(duplicated(mbw_match$rsID) == TRUE)]

#Now we can get the chromosome and position for those that are duplicated and for those that are not.
#First the non-duplicated, because it is gonna be easier.

mbw_match_not_dupl <- mbw_match[which(!(mbw_match$rsID%in%dupl)),]
bmi_match_not_dupl <- bmi_match[which(bmi_match$RSID%in%mbw_match_not_dupl$rsID),] #so close!! Some variants, while not being duplicates in mbw, they are for BMI.

#That is not problem here, though, we just want the chr_pos37.

bmi_match_not_dupl <- bmi_match_not_dupl[which(duplicated(bmi_match_not_dupl$RSID) == FALSE),] #this did the trick.

#Let's match the data...

bmi_match_not_dupl <- bmi_match_not_dupl[order(match(bmi_match_not_dupl$RSID, mbw_match_not_dupl$rsID)),]

#Let's check that it has worked out...

length(which(bmi_match_not_dupl$RSID == mbw_match_not_dupl$rsID)) #all of them.

mbw_match_not_dupl$chr_pos_37 <- paste("chr", bmi_match_not_dupl$CHR, ":", bmi_match_not_dupl$POS, sep = "")

#Now we should go for the duplicates!!

mbw_match_dupl <- mbw_match[which(mbw_match$rsID%in%dupl),]
mbw_match_dupl <- mbw_match_dupl[order(mbw_match_dupl$rsID),] #this is gonna make us follow the data easier.

mbw_match_1 <- mbw_match_dupl[which(duplicated(mbw_match_dupl$rsID) == FALSE),]
mbw_match_2 <- mbw_match_dupl[which(duplicated(mbw_match_dupl$rsID) == TRUE),] #there are double duplicates, making this a bit more complicated...

#To do this OK, we are going to repeat this, taking into account that SNP as a special case.

special_snp <- mbw_match_2$rsID[which(duplicated(mbw_match_2$rsID) == TRUE)] #there are two different cases, not one with three copies!

mbw_match_dupl_special <- mbw_match_dupl[which(mbw_match_dupl$rsID%in%special_snp),] 

bmi_match_dupl_special <- bmi_match[which(bmi_match$RSID%in%mbw_match_dupl_special$rsID),] #here only 4 but that won't be a problem.

chr_pos_37_special <- c("chr18:3583873", "chr18:3583873", "chr18:3583873")

#And finally let's check what is going on here:

mbw_match_dupl_special$chr_pos_37 <- chr_pos_37_special

#And now we get the dupl, without the special ones.

mbw_match_dupl_normal <- mbw_match_dupl[which(!(mbw_match_dupl$rsID%in%special_snp)),]

#And finally, we get the duplicates in two sets:

mbw_match_1 <- mbw_match_dupl_normal[which(duplicated(mbw_match_dupl_normal$rsID) == FALSE),]
mbw_match_2 <- mbw_match_dupl_normal[which(duplicated(mbw_match_dupl_normal$rsID) == TRUE),] 

#Let's order the data so it is easier to retrieve the data:

mbw_match_1 <- mbw_match_1[order(match(mbw_match_1$rsID, mbw_match_2$rsID)),]

#Now let's get the BMI_match data

bmi_match_dupl_normal <- bmi_match[which(bmi_match$RSID%in%mbw_match_1$rsID),] #A difference of 3000!!! Are they all duplicates in BMI? Jesus.
bmi_match_dupl_normal <- bmi_match_dupl_normal[which(duplicated(bmi_match_dupl_normal$RSID) == FALSE),] #they were all duplicates. Lots of triallelic variants in there...
bmi_match_dupl_normal <- bmi_match_dupl_normal[order(match(bmi_match_dupl_normal$RSID, mbw_match_1$rsID)),]

#Let's check if all matches. If it does, we are almost done here.

length(which(mbw_match_1$rsID == mbw_match_2$rsID)) #they all match
length(which(mbw_match_2$rsID == bmi_match_dupl_normal$RSID)) #they all match

#Now we can put the chr_pos_37 in these fellas too:

mbw_match_1$chr_pos_37 <- paste("chr",bmi_match_dupl_normal$CHR, ":", bmi_match_dupl_normal$POS, sep = "")
mbw_match_2$chr_pos_37 <- paste("chr",bmi_match_dupl_normal$CHR, ":", bmi_match_dupl_normal$POS, sep = "")

#Finally, let's match all of these dataframes with the chr_pos_37

mbw_match_end <- rbind(mbw_match_1, mbw_match_2, mbw_match_dupl_special, mbw_match_not_dupl)

#And now we are gonna save these fellas, cuz it took a long time to generate them:

fwrite(mbw_match_end, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/OUTPUT/Birthweight_offspring_mothers2021_with_SE_chr_pos_37_with_Pulit.txt")

######################################################################
#And now let's gonna run liftover for the data that we could not find#
######################################################################

mbw_mismatch <- mbw[which(!(mbw$rsID%in%bmi$RSID)),] #this takes also the ones that do not have RSID.

dim(mbw_mismatch)[1] + dim(mbw_match)[1] #perfect split!!

#With 5M variants..., I think we can run this directly on liftover if we are lucky!

options(scipen = 100) #important to not make liftover crazy!

pos1 <- as.numeric(mbw_mismatch$Pos)-1
pos2 <- as.numeric(mbw_mismatch$Pos)+1

mbw_mismatch$build38_ranges <- paste(mbw_mismatch$Chr, ":", pos1, "-", pos2, sep = "")

fwrite(as.data.frame(mbw_mismatch$build38_ranges), "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/MBW/mbw_build38.txt", sep = " ", col.names = FALSE)

failures <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/MBW/Failed_Conversions.txt")

#Now let's repeat the conversion without the failed conversions:

mbw_mismatch_no_fail <- mbw_mismatch[which(!(mbw_mismatch$build38_ranges%in%failures$`#Deleted in new`)),]
mbw_mismatch_fail <- mbw_mismatch[which(mbw_mismatch$build38_ranges%in%failures$`#Deleted in new`),] #the numbers match with liftover!

fwrite(as.data.frame(mbw_mismatch_no_fail$build38_ranges), "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/MBW/mbw_build38_clean.txt", sep = " ", col.names = FALSE)

#The cleaned data was obtain in liftover and saved as mbw_build38_recovered.bed.
#Which is used in the next piece of code to obtain the positions in build 37.

#####################################
#Let's retrieve the data in build 38#
#####################################

data_37 <- as.data.frame(read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/MBW/mbw_build37_raw.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

mbw_mismatch_no_fail$chr_pos_37 <- data_37$V1
mbw_mismatch_fail$chr_pos_37 <- "-"

fwrite(mbw_mismatch_no_fail, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/MBW/mbw_build37_raw.txt")

############################################################################
#Now in Computerome we are going to clean this data with the following code#
############################################################################

data_37 <- fread("/home/projects/ku_00095/people/marure/BW_project/RAW_DATA/mbw_build37_raw.txt")

data_37$chr_pos_37_clean <- as.character(unlist(sapply(data_37$chr_pos_37,clean_liftover_data)))

fwrite(as.data.frame(data_37), "/home/projects/ku_00095/people/marure/BW_project/CURATED_DATA/mbw_build37_clean.txt", sep = " ", col.names = FALSE)

#################################################################
#And now we just need to load this data in our current R session#
#################################################################

mbw_mismatch_no_fail_clean <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/MBW/mbw_build37_clean.txt", fill = TRUE)

#FINALLY:

mbw_mismatch_no_fail$chr_pos_37_clean <- mbw_mismatch_no_fail_clean$V18

#It seems like this worked like a charm, but some of them could not be properly computed for some reason.
#Let's check how many and why:

length(which(mbw_mismatch_no_fail$chr_pos_37_clean == "")) #77830

#That is weird, but with it seems that the rest of conversions worked fine.
#That means that we can split up and try again:

mbw_mismatch_no_fail_missed <- mbw_mismatch_no_fail[which(mbw_mismatch_no_fail$chr_pos_37_clean == ""),]
mbw_mismatch_no_fail_OK <- mbw_mismatch_no_fail[which(mbw_mismatch_no_fail$chr_pos_37_clean != ""),]

#And finally...

mbw_mismatch_no_fail_missed$chr_pos_37_clean <- as.character(unlist(sapply(mbw_mismatch_no_fail_missed$chr_pos_37, clean_liftover_data)))

#It worked like a charm:

mbw_mismatch_no_fail_curated <- rbind(mbw_mismatch_no_fail_missed, mbw_mismatch_no_fail_OK)

mbw_mismatch_no_fail_curated <- mbw_mismatch_no_fail_curated %>%
  dplyr::select("Chr",              "Pos",              "rsID",             "A0",               "A1",               "IS-frq",          
         "IS-info",          "EGG-frq",          "Beta-A1",          "P",                "I2",               "P-het",           
         "SE",               "merged_rsID",      "chr_pos_38",       "build38_ranges",   "chr_pos_37_clean")

colnames(mbw_mismatch_no_fail_curated) <- c("Chr",              "Pos",              "rsID",             "A0",               "A1",               "IS-frq",          
                "IS-info",          "EGG-frq",          "Beta-A1",          "P",                "I2",               "P-het",           
                "SE",               "merged_rsID",      "chr_pos_38",       "build38_ranges",   "chr_pos_37")

####################################################
#NOW WE CAN MERGE ALL THE DATA POSSIBLE, THANKS GOD#
####################################################

mbw_mismatch_full <- rbind(mbw_mismatch_no_fail_curated, mbw_mismatch_fail)

#We can erase the ranges now:

mbw_mismatch_full_2 <- mbw_mismatch_full %>%
  dplyr::select(-c(build38_ranges))

mbw_curated_end <- rbind(mbw_mismatch_full_2, mbw_match_end) 

fwrite(mbw_curated_end, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/OUTPUT/MBW/Birthweight_offspring_mothers2021_curated_FULL.txt")

check <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/OUTPUT/MBW/Birthweight_offspring_mothers2021_curated_FULL.txt")
