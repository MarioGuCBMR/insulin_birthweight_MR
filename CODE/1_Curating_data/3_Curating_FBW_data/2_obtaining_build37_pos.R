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

fbw <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Icelandic_BW_data/Birthweight2021_with_SE.txt")

mbw_curated <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/OUTPUT/MBW/Birthweight_offspring_mothers2021_curated_FULL.txt")

###################
#Matching MBW data#
###################

mbw_match <- mbw_curated[which(mbw_curated$chr_pos_38%in%fbw$chr_pos_38),]
mbw_match_no_dupl <- mbw_match[which(duplicated(mbw_match$chr_pos_38) == FALSE),]

#We have duplicates, but we do not care about that!
#We care about the positions, not about the alleles.

fbw_match <- fbw[which(fbw$chr_pos_38%in%mbw_match_no_dupl$chr_pos_38),]

#It seems we also have duplicates in fbw.
#We will take care of those, of course.

fbw_match_non_dupl <- fbw_match[which(duplicated(fbw_match$chr_pos_38) == FALSE),] #that is a perfect match!!!
fbw_match_dupl <- fbw_match[which(duplicated(fbw_match$chr_pos_38) == TRUE),] #that is a perfect match!!!

#Let's first divide the data like this:

fbw_match_non_dupl_1 <- fbw_match_dupl[which(duplicated(fbw_match_dupl$chr_pos_38) == FALSE),]
fbw_match_dupl_1 <- fbw_match_dupl[which(duplicated(fbw_match_dupl$chr_pos_38) == TRUE),]

#Now let's repeat:

fbw_match_non_dupl_2 <- fbw_match_dupl_1[which(duplicated(fbw_match_dupl_1$chr_pos_38) == FALSE),]
fbw_match_dupl_2 <- fbw_match_dupl_1[which(duplicated(fbw_match_dupl_1$chr_pos_38) == TRUE),] 

#This ends here. This means that fbw_match_dupl_1 is the last one we need! 
#We can use the non_dupl_1 and non_dupl_2.

############################################
#Let's match the data for non-dupl and dupl#
############################################

#We are going to first do this one.
#Then take the duplicates.
#And then takes those that mismatch.
#Easy peasy lemon squeazy.

fbw_match_non_dupl <- fbw_match_non_dupl[order(match(fbw_match_non_dupl$chr_pos_38, mbw_match_no_dupl$chr_pos_38)),]

length(which(fbw_match_non_dupl$chr_pos_38 == mbw_match_no_dupl$chr_pos_38)) #all of them!

fbw_match_non_dupl$chr_pos_37 <- mbw_match_no_dupl$chr_pos_37 #awesome.

#Now let's take care of the duplicates.

mbw_4_dupl_1 <- mbw_match_no_dupl[which(mbw_match_no_dupl$chr_pos_38%in%fbw_match_non_dupl_1$chr_pos_38),] #perfect match
mbw_4_dupl_2 <- mbw_match_no_dupl[which(mbw_match_no_dupl$chr_pos_38%in%fbw_match_non_dupl_2$chr_pos_38),] #perfect match

#Okay, we only need to order them, take the chr_pos and leave~

fbw_match_non_dupl_1 <- fbw_match_non_dupl_1[order(match(fbw_match_non_dupl_1$chr_pos_38, mbw_4_dupl_1$chr_pos_38)),]
fbw_match_non_dupl_2 <- fbw_match_non_dupl_2[order(match(fbw_match_non_dupl_2$chr_pos_38, mbw_4_dupl_2$chr_pos_38)),]

#Let's bind the data now, since they are all ordered they are bound to have no problem:

fbw_dupl_end <- rbind(fbw_match_non_dupl_1, fbw_match_non_dupl_2)

mbw_dupl_end <- rbind(mbw_4_dupl_1, mbw_4_dupl_2)

#Checking that all is good:

length(which(fbw_dupl_end$chr_pos_38 == mbw_dupl_end$chr_pos_38)) #perfect match

fbw_dupl_end$chr_pos_37 <- mbw_dupl_end$chr_pos_37

###############################################################################
#Using Pulit data to get SNPs in build 37 on those that could not be recovered#
###############################################################################

fbw_mismatch <- fbw[which(!(fbw$chr_pos_38%in%mbw_match_no_dupl$chr_pos_38)),]

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

bmi_match <- bmi[which(bmi$RSID%in%fbw_mismatch$rsID),]
fbw_mismatch_match <- fbw_mismatch[which(fbw_mismatch$rsID%in%bmi$RSID),]

#We have duplicates, probably in both.
#Since we want all the data on fbw_data..., we are gonna divide the data into fbw_dupl and fbw_non_dupl and match accordingly.

#In order to this properly we need to divide this between duplicated or not duplicated:

dupl <- fbw_mismatch_match$rsID[which(duplicated(fbw_mismatch_match$rsID) == TRUE)]

#Now we can get the chromosome and position for those that are duplicated and for those that are not.
#First the non-duplicated, because it is gonna be easier.

fbw_mismatch_match_not_dupl <- fbw_mismatch_match[which(!(fbw_mismatch_match$rsID%in%dupl)),]
bmi_match_not_dupl <- bmi_match[which(bmi_match$RSID%in%fbw_mismatch_match_not_dupl$rsID),] #so close!! Some variants, while not being duplicates in fbw, they are for BMI.

#That is not problem here, though, we just want the chr_pos37.

bmi_match_not_dupl <- bmi_match_not_dupl[which(duplicated(bmi_match_not_dupl$RSID) == FALSE),] #this did the trick.

#Let's match the data...

bmi_match_not_dupl <- bmi_match_not_dupl[order(match(bmi_match_not_dupl$RSID, fbw_mismatch_match_not_dupl$rsID)),]

#Let's check that it has worked out...

length(which(bmi_match_not_dupl$RSID == fbw_mismatch_match_not_dupl$rsID)) #all of them.

fbw_mismatch_match_not_dupl$chr_pos_37 <- paste("chr", bmi_match_not_dupl$CHR, ":", bmi_match_not_dupl$POS, sep = "")

#Now we should go for the duplicates!!

fbw_mismatch_match_dupl <- fbw_mismatch_match[which(fbw_mismatch_match$rsID%in%dupl),]
fbw_mismatch_match_dupl <- fbw_mismatch_match_dupl[order(fbw_mismatch_match_dupl$rsID),] 

fbw_mismatch_match_1 <- fbw_mismatch_match_dupl[which(duplicated(fbw_mismatch_match_dupl$rsID) == FALSE),]
fbw_mismatch_match_2 <- fbw_mismatch_match_dupl[which(duplicated(fbw_mismatch_match_dupl$rsID) == TRUE),] #no double duplicates. This one is gonna be easy.

#We just need to match bmi with these fellas once:

bmi_match_dupl <- bmi_match[which(bmi_match$RSID%in%fbw_mismatch_match_1$rsID),]
bmi_match_dupl <- bmi_match_dupl[which(duplicated(bmi_match_dupl$RSID) == FALSE),]
bmi_match_dupl <- bmi_match_dupl[order(match(bmi_match_dupl$RSID, fbw_mismatch_match_1$rsID)),]

#Let's check if all matches. If it does, we are almost done here.

length(which(fbw_mismatch_match_1$rsID == fbw_mismatch_match_2$rsID)) #they all match
length(which(fbw_mismatch_match_2$rsID == bmi_match_dupl$RSID)) #they all match

#Now we can put the chr_pos_37 in these fellas too:

fbw_mismatch_match_1$chr_pos_37 <- paste("chr",bmi_match_dupl$CHR, ":", bmi_match_dupl$POS, sep = "")
fbw_mismatch_match_2$chr_pos_37 <- paste("chr",bmi_match_dupl$CHR, ":", bmi_match_dupl$POS, sep = "")

#Finally, let's match all of these dataframes with the chr_pos_37

fbw_mismatch_match_end <- rbind(fbw_mismatch_match_1, fbw_mismatch_match_2, fbw_mismatch_match_not_dupl)

######################################################################
#And now let's gonna run liftover for the data that we could not find#
######################################################################

fbw_mismatch_2 <- fbw_mismatch[which(!(fbw_mismatch$rsID%in%fbw_mismatch_match_end$rsID)),] #this takes also the ones that do not have RSID.

dim(fbw_mismatch_2) #24510

#With 5M variants..., I think we can run this directly on liftover if we are lucky!

options(scipen = 100) #important to not make liftover crazy!

pos1 <- as.numeric(fbw_mismatch_2$Pos)-1
pos2 <- as.numeric(fbw_mismatch_2$Pos)+1

fbw_mismatch_2$build38_ranges <- paste(fbw_mismatch_2$Chr, ":", pos1, "-", pos2, sep = "")

fwrite(as.data.frame(fbw_mismatch_2$build38_ranges), "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/FBW/fbw_build38.txt", sep = " ", col.names = FALSE)

#####################################
#Let's retrieve the data in build 37#
#####################################

data_37 <- as.data.frame(read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/FBW/fbw_build37_raw.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

#We miss one and we know which on is:

fbw_mismatch_2_fail <- fbw_mismatch_2[which(fbw_mismatch_2$build38_ranges == "chr9:66337936-66337938"),]
fbw_mismatch_2_ok <- fbw_mismatch_2[which(!(fbw_mismatch_2$build38_ranges == "chr9:66337936-66337938")),]

#Let's run liftover again, just in case:

fwrite(as.data.frame(fbw_mismatch_2_ok$build38_ranges), "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/FBW/fbw_build38_clean.txt", sep = " ", col.names = FALSE)

#Let's retrieve the data one more time with feeling:

data_37 <- as.data.frame(read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/INPUT/FBW/fbw_build37_raw_clean.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

fbw_mismatch_2_ok$chr_pos_37 <- as.character(unlist(sapply(data_37$V1,clean_liftover_data)))
fbw_mismatch_2_fail$chr_pos_37 <- "-"

fbw_mismatch_2 <- rbind(fbw_mismatch_2_ok, fbw_mismatch_2_fail)

####################################################
#NOW WE CAN MERGE ALL THE DATA POSSIBLE, THANKS GOD#
####################################################

fbw_mismatch_2_clean <- fbw_mismatch_2 %>%
  dplyr::select(-c(build38_ranges))

fbw_mismatch_full <- rbind(fbw_mismatch_match_end, fbw_mismatch_2_clean)

#We can merge the final data now:

fbw_curated_end <- rbind(fbw_mismatch_full, fbw_dupl_end, fbw_match_non_dupl) 

fwrite(fbw_curated_end, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/OUTPUT/FBW/Birthweight2021_Curated.txt")

##########
#Checking#
##########

check <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/OUTPUT/FBW/Birthweight2021_Curated.txt")
