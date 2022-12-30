##############
#INTRODUCTION#
##############

#This code is to curate the MBWadj data!

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

mbwadj <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Warrington_BW_data/Maternal_Effect_European_meta_NG2019.txt")

#############################
#Let's curate those variants#
#############################

mbwadj$chr_pos <- paste("chr", mbwadj$chr, ":", mbwadj$pos, sep = "")

mbwadj$ea <- as.character(unlist(sapply(mbwadj$ea, toupper)))
mbwadj$nea <- as.character(unlist(sapply(mbwadj$nea, toupper)))

###########
#CLEAN MAF#
###########

summary(mbwadj$eaf)

mbwadj_ <- mbwadj[which(mbwadj$eaf < 0.99 & mbwadj$eaf > 0.01),]

summary(mbwadj_$eaf)

###########
#CLEAN MHC#
###########

mbwadj_mhc <- mbwadj_[which(mbwadj_$chr == 6 & mbwadj_$pos >= 26000000 & mbwadj_$pos <= 34000000),]

summary(mbwadj_mhc$chr) #perfect
summary(mbwadj_mhc$pos) #perfect

mbwadj_mhc_clean <- mbwadj_[which(!(mbwadj_$chr_pos%in%mbwadj_mhc$chr_pos)),]

length(mbwadj_mhc$chr_pos)-length(mbwadj_$chr_pos) #perfect

#THE QC IS DONE!!!

fwrite(mbwadj_mhc_clean, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/MBWadj/MBWadj_curated.txt")
