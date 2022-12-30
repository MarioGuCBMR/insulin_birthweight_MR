##############
#INTRODUCTION#
##############

#This code is to curate the fbwadj data!

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

fbwadj <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Warrington_BW_data/Fetal_Effect_European_meta_NG2019.txt")

#############################
#Let's curate those variants#
#############################

fbwadj$chr_pos <- paste("chr", fbwadj$chr, ":", fbwadj$pos, sep = "")

fbwadj$ea <- as.character(unlist(sapply(fbwadj$ea, toupper)))
fbwadj$nea <- as.character(unlist(sapply(fbwadj$nea, toupper)))

###########
#CLEAN MAF#
###########

summary(fbwadj$eaf)

fbwadj_ <- fbwadj[which(fbwadj$eaf < 0.99 & fbwadj$eaf > 0.01),]

summary(fbwadj_$eaf)

###########
#CLEAN MHC#
###########

fbwadj_mhc <- fbwadj_[which(fbwadj_$chr == 6 & fbwadj_$pos >= 26000000 & fbwadj_$pos <= 34000000),]

summary(fbwadj_mhc$chr) #perfect
summary(fbwadj_mhc$pos) #perfect

fbwadj_mhc_clean <- fbwadj_[which(!(fbwadj_$chr_pos%in%fbwadj_mhc$chr_pos)),]

length(fbwadj_mhc$chr_pos)-length(fbwadj_$chr_pos) #perfect

#THE QC IS DONE!!!

fwrite(fbwadj_mhc_clean, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/FBWadj/FBWAdj_curated.txt")
