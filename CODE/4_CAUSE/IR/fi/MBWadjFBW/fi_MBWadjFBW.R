##############
#INTRODUCTION#
##############

#This is code to run CAUSE for fi -> MBWadjFBW

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

##############
#Loading data#
##############

fi <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/FI/fi_combined_curated.txt")
MBWadjFBW <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/MBWAdj/MBWAdj_curated.txt")

##################
#Cleaning fi data#
##################

#Let's first check the usual:
#The sample size seems weird, let's check whether we have no issue:

summary(fi$n) #GREAT. We could erase some of them, but who cares.

#Let's check the fi:

summary(fi$eaf) #PERFECT

#Let's run the summary:

yes_vect <- c("A", "G", "C", "T")

fi <- fi[which(fi$a2%in%yes_vect),]
fi <- fi[which(fi$a1%in%yes_vect),]

#Let's remove the fellas that are not in autosomal data:

fi$chr <- as.character(unlist(sapply(fi$chr_pos_37, chr_parser)))
fi$pos <- as.character(unlist(sapply(fi$chr_pos_37, pos_parser)))

#And now we remove the data in the MHC region:

fi_mhc <- fi[which(as.numeric(fi$chr) == 6 & as.numeric(fi$pos) > 26000000 & as.numeric(fi$pos) < 34000000),] #so we are gonna approximate it this way.

summary(as.numeric(fi_mhc$chr)) #perfect
summary(as.numeric(fi_mhc$pos)) #perfect

fi_ <- fi[which(!(fi$chr_pos_37%in%fi_mhc$chr_pos_37)),]

length(which(as.numeric(fi_$chr) == 6 & as.numeric(fi_$pos) > 26000000 & as.numeric(fi_$pos) < 34000000)) #we erased them all!!

#Perfect!

rm(fi)
rm(fi_mhc)

#########################
#Cleaning MBWadjFBW data#
#########################

summary(MBWadjFBW$eaf) #Perfect

#Now let's merge the data:

MBWadjFBW_maf <- MBWadjFBW

#No other info is here!!

MBWadjFBW_ <- MBWadjFBW_maf

rm(MBWadjFBW)
rm(MBWadjFBW_maf)

########################################
#Matching data in the best way possible#
########################################

#We have some variants in fi that have no chr_pos...
#We are going to do the matching with chromosome and position, so we are gonna be really careful.

fi_missing <- fi_[which(fi_$chr_pos_37 == "-"),] #996

#Do all of these have rsids???

length(which(str_detect(fi_missing$rsid, "rs") == TRUE)) #Nope, there are some that we do not!!

fi_missing_ok <- fi_missing[which(str_detect(fi_missing$rsid, "rs") == TRUE),]
fi_missing_fail <- fi_missing[which(str_detect(fi_missing$rsid, "rs") == FALSE),] #these are those in build hg18 that could not be liftovered. We remove them.

summary(fi_missing_fail$`p-value`) #they do not go through the P< 0.001 threshold so...

fi_missing <- fi_missing_ok

#Let's check if the rsids are either in the Father rsIDs or in the merged rsIDs that they also share, which is supercool.

fi_missing_match_1 <- fi_missing[which(fi_missing$rsid%in%MBWadjFBW_$RSID),] #22! Not too many
MBWadjFBW_missing_match_1 <- MBWadjFBW_[which(MBWadjFBW_$RSID%in%fi_missing$rsid),]

#Here we do not have merged RSIDs to be 100% sure.
#We need to strategy from fi, fiadjISI and ISI:

#Let's do the rest:

fi_missing_rest <- fi_missing[which(!(fi_missing$rsid%in%fi_missing_match_1$rsid)),]

#Let's play with merged variants with our reference panel because we are going to need it.

merged_rs <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Romain_Mario/RAW_DATA/MergedRs/RsMergeArch.bcp.gz")

merged_rs$old_rs <- paste("rs", merged_rs$V1, sep = "")
merged_rs$new_rs <- paste("rs", merged_rs$V2, sep = "")

###################
#Now we match them#
###################

index_merged_old <- which(merged_rs$old_rs%in%fi_missing_rest$rsid) #59 
index_merged_new <- which(merged_rs$new_rs%in%fi_missing_rest$rsid) #We have 1516. That means that we have A LOT of duplicates. All of them were too new it seems!! 

swapping_snps_old_2_new <- merged_rs[index_merged_old,]
swapping_snps_new_2_old <- merged_rs[index_merged_new,]

#Now we are gonna check which of these can be matched.

which(swapping_snps_old_2_new$new_rs%in%MBWadjFBW_$RSID) #several need to be changed.
which(swapping_snps_new_2_old$old_rs%in%MBWadjFBW_$RSID) #we need to change some of these fellas too.

mismatch <- swapping_snps_old_2_new[which(swapping_snps_old_2_new$new_rs%in%MBWadjFBW_$RSID),]

fi_missing_match_2 <- fi_missing_rest[which(fi_missing_rest$rsid%in%mismatch$old_rs),]
fi_missing_match_2 <- fi_missing_match_2[order(match(fi_missing_match_2$rsid, mismatch$old_rs)),]

length(which(fi_missing_match_2$rsid == mismatch$old_rs)) #all of them. We can change them.

fi_missing_match_2$rsid <- mismatch$new_rs #CHANGED

MBWadjFBW_missing_match_2 <- MBWadjFBW_[which(MBWadjFBW_$RSID%in%fi_missing_match_2$rsid),]

#Since here we do have some mismatches with MBWadjFBW, we are going to do it this way:

mismatch <- swapping_snps_new_2_old[which(swapping_snps_new_2_old$old_rs%in%MBWadjFBW_$RSID),]

fi_missing_match_3 <- fi_missing_rest[which(fi_missing_rest$rsid%in%mismatch$new_rs),]
fi_missing_match_3 <- fi_missing_match_3[order(match(fi_missing_match_3$rsid, mismatch$new_rs)),]

length(which(fi_missing_match_3$rsid == mismatch$new_rs)) #all of them. We can change them.

fi_missing_match_3$rsid <- mismatch$old_rs #CHANGED

MBWadjFBW_missing_match_3 <- MBWadjFBW_[which(MBWadjFBW_$RSID%in%fi_missing_match_3$rsid),]

################################
#AWESOME!!! THIS SHOULD DO IT!!#
################################

#That means that we have to do 3 merges in the end.

fi_not_missing <- fi_[which(fi_$chr_pos_37 != "-"),] 

#Now I understand, these are only available in build 38 or 18.
#Well that makes a lot of sense. 
#Let's see if we have some match with merged_rsID:

#The first merge with the non_missing fi_

library(cause)

#Careful, we need to change the columns of fi_not_missing: they cannot have "snp".
#It confused CAUSE.

X <- gwas_merge(fi_not_missing, MBWadjFBW_, snp_name_cols = c("chr_pos_37", "chr_pos"), 
                beta_hat_cols = c("beta", "beta"), 
                se_cols = c("se", "se"), 
                A1_cols = c("a2", "ea"), 
                A2_cols = c("a1", "nea"))

#And now we do it with the missing:

X_missing_1 <- gwas_merge(fi_missing_match_1, MBWadjFBW_missing_match_1, snp_name_cols = c("rsid", "RSID"), 
                          beta_hat_cols = c("beta", "beta"), 
                          se_cols = c("se", "se"), 
                          A1_cols = c("a2", "ea"), 
                          A2_cols = c("a1", "nea"))

X_missing_2 <- gwas_merge(fi_missing_match_2, MBWadjFBW_missing_match_2, snp_name_cols = c("rsid", "RSID"), 
                          beta_hat_cols = c("beta", "beta"), 
                          se_cols = c("se", "se"), 
                          A1_cols = c("a2", "ea"), 
                          A2_cols = c("a1", "nea"))

X_missing_3 <- gwas_merge(fi_missing_match_3, MBWadjFBW_missing_match_3, snp_name_cols = c("rsid", "RSID"), 
                          beta_hat_cols = c("beta", "beta"), 
                          se_cols = c("se", "se"), 
                          A1_cols = c("a2", "ea"), 
                          A2_cols = c("a1", "nea"))

#Now let's retrieve the rsID for the chr_pos data set.

fi_match <- fi_[which(fi_$chr_pos_37%in%X$snp),]

fi_match <- fi_match[order(match(fi_match$chr_pos_37, X$snp)),]

length(fi_match$chr_pos_37 == X$snp) #all of them

#Let's check how are we going to work with this:

check <- fi_match[order(fi_match$rsid),] #do we have any SNPs without RSID?

#In this case the issue comes from those that have as RSID those that have chr_pos_18-

head(check$rsid) #yes
tail(check$rsid) #ok

#Let's see how many:

length(which(str_detect(fi_match$rsid, "rs") == FALSE)) #quite many!

fi_match_no_rsid <- fi_match[which(str_detect(fi_match$rsid, "rs") == FALSE),]
fi_match_rsid <- fi_match[which(str_detect(fi_match$rsid, "rs") == TRUE),]

#Let's try to save these fellas:

MBWadjFBW_match_no_rsid <- MBWadjFBW_[which(MBWadjFBW_$chr_pos%in%fi_match_no_rsid$chr_pos_37),] #WE HAVE THEM.

#Now let's order them and check which ones can we retrieve there.

MBWadjFBW_match_no_rsid <- MBWadjFBW_match_no_rsid[order(match(MBWadjFBW_match_no_rsid$chr_pos, fi_match_no_rsid$chr_pos_37)),]

length(which(MBWadjFBW_match_no_rsid$chr_pos == fi_match_no_rsid$chr_pos_37)) #DONE.

#Let's be wary and get an index of all of those that have rsids:

index_rsid <- which(str_detect(MBWadjFBW_match_no_rsid$RSID, "rs") == TRUE) #ALL OF THEM

length(which(MBWadjFBW_match_no_rsid$chr_pos == fi_match_no_rsid$chr_pos_37)) #DONE. 

fi_match_no_rsid$rsid <- MBWadjFBW_match_no_rsid$RSID

#Now we can merge all of them without any issues:

fi_match_end <- rbind(fi_match_no_rsid, fi_match_rsid)

#And finally we can do this:

fi_match_end <- fi_match_end[order(match(fi_match_end$chr_pos_37, X$snp)),]

length(fi_match_end$chr_pos_37 == X$snp) #perfect

#We can do this!

X$snp <- fi_match_end$rsid
X$p1 <- fi_match_end$`p-value`

###############################################
#NOW LET'S DO ALL THE CASES FOR THE MISMATCHES#
###############################################

fi_missing_match_1_match <- fi_missing_match_1[which(fi_missing_match_1$rsid%in%X_missing_1$snp),]
fi_missing_match_1_match <- fi_missing_match_1_match[order(match(fi_missing_match_1_match$rsid, X_missing_1$snp)),]

length(fi_missing_match_1_match$rsid == X_missing_1$snp) #all of them

X_missing_1$p1 <- fi_missing_match_1_match$`p-value`

#And now we do the same for fiadjISI match 2:

fi_missing_match_2_match <- fi_missing_match_2[which(fi_missing_match_2$rsid%in%X_missing_2$snp),]
fi_missing_match_2_match <- fi_missing_match_2_match[order(match(fi_missing_match_2_match$rsid, X_missing_2$snp)),]

length(fi_missing_match_2_match$rsid == X_missing_2$snp) #all of them

X_missing_2$p1 <- fi_missing_match_2_match$`p-value`

#And 3:

fi_missing_match_3_match <- fi_missing_match_3[which(fi_missing_match_3$rsid%in%X_missing_3$snp),]
fi_missing_match_3_match <- fi_missing_match_3_match[order(match(fi_missing_match_3_match$rsid, X_missing_3$snp)),]

length(fi_missing_match_3_match$rsid == X_missing_3$snp) #all of them

X_missing_3$p1 <- fi_missing_match_3_match$`p-value`

X_end <- rbind(X, X_missing_1, X_missing_2, X_missing_3)

############################################
#Calculating Nuisance and saving dataframes#
############################################

set.seed(100)
varlist <- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(X_end, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/CAUSE/dataframes/X_BodyFatPerc_ST.rds")
saveRDS(params, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/CAUSE/dataframes/params_BodyFatPerc_ST_ORIGINAL.rds")

########################
#PERFORMING LD-CLUMPING#
########################

library(tidyverse)

variants <- X_end %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))

#1

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr1_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr1_AF0.05_snpdata.RDS")

pruned_1 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#2

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr2_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr2_AF0.05_snpdata.RDS")

pruned_2 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#3

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr3_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr3_AF0.05_snpdata.RDS")

pruned_3 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#4

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr4_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr4_AF0.05_snpdata.RDS")

pruned_4 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#5

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr5_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr5_AF0.05_snpdata.RDS")

pruned_5 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#6

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr6_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr6_AF0.05_snpdata.RDS")

pruned_6 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#7

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr7_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr7_AF0.05_snpdata.RDS")

pruned_7 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#8

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr8_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr8_AF0.05_snpdata.RDS")

pruned_8 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#9

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr9_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr9_AF0.05_snpdata.RDS")

pruned_9 <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))

#10

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr10_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr10_AF0.05_snpdata.RDS")

pruned_10 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#11

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr11_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr11_AF0.05_snpdata.RDS")

pruned_11 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#12

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr12_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr12_AF0.05_snpdata.RDS")

pruned_12 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#13

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr13_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr13_AF0.05_snpdata.RDS")

pruned_13 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#14

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr14_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr14_AF0.05_snpdata.RDS")

pruned_14 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#15

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr15_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr15_AF0.05_snpdata.RDS")

pruned_15 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#16

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr16_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr16_AF0.05_snpdata.RDS")

pruned_16 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))


#17

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr17_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr17_AF0.05_snpdata.RDS")

pruned_17 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#18

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr18_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr18_AF0.05_snpdata.RDS")

pruned_18 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#19

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr19_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr19_AF0.05_snpdata.RDS")

pruned_19 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))


#20

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr20_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr20_AF0.05_snpdata.RDS")

pruned_20 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

#21

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr21_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr21_AF0.05_snpdata.RDS")

pruned_21 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))
#22

#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_0.1.RDS?download=1", destfile = "C:/Users/zlc436/Desktop/Leisure_Project/chr22_AF0.05_0.1.RDS")
#download.file("https://zenodo.org/record/1464357/files/chr22_AF0.05_snpdata.RDS?download=1", destfile="C:/Users/zlc436/Desktop/Leisure_Project/chr22_AF0.05_snpdata.RDS")

ld <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr22_AF0.05_0.1.RDS")
snp_info <- readRDS("C:/Users/zlc436/Desktop/Leisure_Project/chr22_AF0.05_snpdata.RDS")

#variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
pruned_22 <- ld_prune(variants = variants, 
                      ld = ld, total_ld_variants = snp_info$SNP, 
                      pval_cols = c("pval1"), 
                      pval_thresh = c(1e-3))

pruned_all <- c(pruned_1, pruned_2, pruned_3, pruned_4, pruned_5, pruned_6, pruned_7, pruned_8, pruned_9, pruned_10, pruned_11, pruned_12, pruned_13, pruned_14, pruned_15, pruned_16, pruned_17, pruned_18, pruned_19, pruned_20, pruned_21, pruned_22)

################
#Check:
BodyFatPerc_22 <- BodyFatPerc_match[which(BodyFatPerc_match$RSID%in%pruned_22),]
summary(BodyFatPerc_22$CHR) #makes absolute sense.
################

BodyFatPerc_pruned <- BodyFatPerc_match[which(BodyFatPerc_match$SNP%in%pruned_all),]

write.table(BodyFatPerc_pruned, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/CAUSE/dataframes/BodyFatPerc_pruned_ST.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

######################################################

top_BodyFatPerc_pruned_df <- read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/CAUSE/dataframes/BodyFatPerc_pruned_ST.txt", header = TRUE, stringsAsFactors = FALSE)
top_BodyFatPerc_pruned_snps <- top_BodyFatPerc_pruned_df$SNP

###############
#RUNNING CAUSE#
###############

res <- cause(X=X_end, variants = top_BodyFatPerc_pruned_snps, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/CAUSE/dataframes/RES_BodyFatPerc_to_ST_ORIGINAL.rds")

res <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/CAUSE/dataframes/RES_BodyFatPerc_to_ST_ORIGINAL.rds")

res$elpd

#model1  model2  delta_elpd se_delta_elpd         z
#1    null sharing -143.037627    15.8592954 -9.019167
#2    null  causal -144.251749    16.1337127 -8.941014
#3 sharing  causal   -1.214122     0.3555328 -3.414936

summary(res)

#p-value testing that causal model is a better fit:  0.00032 
#Posterior medians and  95 % credible intervals:
#  model     gamma               eta                q                  
#[1,] "Sharing" NA                  "0.2 (0.17, 0.24)" "0.91 (0.76, 0.99)"
#[2,] "Causal"  "0.19 (0.12, 0.25)" "0 (-0.33, 0.31)"  "0.18 (0, 0.86)"   

plot(res, type="data")
