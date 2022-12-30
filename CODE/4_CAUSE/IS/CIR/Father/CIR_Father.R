##############
#INTRODUCTION#
##############

#This is code to run CAUSE for CIR -> Father

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

CIR <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/CIR/CIR_combined_Curated.txt")
Father <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/father/Birthweight_offspring_fathers2021_curated_FULL.txt")

###################
#Cleaning CIR data#
###################

#Let's first check the usual:

#EAF/MAF, info or sample size.

#NOTE we only have MAF.

summary(CIR$maf) #min 0.0110 = all good.

#We do not have the rest, so we skip that..., but we do have alleles! CAUSE will only detect
#standard alleles

yes_vect <- c("A", "G", "C", "T")

CIR <- CIR[which(CIR$effect_allele%in%yes_vect),]
CIR <- CIR[which(CIR$other_allele%in%yes_vect),]

#And now we remove the data in the MHC region:

options(scipen=999) #to avoid any weird transformations.

CIR$CHR <- as.numeric(as.character(unlist(sapply(CIR$chr_pos_37, chr_parser)))) #warnings due to NAs, no worries.
CIR$POS <- as.numeric(as.character(unlist(sapply(CIR$chr_pos_37, pos_parser)))) #warnings due to NAs, no worries.

CIR_mhc <- CIR[which(CIR$CHR == 6 & CIR$POS > 26000000 & CIR$POS < 34000000),] #so we are gonna approximate it this way.

summary(CIR_mhc$CHR) #perfect
summary(CIR_mhc$POS) #perfect

CIR_ <- CIR[which(!(CIR$chr_pos_37%in%CIR_mhc$chr_pos_37)),]

length(which(CIR_$CHR == 6 & CIR_$POS > 26000000 & CIR_$POS < 34000000)) #we erased them all!!

#Perfect!

rm(CIR)
rm(CIR_mhc)

######################
#Cleaning Father data#
######################

summary(Father$`IS-frq`) #Min 1.01 and max 98.99. Perfect

#Now let's merge the data:

Father_maf <- Father #just to be consistent with nomenclature with other combinations.

rm(Father)

#Now let's check the info:

summary(Father_maf$`IS-info`) #all > 0.80

yes_vect <- c("A", "G", "C", "T")

Father_maf <- Father_maf[which(Father_maf$A0%in%yes_vect),]
Father_maf <- Father_maf[which(Father_maf$A1%in%yes_vect),]

Father_ <- Father_maf
rm(Father_maf)

########################################
#Matching data in the best way possible#
########################################

#We have some variants in CIR that have no chr_pos...
#We are going to do the matching with chromosome and position, so we are gonna be really careful.

CIR_missing <- CIR_[which(CIR_$chr_pos_37 == "chrNA:NA"),] #347

#Let's check if the rsids are either in the Father rsIDs or in the merged rsIDs that they also share, which is supercool.

CIR_missing_match_1 <- CIR_missing[which(CIR_missing$snp%in%Father_$rsID),] #150!!! We find one more than with the rest. Probably cuz there is no EGG data here to filter with.
Father_missing_match_1 <- Father_[which(Father_$rsID%in%CIR_missing$snp),]

#Now I understand, these are only available in build 38 or 18.
#Well that makes a lot of sense. 
#Let's see if we have some match with merged_rsID:

CIR_missing_match_2 <- CIR_missing[which(CIR_missing$snp%in%Father_$merged_rsID),] #0
Father_missing_match_2 <- Father_[which(Father_$merged_rsID%in%CIR_missing$snp),] #0

#That means that we have to do 2 merges:

CIR_not_missing <- CIR_[which(CIR_$chr_pos_37 != "chrNA:NA"),] 

#The first merge with the non_missing CIR_

library(cause)

#Careful, we need to change the columns of CIR_not_missing: they cannot have "snp".
#It confused CAUSE.

colnames(CIR_not_missing) <- c("SNP", "effect_allele", "other_allele",  "maf", "effect", "stderr", "pvalue", "chr_pos_37", "CHR", "POS")

X <- gwas_merge(CIR_not_missing, Father_, snp_name_cols = c("chr_pos_37", "chr_pos_37"), 
                beta_hat_cols = c("effect", "Beta-A1"), 
                se_cols = c("stderr", "SE"), 
                A1_cols = c("effect_allele", "A1"), 
                A2_cols = c("other_allele", "A0"))

#And now we do it with the missing:

colnames(CIR_missing_match_1) <- c("SNP", "effect_allele", "other_allele",  "maf", "effect", "stderr", "pvalue", "chr_pos_37", "CHR", "POS")

X_missing <- gwas_merge(CIR_missing_match_1, Father_missing_match_1, snp_name_cols = c("SNP", "rsID"), 
                beta_hat_cols = c("effect", "Beta-A1"), 
                se_cols = c("stderr", "SE"), 
                A1_cols = c("effect_allele", "A1"), 
                A2_cols = c("other_allele", "A0"))

#Now let's retrieve the rsID for the chr_pos data set.

Father_match <- Father_[which(Father_$chr_pos_37%in%X$snp),]

#Let's check the rsID here.

Father_match <- Father_match[order(Father_match$rsID),]

head(Father_match, 100) #we do not have most of them
tail(Father_match) 

#Hence, we have to use the rsIDs from CIR.
#Far from perfect, but hey.

CIR_match <- CIR_[which(CIR_$chr_pos_37%in%X$snp),]

CIR_match <- CIR_match[order(match(CIR_match$chr_pos_37, X$snp)),]

length(CIR_match$chr_pos_37 == X$snp) #all of them

X$snp <- CIR_match$snp

#Finally, let's combine it with the missing ones:

X_end <- rbind(X, X_missing)

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
