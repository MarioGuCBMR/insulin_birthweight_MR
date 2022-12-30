##############
#INTRODUCTION#
##############

#This is code to run CAUSE for CIRadjISI -> MBWadjFBW

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

CIRadjISI <- fread("/home/projects/ku_00095/people/marure/BW_project/CURATED_DATA/DATA_4_CAUSE_&_LHC_MR/CIRadjISI_combined_Curated.txt")
MBWadjFBW <- fread("/home/projects/ku_00095/people/marure/BW_project/CURATED_DATA/DATA_4_CAUSE_&_LHC_MR/MBWadj_curated.txt")

CIRadjISI<- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/CIRadjISI/CIRadjISI_combined_Curated.txt")
MBWadjFBW <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/MBWAdj/MBWAdj_curated.txt")

###################
#Cleaning CIRadjISI data#
###################

#Let's first check the usual:

#EAF/MAF, info or sample size.

#NOTE we only have MAF.

summary(CIRadjISI$maf) #min 0.0110 = all good.

#We do not have the rest, so we skip that..., but we do have alleles! CAUSE will only detect
#standard alleles

yes_vect <- c("A", "G", "C", "T")

CIRadjISI <- CIRadjISI[which(CIRadjISI$effect_allele%in%yes_vect),]
CIRadjISI <- CIRadjISI[which(CIRadjISI$other_allele%in%yes_vect),]

#And now we remove the data in the MHC region:

options(scipen=999) #to avoid any weird transformations.

CIRadjISI$CHR <- as.numeric(as.character(unlist(sapply(CIRadjISI$chr_pos_37, chr_parser)))) #warnings due to NAs, no worries.
CIRadjISI$POS <- as.numeric(as.character(unlist(sapply(CIRadjISI$chr_pos_37, pos_parser)))) #warnings due to NAs, no worries.

CIRadjISI_mhc <- CIRadjISI[which(CIRadjISI$CHR == 6 & CIRadjISI$POS > 26000000 & CIRadjISI$POS < 34000000),] #so we are gonna approximate it this way.

summary(CIRadjISI_mhc$CHR) #perfect
summary(CIRadjISI_mhc$POS) #perfect

CIRadjISI_ <- CIRadjISI[which(!(CIRadjISI$chr_pos_37%in%CIRadjISI_mhc$chr_pos_37)),]

length(which(CIRadjISI_$CHR == 6 & CIRadjISI_$POS > 26000000 & CIRadjISI_$POS < 34000000)) #we erased them all!!

#Perfect!

rm(CIRadjISI)
rm(CIRadjISI_mhc)

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

#We have some variants in CIRadjISI that have no chr_pos...
#We are going to do the matching with chromosome and position, so we are gonna be really careful.

CIRadjISI_missing <- CIRadjISI_[which(CIRadjISI_$chr_pos_37 == "chrNA:NA"),] #347

#Let's check if the rsids are either in the MBWadjFBW rsIDs or in the merged rsIDs that they also share, which is supercool.

CIRadjISI_missing_match_1 <- CIRadjISI_missing[which(CIRadjISI_missing$snp%in%MBWadjFBW_$RSID),] #4! Not too many.
MBWadjFBW_missing_match_1 <- MBWadjFBW_[which(MBWadjFBW_$RSID%in%CIRadjISI_missing$snp),]

#Let's do the rest:

CIRadjISI_missing_rest <- CIRadjISI_missing[which(!(CIRadjISI_missing$snp%in%CIRadjISI_missing_match_1$snp)),]

#Let's play with merged variants with our reference panel because we are going to need it.

merged_rs <- fread("/home/projects/ku_00095/people/marure/BW_project/RAW_DATA/RsMergeArch.bcp.gz")

merged_rs$old_rs <- paste("rs", merged_rs$V1, sep = "")
merged_rs$new_rs <- paste("rs", merged_rs$V2, sep = "")

###################
#Now we match them#
###################

index_merged_old <- which(merged_rs$old_rs%in%CIRadjISI_missing_rest$snp) #Only 6! 
index_merged_new <- which(merged_rs$new_rs%in%CIRadjISI_missing_rest$snp) #We have 714. That means that we have A LOT of duplicates. All of them were too new it seems!! 

swapping_snps_old_2_new <- merged_rs[index_merged_old,]
swapping_snps_new_2_old <- merged_rs[index_merged_new,]

#Now we are gonna check which of these can be matched.

which(swapping_snps_old_2_new$new_rs%in%MBWadjFBW_$RSID) #1 needs to be changed.
which(swapping_snps_new_2_old$old_rs%in%MBWadjFBW_$RSID) #0: that means that we don't need to change these either.

mismatch <- swapping_snps_old_2_new[which(swapping_snps_old_2_new$new_rs%in%MBWadjFBW_$RSID),]

CIRadjISI_missing_match_2 <- CIRadjISI_missing_rest[which(CIRadjISI_missing_rest$snp == mismatch$old_rs),]

CIRadjISI_missing_match_2$snp <- mismatch$new_rs #CHANGED

MBWadjFBW_missing_match_2 <- MBWadjFBW_[which(MBWadjFBW_$RSID%in%CIRadjISI_missing_match_2$snp),]

#That means that we have to do 3 merges in the end.

CIRadjISI_not_missing <- CIRadjISI_[which(CIRadjISI_$chr_pos_37 != "chrNA:NA"),] 

#The first merge with the non_missing CIRadjISI_

library(cause)

#Careful, we need to change the columns of CIRadjISI_not_missing: they cannot have "snp".
#It confused CAUSE.

colnames(CIRadjISI_not_missing) <- c("SNP", "effect_allele", "other_allele",  "maf", "effect", "stderr", "pvalue", "chr_pos_37", "CHR", "POS")

X <- gwas_merge(CIRadjISI_not_missing, MBWadjFBW_, snp_name_cols = c("chr_pos_37", "chr_pos"), 
                beta_hat_cols = c("effect", "beta"), 
                se_cols = c("stderr", "se"), 
                A1_cols = c("effect_allele", "ea"), 
                A2_cols = c("other_allele", "nea"))

#And now we do it with the missing:

colnames(CIRadjISI_missing_match_1) <- c("SNP", "effect_allele", "other_allele",  "maf", "effect", "stderr", "pvalue", "chr_pos_37", "CHR", "POS")

X_missing_1 <- gwas_merge(CIRadjISI_missing_match_1, MBWadjFBW_missing_match_1, snp_name_cols = c("SNP", "RSID"), 
                beta_hat_cols = c("effect", "beta"), 
                se_cols = c("stderr", "se"), 
                A1_cols = c("effect_allele", "ea"), 
                A2_cols = c("other_allele", "nea"))

colnames(CIRadjISI_missing_match_2) <- c("SNP", "effect_allele", "other_allele",  "maf", "effect", "stderr", "pvalue", "chr_pos_37", "CHR", "POS")

X_missing_2 <- gwas_merge(CIRadjISI_missing_match_2, MBWadjFBW_missing_match_2, snp_name_cols = c("SNP", "RSID"), 
                        beta_hat_cols = c("effect", "beta"), 
                        se_cols = c("stderr", "se"), 
                        A1_cols = c("effect_allele", "ea"), 
                        A2_cols = c("other_allele", "nea"))

#Now let's retrieve the rsID for the chr_pos data set.

MBWadjFBW_match <- MBWadjFBW_[which(MBWadjFBW_$chr_pos%in%X$snp),]

#Let's check the rsID here.

MBWadjFBW_match <- MBWadjFBW_match[order(MBWadjFBW_match$RSID),]

head(MBWadjFBW_match, 100) #we do not have most of them
tail(MBWadjFBW_match) 

#Hence, we have to use the rsIDs from CIRadjISI.
#Far from perfect, but hey.

CIRadjISI_match <- CIRadjISI_[which(CIRadjISI_$chr_pos_37%in%X$snp),]

CIRadjISI_match <- CIRadjISI_match[order(match(CIRadjISI_match$chr_pos_37, X$snp)),]

length(CIRadjISI_match$chr_pos_37 == X$snp) #all of them

X$snp <- CIRadjISI_match$snp
X$p1 <- CIRadjISI_match$pvalue

#We have to do the same for CIRadjISIadjISI match 1:

CIRadjISI_missing_match_1_match <- CIRadjISI_missing_match_1[which(CIRadjISI_missing_match_1$SNP%in%X_missing_1$snp),]
CIRadjISI_missing_match_1_match <- CIRadjISI_missing_match_1_match[order(match(CIRadjISI_missing_match_1_match$SNP, X_missing_1$snp)),]

length(CIRadjISI_missing_match_1_match$SNP == X_missing_1$snp) #all of them

X_missing_1$p1 <- CIRadjISI_missing_match_1_match$pvalue

#And now we do the same for CIRadjISIadjISI match 2:

CIRadjISI_missing_match_2_match <- CIRadjISI_missing_match_2[which(CIRadjISI_missing_match_2$SNP%in%X_missing_2$snp),]
CIRadjISI_missing_match_2_match <- CIRadjISI_missing_match_2_match[order(match(CIRadjISI_missing_match_2_match$SNP, X_missing_2$snp)),]

length(CIRadjISI_missing_match_2_match$SNP == X_missing_2$snp) #all of them

X_missing_2$p1 <- CIRadjISI_missing_match_2_match$pvalue

#Finally, let's combine it with the missing ones:

X_end <- rbind(X, X_missing_1, X_missing_2)


############################################
#Calculating Nuisance and saving dataframes#
############################################

set.seed(100)
varlist <- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(X_end, file = "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IS/CIRadjISI/MBWadjFBW/dataframes/X_CIRadjISI_MBWadjFBW.rds")
saveRDS(params, file = "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IS/CIRadjISI/MBWadjFBW/dataframes/params_CIRadjISI_MBWadjFBW.rds")

########################
#PERFORMING LD-CLUMPING#
########################

library(tidyverse)

r2_thresh = 0.01
pval_thresh = 1e-3

X_clump <- X_end %>%
  rename(rsid = snp,
         pval = p1) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = "/services/tools/plink2/1.90beta6.24/plink", 
                     pop = "EUR")

top_vars <- X_clump$rsid

saveRDS(top_vars, "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IS/CIRadjISI/MBWadjFBW/dataframes/pruned_CIRadjISI_MBWadjFBW_res.rds")

###############
#RUNNING CAUSE#
###############

res <- cause(X=X_end, variants = top_vars, param_ests = params)

saveRDS(res, file = "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IS/CIRadjISI/MBWadjFBW/output/res_CIRadjISI_MBWadjFBW.rds")

res_strict <- cause(X=X_end, variants = top_vars, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IS/CIRadjISI/MBWadjFBW/output/res_CIRadjISI_MBWadjFBW_strict.rds")

print(res$elpd)

summary(res, ci_size=0.95)

#And now let's plot these bad fellas:

tiff("/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IS/CIRadjISI/MBWadjFBW/output/res_CIRadjISI_MBWadjFBW_plot.tiff", units="in", width=1200, height=500, res=300)
plot(res, type="data")
dev.off()
