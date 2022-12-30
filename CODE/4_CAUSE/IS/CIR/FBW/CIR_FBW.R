##############
#INTRODUCTION#
##############

#This is code to run CAUSE for CIR -> FBW

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
FBW <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/OUTPUT/FBW/Birthweight2021_Curated.txt")

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

###################
#Cleaning FBW data#
###################

summary(FBW$`IS-frq`) #Min 1.01 and max 98.99. Perfect
summary(as.numeric(FBW$`EGG-frq`)) #Lots of NAs, but we have frequencies that are a tad weird.

FBW_maf_no_NA <- FBW[which(is.na(as.numeric(FBW$`EGG-frq`)) == FALSE),] 
FBW_maf_NA <- FBW[which(is.na(as.numeric(FBW$`EGG-frq`)) == TRUE),]

#Now let's clean them:

summary(as.numeric(FBW_maf_no_NA$`EGG-frq`)) #CLEANED

FBW_maf_no_NA <- FBW_maf_no_NA[which(as.numeric(FBW_maf_no_NA$`EGG-frq`) > 1),]
FBW_maf_no_NA <- FBW_maf_no_NA[which(as.numeric(FBW_maf_no_NA$`EGG-frq`) < 99),]

summary(as.numeric(FBW_maf_no_NA$`EGG-frq`)) #CLEANED

#Now let's merge the data:

FBW_maf <- rbind(FBW_maf_no_NA, FBW_maf_NA)

rm(FBW_maf_NA)
rm(FBW_maf_no_NA)
rm(FBW)

#Now let's check the info:

summary(FBW_maf$`IS-info`) #all > 0.80

yes_vect <- c("A", "G", "C", "T")

FBW_maf <- FBW_maf[which(FBW_maf$A0%in%yes_vect),]
FBW_maf <- FBW_maf[which(FBW_maf$A1%in%yes_vect),]

FBW_ <- FBW_maf
rm(FBW_maf)

########################################
#Matching data in the best way possible#
########################################

#We have some variants in CIR that have no chr_pos...
#We are going to do the matching with chromosome and position, so we are gonna be really careful.

CIR_missing <- CIR_[which(CIR_$chr_pos_37 == "chrNA:NA"),] #347

#Let's check if the rsids are either in the FBW rsIDs or in the merged rsIDs that they also share, which is supercool.

CIR_missing_match_1 <- CIR_missing[which(CIR_missing$snp%in%FBW_$rsID),] #149!!!Quite many. That is great
FBW_missing_match_1 <- FBW_[which(FBW_$rsID%in%CIR_missing$snp),]

#Now I understand, these are only available in build 38 or 18.
#Well that makes a lot of sense. 
#Let's see if we have some match with merged_rsID:

CIR_missing_match_2 <- CIR_missing[which(CIR_missing$snp%in%FBW_$merged_rsID),] #0
FBW_missing_match_2 <- FBW_[which(FBW_$merged_rsID%in%CIR_missing$snp),] #0

#That means that we have to do 2 merges:

CIR_not_missing <- CIR_[which(CIR_$chr_pos_37 != "chrNA:NA"),] 

#The first merge with the non_missing CIR_

library(cause)

#Careful, we need to change the columns of CIR_not_missing: they cannot have "snp".
#It confused CAUSE.

colnames(CIR_not_missing) <- c("SNP", "effect_allele", "other_allele",  "maf", "effect", "stderr", "pvalue", "chr_pos_37", "CHR", "POS")

X <- gwas_merge(CIR_not_missing, FBW_, snp_name_cols = c("chr_pos_37", "chr_pos_37"), 
                beta_hat_cols = c("effect", "Beta-A1"), 
                se_cols = c("stderr", "SE"), 
                A1_cols = c("effect_allele", "A1"), 
                A2_cols = c("other_allele", "A0"))

#And now we do it with the missing:

colnames(CIR_missing_match_1) <- c("SNP", "effect_allele", "other_allele",  "maf", "effect", "stderr", "pvalue", "chr_pos_37", "CHR", "POS")

X_missing <- gwas_merge(CIR_missing_match_1, FBW_missing_match_1, snp_name_cols = c("SNP", "rsID"), 
                beta_hat_cols = c("effect", "Beta-A1"), 
                se_cols = c("stderr", "SE"), 
                A1_cols = c("effect_allele", "A1"), 
                A2_cols = c("other_allele", "A0"))

#Now let's retrieve the rsID for the chr_pos data set.

FBW_match <- FBW_[which(FBW_$chr_pos_37%in%X$snp),]

#Let's check the rsID here.

FBW_match <- FBW_match[order(FBW_match$rsID),]

head(FBW_match, 100) #we do not have most of them
tail(FBW_match) 

#Hence, we have to use the rsIDs from CIR.
#Far from perfect, but hey.

CIR_match <- CIR_[which(CIR_$chr_pos_37%in%X$snp),]

CIR_match <- CIR_match[order(match(CIR_match$chr_pos_37, X$snp)),]

length(CIR_match$chr_pos_37 == X$snp) #all of them

X$snp <- CIR_match$snp

#Be careful, and let's get the p-values on the merged dataframes for CIR.
#this is only necessary for the clumping:

X$p1 <- CIR_match$pvalue

#We have to do the same for CIRadjISI:

CIR_missing_match_1_match <- CIR_missing_match_1[which(CIR_missing_match_1$SNP%in%X_missing$snp),]
CIR_missing_match_1_match <- CIR_missing_match_1_match[order(match(CIR_missing_match_1_match$SNP, X_missing$snp)),]

length(CIR_missing_match_1_match$SNP == X_missing$snp) #all of them

X_missing$p1 <- CIR_missing_match_1_match$pvalue

#Finally, let's combine it with the missing ones:

X_end <- rbind(X, X_missing)

############################################
#Calculating Nuisance and saving dataframes#
############################################

set.seed(100)
varlist <- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(X_end, file = "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IS/CIR/FBW/dataframes/X_CIR_FBW.rds")
saveRDS(params, file = "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IS/CIR/FBW/dataframes/params_CIR_FBW.rds")

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

saveRDS(top_vars, "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IS/CIR/FBW/dataframes/pruned_CIR_FBW.res")

###############
#RUNNING CAUSE#
###############

res <- cause(X=X_end, variants = top_vars, param_ests = params)

saveRDS(res, file = "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IS/CIR/FBW/output/res_CIR_FBW.rds")

res_strict <- cause(X=X_end, variants = top_vars, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IS/CIR/FBW/output/res_CIR_FBW_strict.rds")

print(res$elpd)

summary(res, ci_size=0.95)

#And now let's plot these bad fellas:

tiff("/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IS/CIR/FBW/output/res_CIR_FBW_plot.tiff", units="in", width=1200, height=500, res=300)
plot(res, type="data")
dev.off()
