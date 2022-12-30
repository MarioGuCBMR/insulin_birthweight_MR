##############
#INTRODUCTION#
##############

#This is code to run CAUSE for fi -> MBW

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

fi <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/FI/FI_combined_curated.txt")
MBW <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/OUTPUT/MBW/Birthweight_offspring_mothers2021_curated_FULL.txt")

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

###################
#Cleaning MBW data#
###################

summary(MBW$`IS-frq`) #Min 1.01 and max 98.99. Perfect
summary(as.numeric(MBW$`EGG-frq`)) #Lots of NAs, but we the frequencies are OK.

MBW_maf_no_NA <- MBW[which(is.na(as.numeric(MBW$`EGG-frq`)) == FALSE),] 
MBW_maf_NA <- MBW[which(is.na(as.numeric(MBW$`EGG-frq`)) == TRUE),]

#Now let's clean them:

summary(as.numeric(MBW_maf_no_NA$`EGG-frq`)) #CLEANED before

#Now let's merge the data:

MBW_maf <- rbind(MBW_maf_no_NA, MBW_maf_NA)

rm(MBW_maf_NA)
rm(MBW_maf_no_NA)
rm(MBW)

#Now let's check the info:

summary(MBW_maf$`IS-info`) #all > 0.80

yes_vect <- c("A", "G", "C", "T")

MBW_maf <- MBW_maf[which(MBW_maf$A0%in%yes_vect),]
MBW_maf <- MBW_maf[which(MBW_maf$A1%in%yes_vect),]

MBW_ <- MBW_maf
rm(MBW_maf)

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

#Let's check if the rsids are either in the MBW rsIDs or in the merged rsIDs that they also share, which is supercool.

fi_missing_match_1 <- fi_missing[which(fi_missing$rsid%in%MBW_$rsID),] #311!!!Quite many. That is great
MBW_missing_match_1 <- MBW_[which(MBW_$rsID%in%fi_missing$rsid),]

#Now I understand, these are only available in build 38 or 18.
#Well that makes a lot of sense. 
#Let's see if we have some match with merged_rsID:

fi_missing_match_2 <- fi_missing[which(fi_missing$rsid%in%MBW_$merged_rsID),] #3
MBW_missing_match_2 <- MBW_[which(MBW_$merged_rsID%in%fi_missing$rsid),] #3

#That means that we have to do 2 merges:

fi_not_missing <- fi_[which(fi_$chr_pos_37 != "-"),] 

#The first merge with the non_missing fi_

library(cause)

#Careful, we need to change the columns of fi_not_missing: they cannot have "snp".
#It confused CAUSE.

X <- gwas_merge(fi_not_missing, MBW_, snp_name_cols = c("chr_pos_37", "chr_pos_37"), 
                beta_hat_cols = c("beta", "Beta-A1"), 
                se_cols = c("se", "SE"), 
                A1_cols = c("a2", "A1"), 
                A2_cols = c("a1", "A0"))

#And now we do it with the missing:

X_missing_1 <- gwas_merge(fi_missing_match_1, MBW_missing_match_1, snp_name_cols = c("rsid", "rsID"), 
                          beta_hat_cols = c("beta", "Beta-A1"), 
                          se_cols = c("se", "SE"), 
                          A1_cols = c("a2", "A1"), 
                          A2_cols = c("a1", "A0"))

X_missing_2 <- gwas_merge(fi_missing_match_2, MBW_missing_match_2, snp_name_cols = c("rsid", "rsID"), 
                          beta_hat_cols = c("beta", "Beta-A1"), 
                          se_cols = c("se", "SE"), 
                          A1_cols = c("a2", "A1"), 
                          A2_cols = c("a1", "A0"))


#Now let's retrieve the rsID for the chr_pos data set.

MBW_match <- MBW_[which(MBW_$chr_pos_37%in%X$snp),]

#Let's check the rsID here.

MBW_match <- MBW_match[order(MBW_match$rsID),]

head(MBW_match, 100) #we do not have most of them
tail(MBW_match) 

#Hence, we have to use the rsIDs from fi.
#Far from perfect, but hey.

fi_match <- fi_[which(fi_$chr_pos_37%in%X$snp),]

fi_match <- fi_match[order(match(fi_match$chr_pos_37, X$snp)),]

length(fi_match$chr_pos_37 == X$snp) #all of them

X$snp <- fi_match$rsid

#Be careful, and let's get the p-values on the merged dataframes for fi.
#this is only necessary for the clumping:

X$p1 <- fi_match$`p-value`

#We have to do the same for fiadjISI:

fi_missing_match_1_match <- fi_missing_match_1[which(fi_missing_match_1$rsid%in%X_missing_1$snp),]
fi_missing_match_1_match <- fi_missing_match_1_match[order(match(fi_missing_match_1_match$rsid, X_missing_1$snp)),]

length(fi_missing_match_1_match$rsid == X_missing_1$snp) #all of them

X_missing_1$p1 <- fi_missing_match_1_match$`p-value`

#Finally, let's combine it with the missing ones:

X_end <- rbind(X, X_missing_1)

############################################
#Calculating Nuisance and saving dataframes#
############################################

set.seed(100)
varlist <- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(X_end, file = "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IR/fi/MBW/dataframes/X_fi_MBW.rds")
saveRDS(params, file = "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IR/fi/MBW/dataframes/params_fi_MBW.rds")

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

saveRDS(top_vars, "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IR/fi/MBW/dataframes/pruned_fi_MBW.res")

###############
#RUNNING CAUSE#
###############

res <- cause(X=X_end, variants = top_vars, param_ests = params)

saveRDS(res, file = "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IR/fi/MBW/output/res_fi_MBW.rds")

res_strict <- cause(X=X_end, variants = top_vars, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IR/fi/MBW/output/res_fi_MBW_strict.rds")

print(res$elpd)

summary(res, ci_size=0.95)

#And now let's plot these bad fellas:

tiff("/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IR/fi/MBW/output/res_fi_MBW_plot.tiff", units="in", width=1200, height=500, res=300)
plot(res, type="data")
dev.off()
