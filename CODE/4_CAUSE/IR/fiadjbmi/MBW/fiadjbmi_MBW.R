##############
#INTRODUCTION#
##############

#This is code to run CAUSE for fiadjbmi -> MBW

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

fiadjbmi <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/FIadjBMI/FIadjBMI_curated.txt")
MBW <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/OUTPUT/MBW/Birthweight_offspring_mothers2021_curated_FULL.txt")

########################
#Cleaning fiadjbmi data#
########################

#Let's first check the usual:
#The sample size seems weird, let's check whether we have no issue:

summary(fiadjbmi$sample_size) #GREAT. We could erase some of them, but who cares.

#Let's check the fiadjbmi:

summary(fiadjbmi$effect_allele_frequency) #I know that we have NAs, so we are going to still take them and remove them if in the outcome the eaf is weird.

fiadjbmi_NA <- fiadjbmi[which(is.na(fiadjbmi$effect_allele_frequency) == TRUE),]
fiadjbmi_no_NA <- fiadjbmi[which(is.na(fiadjbmi$effect_allele_frequency) == FALSE),]

#Let's run the summary:

summary(fiadjbmi_no_NA$effect_allele_frequency) #let's check these fellas.

fiadjbmi_no_NA <- fiadjbmi_no_NA[which(fiadjbmi_no_NA$effect_allele_frequency > 0.01 & fiadjbmi_no_NA$effect_allele_frequency < 0.99),]

summary(fiadjbmi_no_NA$effect_allele_frequency) #AMAZING

fiadjbmi <- rbind(fiadjbmi_no_NA, fiadjbmi_NA)

#We do not have the rest, so we skip that..., but we do have alleles! CAUSE will only detect
#standard alleles

yes_vect <- c("A", "G", "C", "T")

fiadjbmi <- fiadjbmi[which(fiadjbmi$effect_allele%in%yes_vect),]
fiadjbmi <- fiadjbmi[which(fiadjbmi$other_allele%in%yes_vect),]

#Let's remove the fellas that are not in autosomal data:

yes_chr <- seq(1,22)

fiadjbmi <- fiadjbmi[which(fiadjbmi$chromosome%in%yes_chr),]

#And now we remove the data in the MHC region:

fiadjbmi_mhc <- fiadjbmi[which(as.numeric(fiadjbmi$chromosome) == 6 & as.numeric(fiadjbmi$base_pair_location) > 26000000 & as.numeric(fiadjbmi$base_pair_location) < 34000000),] #so we are gonna approximate it this way.

summary(as.numeric(fiadjbmi_mhc$chromosome)) #perfect
summary(as.numeric(fiadjbmi_mhc$base_pair_location)) #perfect

fiadjbmi_ <- fiadjbmi[which(!(fiadjbmi$chr_pos%in%fiadjbmi_mhc$chr_pos)),]

length(which(as.numeric(fiadjbmi_$chromosome) == 6 & as.numeric(fiadjbmi_$base_pair_location) > 26000000 & as.numeric(fiadjbmi_$base_pair_location) < 34000000)) #we erased them all!!

#Perfect!

rm(fiadjbmi)
rm(fiadjbmi_mhc)

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

#The first merge with the non_missing fiadjbmi_

library(cause)

#Careful, we need to change the columns of fiadjbmi_not_missing: they cannot have "snp".
#It confused CAUSE.

#It didn't work with chr_pos and chr_pos_37, so we are gonna change the columns:

colnames(fiadjbmi_) <- c("chromosome",              "base_pair_location",      "effect_allele",          
                         "other_allele",            "effect_allele_frequency", "beta",                   
                         "standard_error",          "p_value",                 "sample_size",            
                         "chr_pos_37",                 "SNP")

fiadjbmi_$chr_pos_37 <- paste("chr", fiadjbmi_$chr_pos_37, sep = "")

X <- gwas_merge(fiadjbmi_, MBW_, snp_name_cols = c("chr_pos_37", "chr_pos_37"), 
                beta_hat_cols = c("beta", "Beta-A1"), 
                se_cols = c("standard_error", "SE"), 
                A1_cols = c("effect_allele", "A1"), 
                A2_cols = c("other_allele", "A0"))


#Now let's retrieve the rsID for the chr_pos data set.

fiadjbmi_match <- fiadjbmi_[which(fiadjbmi_$chr_pos_37%in%X$snp),]

fiadjbmi_match <- fiadjbmi_match[order(match(fiadjbmi_match$chr_pos_37, X$snp)),]

length(fiadjbmi_match$chr_pos_37 == X$snp) #all of them

#Let's check how are we going to work with this:

check <- fiadjbmi_match[order(fiadjbmi_match$SNP),] #do we have any SNPs without RSID?

head(check$SNP) #yes
tail(check$SNP) #ok

#Let's see how many:

length(which(fiadjbmi_match$SNP == "-")) #not too many!!

fiadjbmi_match_no_rsid <- fiadjbmi_match[which(fiadjbmi_match$SNP == "-"),]
fiadjbmi_match_rsid <- fiadjbmi_match[which(fiadjbmi_match$SNP != "-"),]

length(which(str_detect(fiadjbmi_match_rsid$SNP, "rs") == TRUE)) #not all of them, but most.

weird_ones <- fiadjbmi_match_rsid[which(str_detect(fiadjbmi_match_rsid$SNP, "rs") == FALSE),] #expected ones. Those from Hb1AC.

fiadjbmi_match_rsid <- fiadjbmi_match_rsid[which(str_detect(fiadjbmi_match_rsid$SNP, "rs") == TRUE),] #expected ones. Those from Hb1AC.

#Let's check these SNPs in MBW_, maybe they have RSID there:

fiadjbmi_match_no_rsid <- rbind(fiadjbmi_match_no_rsid, weird_ones)

MBW_match_no_rsid <- MBW_[which(MBW_$chr_pos_37%in%fiadjbmi_match_no_rsid$chr_pos_37),] #WE HAVE THEM.

#Now let's order them and check which ones can we retrieve there.

MBW_match_no_rsid <- MBW_match_no_rsid[order(match(MBW_match_no_rsid$chr_pos_37, fiadjbmi_match_no_rsid$chr_pos_37)),]

length(which(MBW_match_no_rsid$chr_pos_37 == fiadjbmi_match_no_rsid$chr_pos_37)) #DONE.

#Let's be wary and get an index of all of those that have rsids:

index_rsid <- which(str_detect(MBW_match_no_rsid$rsID, "rs") == TRUE) #most of them!!

MBW_match_no_rsid_recovered <- MBW_match_no_rsid[index_rsid,]
fiadjbmi_match_no_rsid_recovered <- fiadjbmi_match_no_rsid[index_rsid,]

length(which(MBW_match_no_rsid_recovered$chr_pos_37 == fiadjbmi_match_no_rsid_recovered$chr_pos_37)) #DONE. 

fiadjbmi_match_no_rsid_recovered$SNP <- MBW_match_no_rsid_recovered$rsID

#AWESOME: now let's be wary with the others:

MBW_match_no_rsid_missed <- MBW_match_no_rsid[-index_rsid,]
fiadjbmi_match_no_rsid_missed <- fiadjbmi_match_no_rsid[-index_rsid,]

#Let's try with phenoscanner:

ps_results_1 <- phenoscanner::phenoscanner(MBW_match_no_rsid_missed$chr_pos_37[seq(1,100)])
ps_results_2 <- phenoscanner::phenoscanner(MBW_match_no_rsid_missed$chr_pos_37[seq(101,119)])

ps_end <- rbind(ps_results_1$snps, ps_results_2$snps) #83/120! not bad.

#We can use these fellas.
#As for the other SNPs, I think it is OK if we miss 40 SNPs, honestly.
#We did as great as we could!!

MBW_match_no_rsid_missed_recovered <- MBW_match_no_rsid_missed[which(MBW_match_no_rsid_missed$chr_pos_37%in%ps_end$hg19_coordinates),]
fiadjbmi_match_no_rsid_missed_recovered <- fiadjbmi_match_no_rsid_missed[which(fiadjbmi_match_no_rsid_missed$chr_pos_37%in%ps_end$hg19_coordinates),]

ps_end <- ps_end[order(match(ps_end$hg19_coordinates, MBW_match_no_rsid_missed_recovered$chr_pos_37)),]

length(which(ps_end$hg19_coordinates == MBW_match_no_rsid_missed_recovered$chr_pos_37))
length(which(ps_end$hg19_coordinates == fiadjbmi_match_no_rsid_missed_recovered$chr_pos_37))

fiadjbmi_match_no_rsid_missed_recovered$SNP <- ps_end$rsid

#Finally, let's get those that are not there and put a "-"

MBW_match_no_rsid_missed_missed <- MBW_match_no_rsid_missed[which(!(MBW_match_no_rsid_missed$chr_pos_37%in%ps_end$hg19_coordinates)),]
fiadjbmi_match_no_rsid_missed_missed <- fiadjbmi_match_no_rsid_missed[which(!(fiadjbmi_match_no_rsid_missed$chr_pos_37%in%ps_end$hg19_coordinates)),]

#Let's combine fiadjbmi_match!!

fiadjbmi_match_end <- rbind(fiadjbmi_match_no_rsid_missed_missed, fiadjbmi_match_no_rsid_recovered, fiadjbmi_match_no_rsid_missed_recovered, fiadjbmi_match_rsid)

#And finally we can do this:

fiadjbmi_match_end <- fiadjbmi_match_end[order(match(fiadjbmi_match_end$chr_pos_37, X$snp)),]

length(fiadjbmi_match_end$chr_pos_37 == X$snp) #perfect

#We can do this!

X$snp <- fiadjbmi_match_end$SNP
X$p1 <- fiadjbmi_match_end$p_value

X_end <- X #just to keep with all

############################################
#Calculating Nuisance and saving dataframes#
############################################

set.seed(100)
varlist <- with(X_end, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X_end, varlist)

saveRDS(X_end, file = "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IR/fiadjbmi/MBW/dataframes/X_fiadjbmi_MBW.rds")
saveRDS(params, file = "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IR/fiadjbmi/MBW/dataframes/params_fiadjbmi_MBW.rds")

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

saveRDS(top_vars, "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IR/fiadjbmi/MBW/dataframes/pruned_fiadjbmi_MBW.res")

###############
#RUNNING CAUSE#
###############

res <- cause(X=X_end, variants = top_vars, param_ests = params)

saveRDS(res, file = "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IR/fiadjbmi/MBW/output/res_fiadjbmi_MBW.rds")

res_strict <- cause(X=X_end, variants = top_vars, param_ests = params, qalpha = 1, qbeta = 2)

saveRDS(res, file = "/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IR/fiadjbmi/MBW/output/res_fiadjbmi_MBW_strict.rds")

print(res$elpd)

summary(res, ci_size=0.95)

#And now let's plot these bad fellas:

tiff("/home/projects/ku_00095/people/marure/BW_project/CODE/4_CAUSE/IR/fiadjbmi/MBW/output/res_fiadjbmi_MBW_plot.tiff", units="in", width=1200, height=500, res=300)
plot(res, type="data")
dev.off()
