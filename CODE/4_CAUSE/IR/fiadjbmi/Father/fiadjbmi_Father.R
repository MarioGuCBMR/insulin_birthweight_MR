##############
#INTRODUCTION#
##############

#This is code to run CAUSE for fiadjbmi -> Father

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

fiadjbmi <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/fiadjbmi/fiadjbmi_curated.txt")
Father <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/father/Birthweight_offspring_fathers2021_curated_FULL.txt")

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

X <- gwas_merge(fiadjbmi_, Father_, snp_name_cols = c("chr_pos_37", "chr_pos_37"), 
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

#Let's check these SNPs in Father_, maybe they have RSID there:

fiadjbmi_match_no_rsid <- rbind(fiadjbmi_match_no_rsid, weird_ones)

Father_match_no_rsid <- Father_[which(Father_$chr_pos_37%in%fiadjbmi_match_no_rsid$chr_pos_37),] #WE HAVE THEM.

#Now let's order them and check which ones can we retrieve there.

Father_match_no_rsid <- Father_match_no_rsid[order(match(Father_match_no_rsid$chr_pos_37, fiadjbmi_match_no_rsid$chr_pos_37)),]

length(which(Father_match_no_rsid$chr_pos_37 == fiadjbmi_match_no_rsid$chr_pos_37)) #DONE.

#Let's be wary and get an index of all of those that have rsids:

index_rsid <- which(str_detect(Father_match_no_rsid$rsID, "rs") == TRUE) #most of them!!

Father_match_no_rsid_recovered <- Father_match_no_rsid[index_rsid,]
fiadjbmi_match_no_rsid_recovered <- fiadjbmi_match_no_rsid[index_rsid,]

length(which(Father_match_no_rsid_recovered$chr_pos_37 == fiadjbmi_match_no_rsid_recovered$chr_pos_37)) #DONE. 

fiadjbmi_match_no_rsid_recovered$SNP <- Father_match_no_rsid_recovered$rsID

#AWESOME: now let's be wary with the others:

Father_match_no_rsid_missed <- Father_match_no_rsid[-index_rsid,]
fiadjbmi_match_no_rsid_missed <- fiadjbmi_match_no_rsid[-index_rsid,]

#Let's try with phenoscanner:

ps_results_1 <- phenoscanner::phenoscanner(Father_match_no_rsid_missed$chr_pos_37[seq(1,100)])
ps_results_2 <- phenoscanner::phenoscanner(Father_match_no_rsid_missed$chr_pos_37[seq(101,109)])

ps_end <- rbind(ps_results_1$snps, ps_results_2$snps) #83/120! not bad.

#We can use these fellas.
#As for the other SNPs, I think it is OK if we miss 40 SNPs, honestly.
#We did as great as we could!!

Father_match_no_rsid_missed_recovered <- Father_match_no_rsid_missed[which(Father_match_no_rsid_missed$chr_pos_37%in%ps_end$hg19_coordinates),]
fiadjbmi_match_no_rsid_missed_recovered <- fiadjbmi_match_no_rsid_missed[which(fiadjbmi_match_no_rsid_missed$chr_pos_37%in%ps_end$hg19_coordinates),]

ps_end <- ps_end[order(match(ps_end$hg19_coordinates, Father_match_no_rsid_missed_recovered$chr_pos_37)),]

length(which(ps_end$hg19_coordinates == Father_match_no_rsid_missed_recovered$chr_pos_37))
length(which(ps_end$hg19_coordinates == fiadjbmi_match_no_rsid_missed_recovered$chr_pos_37))

fiadjbmi_match_no_rsid_missed_recovered$SNP <- ps_end$rsid

#Finally, let's get those that are not there and put a "-"

Father_match_no_rsid_missed_missed <- Father_match_no_rsid_missed[which(!(Father_match_no_rsid_missed$chr_pos_37%in%ps_end$hg19_coordinates)),]
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
