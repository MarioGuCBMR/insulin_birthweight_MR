##############
#INTRODUCTION#
##############

#For 2SMR we are going to use the 8 CIR variants reported in Table 1 from Prokopenko study.

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)
library(readxl)

##############
#Loading data#
##############

table_1_meta <- readxl::read_xlsx("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Glycemic_traits_raw_data/CIR_GW_SNPs.xlsx")

table_1_meta$EAF <- 1-as.numeric(table_1_meta$OAF)
table_1_meta$BETA <- as.numeric(table_1_meta$BETA)
table_1_meta$SE <- as.numeric(table_1_meta$SE)
table_1_meta$P <- as.numeric(table_1_meta$P)

#########################################################################
#Let's check if the SNPs meet the clumping requirements for GRS and 2SMR#
#########################################################################

table_1_meta <- as.data.frame(table_1_meta)

table_1_meta$rsid <- table_1_meta$SNP
table_1_meta$pval <- table_1_meta$P

independence_test_end <- ieugwasr::ld_clump_local(table_1_meta, bfile = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/2SMR/1000G/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.1/plinkbinr/bin/plink_Windows.exe", clump_kb = 500, clump_r2 = 0.01, clump_p = 0.05) #the minimum needs to be p = 0.05 for GRS.

#########################
#We have 8 SNPs in total# 
#########################

#The main issue is..., that we do not have EAF, so we are going to load the RAW_DATA.

CIR <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Glycemic_traits_raw_data/MAGIC_INSULIN_SECRETION_CIR_for_release_HMrel27.txt")

CIR_match <- CIR[which(CIR$snp%in%independence_test_end$rsid),]

head(CIR_match)

#We need chr_pos:

CIR_ps <- phenoscanner::phenoscanner(CIR_match$snp)$snp

#Let's check that the SNPs are the same:

length(which(CIR_ps$snp == CIR_match$snp)) #perfect match

#PERFECT

CIR_match$CHR <- CIR_ps$chr
CIR_match$POS <- CIR_ps$pos_hg19

#Let's check if this makese sense:

CIR_match <- CIR_match[order(match(CIR_match$snp, independence_test_end$rsid)),]

#Let's check:

length(which(CIR_match$snp == independence_test_end$rsid))

independence_test_end$chr <- CIR_match$CHR
independence_test_end$pos <- CIR_match$POS

#NOW IT MAKES SENSE.
#I checked with previous data and it all makes sense.
#This info is what we are gonna put in the manuscript.

manuscript_snps <- independence_test_end %>%
  select(SNP, chr, pos, A1, A2, EAF, BETA, SE, P)

write.csv(manuscript_snps, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/CIR/CLEAN_8_CIR_GRS_SNPs_Manuscript.csv", quote = FALSE, row.names = FALSE)

write.table(manuscript_snps, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/CIR/CLEAN_8_CIR_GRS_SNPs_Manuscript.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)