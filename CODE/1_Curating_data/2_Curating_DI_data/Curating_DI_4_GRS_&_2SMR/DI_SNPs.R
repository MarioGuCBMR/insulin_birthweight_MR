##############
#INTRODUCTION#
##############

#For 2SMR we are going to use the 8 DI variants reported in Table 1 from Prokopenko study.

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

DI <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Glycemic_traits_raw_data/MAGIC_INSULIN_SECRETION_DI_for_release_HMrel27.txt")

DI_match <- DI[which(DI$snp%in%table_1_meta$SNP),] 

summary(DI_match$pvalue) #One of the variants is not P < 0.05. We remove it. 

fail_snp <- DI_match$snp[which(DI_match$pvalue > 0.05)]

table_1_meta <- table_1_meta[-which(table_1_meta$SNP == fail_snp),]

DI_match <- DI[which(DI$snp%in%table_1_meta$SNP),] 

#########################################################################
#Let's check if the SNPs meet the clumping requirements for GRS and 2SMR#
#########################################################################

DI_match$rsid <- DI_match$snp
DI_match$pval <- DI_match$pvalue

independence_test_end <- ieugwasr::ld_clump_local(DI_match, bfile = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/2SMR/1000G/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.1/plinkbinr/bin/plink_Windows.exe", clump_kb = 500, clump_r2 = 0.01, clump_p = 0.05) #the minimum needs to be p = 0.05 for GRS.

#########################
#We have 7 SNPs in total# 
#########################

#We need chr_pos:

DI_ps <- phenoscanner::phenoscanner(independence_test_end$snp)$snp

#Let's check that the SNPs are the same:

length(which(DI_ps$snp == independence_test_end$snp)) #perfect match

#PERFECT

independence_test_end$CHR <- DI_ps$chr
independence_test_end$POS <- DI_ps$pos_hg19

independence_test_end <- independence_test_end[order(match(independence_test_end$snp, table_1_meta$SNP)),]

length(which(independence_test_end$snp == table_1_meta$SNP)) #OK
 
#For the allele frequency, I am gonna work with the OAF of the table_1.

independence_test_end$EAF <- ifelse(independence_test_end$effect_allele != table_1_meta$A2, 1-as.numeric(table_1_meta$OAF), as.numeric(table_1_meta$OAF))

#PERFECT

#Let's align the data to the DI-increasing allele

new_a1 <- ifelse(as.numeric(independence_test_end$effect) < 0, as.character(independence_test_end$other_allele), as.character(independence_test_end$effect_allele))
new_a2 <- ifelse(as.numeric(independence_test_end$effect) < 0, as.character(independence_test_end$effect_allele), as.character(independence_test_end$other_allele))
new_beta <- ifelse(as.numeric(independence_test_end$effect) < 0, as.numeric(independence_test_end$effect)*(-1), as.numeric(independence_test_end$effect))
new_eaf <- ifelse(as.numeric(independence_test_end$effect) < 0, as.numeric(1-independence_test_end$EAF), as.numeric(independence_test_end$EAF))

independence_test_end$Final_A1 <- new_a1
independence_test_end$Final_A2 <- new_a2
independence_test_end$Final_BETA <- new_beta
independence_test_end$Final_EAF <- new_eaf

#Let's check if we did this properly:

head(independence_test_end) #worked like a freaking work

manuscript_snps <- independence_test_end %>%
  select(snp, CHR, POS, Final_A1, Final_A2, Final_EAF, Final_BETA, stderr, pvalue)

write.csv(manuscript_snps, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/DI/CLEAN_7_DI_GRS_SNPs_Manuscript.csv", quote = FALSE, row.names = FALSE)

write.table(manuscript_snps, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/DI/CLEAN_7_DI_GRS_SNPs_Manuscript.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)
