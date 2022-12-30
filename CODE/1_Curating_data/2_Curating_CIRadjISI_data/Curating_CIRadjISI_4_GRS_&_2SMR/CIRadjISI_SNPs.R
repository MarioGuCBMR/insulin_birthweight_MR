##############
#INTRODUCTION#
##############

#For 2SMR we are going to use the 8 CIRadjISI variants reported in Table 1 from Prokopenko study.

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)
library(readxl)

##############
#Loading data#
##############

table_1_meta <- readxl::read_xlsx("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Glycemic_traits_raw_data/CIRadjISI_GW_SNPs.xlsx")

CIRadjISI <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Glycemic_traits_raw_data/MAGIC_INSULIN_SECRETION_CIRadjISI_ISI_for_release_HMrel27.txt")

CIRadjISI_match <- CIRadjISI[which(CIRadjISI$snp%in%table_1_meta$SNP),] 

summary(CIRadjISI_match$pvalue) #all of them are significant! We can use them.

#########################################################################
#Let's check if the SNPs meet the clumping requirements for GRS and 2SMR#
#########################################################################

CIRadjISI_match$rsid <- CIRadjISI_match$snp
CIRadjISI_match$pval <- CIRadjISI_match$pvalue

independence_test_end <- ieugwasr::ld_clump_local(CIRadjISI_match, bfile = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/2SMR/1000G/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.1/plinkbinr/bin/plink_Windows.exe", clump_kb = 500, clump_r2 = 0.01, clump_p = 0.05) #the minimum needs to be p = 0.05 for GRS.

#########################
#We have 8 SNPs in total# 
#########################


#We need chr_pos:

CIRadjISI_ps <- phenoscanner::phenoscanner(independence_test_end$snp)$snp

#Let's check that the SNPs are the same:

length(which(CIRadjISI_ps$snp == independence_test_end$snp)) #perfect match

#PERFECT

independence_test_end$CHR <- CIRadjISI_ps$chr
independence_test_end$POS <- CIRadjISI_ps$pos_hg19

independence_test_end <- independence_test_end[order(match(independence_test_end$snp, table_1_meta$SNP)),]

length(which(independence_test_end$snp == table_1_meta$SNP)) #OK
 
#For the allele frequency, I am gonna work with the OAF of the table_1.

independence_test_end$EAF <- ifelse(independence_test_end$effect_allele != table_1_meta$A2, 1-as.numeric(table_1_meta$OAF), as.numeric(table_1_meta$OAF))

#PERFECT

#Let's align the data to the CIRadjISI-increasing allele

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

write.csv(manuscript_snps, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/CIRadjISI/CLEAN_8_CIRadjISI_GRS_SNPs_Manuscript.csv", quote = FALSE, row.names = FALSE)

write.table(manuscript_snps, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/CIRadjISI/CLEAN_8_CIRadjISI_GRS_SNPs_Manuscript.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)
