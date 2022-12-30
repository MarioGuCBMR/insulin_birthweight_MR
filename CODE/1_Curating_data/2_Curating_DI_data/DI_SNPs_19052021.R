##############
#INTRODUCTION#
##############

#This code is to make the DI data for the supplementary tables and for 
#the GRS.

#The code is extracted from the last part of DI_munging.R where I
#already checked that the GRS were done properly.

#Nonetheless, we are gonna do it here again just to organize it
#and make it 100% reproducible.

#The reproducible part is funny because there is a lot of copy paste here
#from supplementary tables and excel sheets from Prokopenko traits.

#Let's go.

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)
library(readxl)

###########
#Functions#
###########

retrieve_beta <- function(beta_se){
  
  beta_ <- as.character(unlist(strsplit(beta_se, "[(]")))[1]
  beta_end <- str_replace_all(string=beta_, pattern=" ", repl="")
  
  return(beta_end)
  
}

retrieve_se <- function(beta_se){
  
  se_ <- as.character(unlist(strsplit(beta_se, "[(]")))[2]
  se_ <- as.character(unlist(strsplit(se_, "[)]")))[1]
  
  se_end <- str_replace_all(string=se_, pattern=" ", repl="")
  
  return(se_end)
  
}

##############
#Loading data#
##############

#Since in Prokopenko et al they didn't focus in the genome-wide significant SNPs for DI,
#we are going to have to take the discovery GWAS for DI and the supplementary table 2A and
#get it done this way.

#1) First we obtain the discovery GWAS from the MAGIC consortium:

DI <- fread("C:/Users/zlc436/Desktop/HERMINA/EMERGENCY_FOLDER/GWAS/Exposure/MAGIC_INSULIN_SECRETION_DI_for_release_HMrel27.txt")

#2) And now we take the supplementary table 2A:

prokopenko_meta_analysis <- read_xlsx("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/CURATED_Work_for_Publication/Rectification_Code/IS_RAW_Data/Prokopenko_MetaAnalysis_SNPs.xlsx")

#And now from the metabochip data from supplementary table 2B:

prokopenko_metabochip <- read.csv("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/CURATED_Work_for_Publication/Rectification_Code/CODE_4_REPLICATION/GRS/INPUT/Metabochip_clean.csv")

prokopenko_metabochip <- prokopenko_metabochip[which(duplicated(prokopenko_metabochip$SNP) == FALSE),]

####################################################################
#Now we add the SNPs from table 1. Let's check where they come from#
####################################################################

#This is done to be able to compare properly between CIR and DI.

test_snps <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/CURATED_Work_for_Publication/Rectification_Code/CODE_4_REPLICATION/GRS/OUTPUT/CLEAN_8_CIR_GRS_SNPs.txt")

##########################################
#Making data for the supplementary tables#
##########################################

#Let's go:

prokopenko_DI <- prokopenko_meta_analysis %>%
  select(SNP, BETA.DI, A1.DI, A2.DI, P.DI)

prokopenko_DI$SE.DI <- as.numeric(as.character(unlist(sapply(prokopenko_DI$BETA.DI, retrieve_se))))
prokopenko_DI$BETA.DI <- as.numeric(as.character(unlist(sapply(prokopenko_DI$BETA.DI, retrieve_beta))))
prokopenko_DI$P.DI <- as.numeric(as.character(prokopenko_DI$P.DI))

prokopenko_DI <- as.data.frame(prokopenko_DI)

new_a1 <- as.character(str_trim(as.character(prokopenko_DI$A1.DI), side = "both"))
new_a2 <- as.character(str_trim(as.character(prokopenko_DI$A2.DI), side = "both"))
new_beta <- as.numeric(as.character(str_trim(as.character(prokopenko_DI$BETA.DI), side = "both")))
new_P <- as.numeric(as.character(str_trim(as.character(prokopenko_DI$P.DI), side = "both")))
new_snp <- as.character(str_trim(as.character(prokopenko_DI$SNP), side = "both"))

prokopenko_DI$A1.DI <- new_a1
prokopenko_DI$A2.DI <- new_a2
prokopenko_DI$BETA.DI <- new_beta
prokopenko_DI$P.DI <- new_P
prokopenko_DI$SNP <- new_snp

#We get only the genome-wide:

prokopenko_DI_gw <- prokopenko_DI[which(prokopenko_DI$P.DI < 0.00000005),] #1 of them

#We also need to add the metabochip data:

prokopenko_metabochip_DI <- prokopenko_metabochip %>%
  select(SNP, BETA.DI, A1, A2, SE.DI, P.DI)

colnames(prokopenko_metabochip_DI) <- c("SNP", "BETA.DI", "A1.DI", "A2.DI", "SE.DI", "P.DI")

new_a1 <- as.character(str_trim(as.character(prokopenko_metabochip_DI$A1.DI), side = "both"))
new_a2 <- as.character(str_trim(as.character(prokopenko_metabochip_DI$A2.DI), side = "both"))
new_beta <- as.numeric(as.character(str_trim(as.character(prokopenko_metabochip_DI$BETA.DI), side = "both")))
new_P <- as.numeric(as.character(str_trim(as.character(prokopenko_metabochip_DI$P.DI), side = "both")))
new_snp <- as.character(str_trim(as.character(prokopenko_metabochip_DI$SNP), side = "both"))

prokopenko_metabochip_DI$A1.DI <- new_a1
prokopenko_metabochip_DI$A2.DI <- new_a2
prokopenko_metabochip_DI$BETA.DI <- new_beta
prokopenko_metabochip_DI$P.DI <- new_P
prokopenko_metabochip_DI$SNP <- new_snp

#########################################################################################
#For the metabochip data we are only going to take care of those that are P< 0.00000005)#
#########################################################################################

prokopenko_metabochip_DI_ind_gw <- prokopenko_metabochip_DI

prokopenko_metabochip_DI_ind_gw <- prokopenko_metabochip_DI_ind_gw[which(prokopenko_metabochip_DI_ind_gw$P.DI < 0.00000005),]

####################################################################
#We also need to add those from Table 1 with the respective values:#
####################################################################

#Now we are gonna merge them:

#To be able to do this, we need to add the SNPs from table 1!
#It is essential to check because we get two SNPs from table 2B that are not there.
#And that might be because they are in the same loci as those in Table 1.

test_snps_2A <- prokopenko_DI[which(prokopenko_DI$SNP%in%test_snps$SNP),] #3 in supl table 2A
test_snps_2B <- prokopenko_metabochip_DI[which(prokopenko_metabochip_DI$SNP%in%test_snps$SNP),] #3 in supl table 2B and 2 are different from supl 2A.
test_snps_discovery <- DI[which(DI$snp%in%test_snps$SNP),] #all in discovery

#We are gonna merge them, but to do that I need...

colnames(test_snps_discovery) <- c("SNP", "A1.DI", "A2.DI", "MAF.DI", "BETA.DI", "SE.DI", "P.DI")
test_snps_2A$MAF.DI <- "-"
test_snps_2B$MAF.DI <- "-"

test_snps_end <- rbind(test_snps_2A, test_snps_2B, test_snps_discovery)

#Now we order by pval and remove duplicates: 

test_snps_end <- test_snps_end[order(test_snps_end$P.DI),]

test_snps_end <- test_snps_end[which(duplicated(test_snps_end$SNP) == FALSE),] #8/8 perfect!

test_snps_end$origin <- "Table 1"

########################
#We merge both datasets#
########################

prokopenko_metabochip_DI_ind_gw$origin <- "Suppl2B"
prokopenko_DI_gw$origin <- "Suppl2A"

prokopenko_DI_end <- rbind(prokopenko_metabochip_DI_ind_gw, prokopenko_DI_gw)

independent_suppl <- prokopenko_DI_end

independent_suppl$rsid <- independent_suppl$SNP
independent_suppl$pval <- independent_suppl$P.DI

independent_suppl_ <- independent_suppl %>%
  select(rsid, pval, origin)

#Now we are going to add the genome-wide significant from the discovery GWAS:

genome_wide_snp <- DI[which(DI$pvalue < 0.00000005),]

#Now the independent ones:

genome_wide_snp$rsid <- genome_wide_snp$snp
genome_wide_snp$pval <- genome_wide_snp$pvalue
genome_wide_snp$origin <- "Discovery"

independent_discovery <- genome_wide_snp %>%
  select(rsid, pval, origin)

#And finally the same for the test_snps from table 1:

independence_table_1 <- test_snps_end %>%
  select(SNP, P.DI, origin)

colnames(independence_table_1) <- c("rsid", "pval", "origin")

independence_test <- rbind(independent_suppl_, independent_discovery, independence_table_1)

#But only those that are p<0.05 for this case...

#We remove duplicates:

independence_test <- independence_test[order(independence_test$pval),]

#Checked, we can remove the duplicates. We will remove the discovery versions or the copies between suppl tables.

independence_test_ <- independence_test[which(duplicated(independence_test$rsid) == FALSE),]
independence_test_ <- independence_test_[which(duplicated(independence_test_$rsid) == FALSE),]
independence_test_ <- independence_test_[which(duplicated(independence_test_$rsid) == FALSE),]

#Just in case...

independence_test_ <- ieugwasr::ld_clump_local(independence_test_, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 1000, clump_r2 = 0.01, clump_p = 1) #we are gonna include them all in this case.

#We got 8 SNPs that pass the filter for this one!!
#The SNPs with the lowest p-value are not exactly the same, 
#but they are in the same loci.

#We have two SNPs that are not the same as CIR, but that is because they are in the same loci:

##############################
#Let's check the effect sizes#
##############################

#If the effect size are similar, then the GRS won't care.

prokopenko_metabochip_DI[which(prokopenko_metabochip_DI$SNP == "rs7756992"),] #0.10
prokopenko_metabochip_DI[which(prokopenko_metabochip_DI$SNP == "rs9368222"),] #0.10

#SAME. It doen't matter which one we use.

#Then, let's go for the other:

prokopenko_metabochip_DI[which(prokopenko_metabochip_DI$SNP == "rs7923866"),] #Not in supplementary tables...
prokopenko_metabochip_DI[which(prokopenko_metabochip_DI$SNP == "rs1111875"),] #0.08

#We are also gonna compare in the originals then...

DI[which(DI$snp == "rs7923866"),] #0.10. Really similar! Let's check in their discovery counterpart. If they are similar too, then we can agree that in the meta-analysis they might also be.
DI[which(DI$snp == "rs1111875"),] #0.09 

#Since what we want is to compare CIR with DI and we only have this difference...
#I think that we are gonna be just fine.

########################################
#CONCLUSIONS: WE ARE GONNA USE CIR DATA#
########################################

#Usually we would use the p<0.05, but we are not gonna because Tuomas said of using them all.

#I leave the code here, just in case.

#test_snps_end <- test_snps_end[which(test_snps_end$P.DI < 0.05),] #this one is not used in the end.

prokopenko_DI_end_ <- test_snps_end

#To continue with the previous code we are going to just do this:

prokopenko_DI <- prokopenko_DI_end_

#Now we use phenoscanner to get the chromosome and the position:

results <- phenoscanner::phenoscanner(prokopenko_DI$SNP)

results <- results$snps

which(results$snp != prokopenko_DI$SNP) #they all match:

prokopenko_DI$CHR.DI <- results$chr
prokopenko_DI$POS.DI <- results$hg19_coordinates

#################################################
#For the MAFs we are gonna use the original data#
#################################################

DI_match <- DI[which(DI$snp%in%prokopenko_DI$SNP),]

DI_match <- DI_match[order(match(DI_match$snp, prokopenko_DI$SNP)),]

#Small check:

which(DI_match$snp != prokopenko_DI$SNP) #perfect.

prokopenko_DI$MAF.DI <- DI_match$maf

#PERFECT.

prokopenko_DI_clean <- prokopenko_DI %>%
  select(SNP, CHR.DI, POS.DI, A1.DI, A2.DI, MAF.DI, BETA.DI, SE.DI, P.DI)

colnames(prokopenko_DI_clean) <- c("SNP", "CHR", "POS", "A1", "A2", "MAF", "BETA", "SE", "P")

#Perfect!!!

#I THINK THAT IS ALL, GODDAMN IT.

manuscript_snps <- prokopenko_DI_clean

#Perfect! Now we order by chromosome:

manuscript_snps <- manuscript_snps[order(manuscript_snps$CHR),]

#Now we switch to the increasing allele:

new_a1 <- ifelse(as.numeric(manuscript_snps$BETA) < 0, manuscript_snps$A2, manuscript_snps$A1)
new_a2 <- ifelse(as.numeric(manuscript_snps$BETA) < 0, manuscript_snps$A1, manuscript_snps$A2)
new_beta <- ifelse(as.numeric(manuscript_snps$BETA) < 0, as.numeric(manuscript_snps$BETA)*(-1), as.numeric(manuscript_snps$BETA))

manuscript_snps$A1 <- new_a1
manuscript_snps$A2 <- new_a2
manuscript_snps$BETA <- new_beta

####################
#NOW IT MAKES SENSE#
####################

write.csv(manuscript_snps, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/CURATED_Work_for_Publication/Rectification_Code/CODE_4_REPLICATION/GRS/OUTPUT/CLEAN_8_DI_GRS_SNPs_Manuscript.csv", quote = FALSE, row.names = FALSE)

write.table(manuscript_snps, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/CURATED_Work_for_Publication/Rectification_Code/CODE_4_REPLICATION/GRS/OUTPUT/CLEAN_8_DI_GRS_SNPs.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
