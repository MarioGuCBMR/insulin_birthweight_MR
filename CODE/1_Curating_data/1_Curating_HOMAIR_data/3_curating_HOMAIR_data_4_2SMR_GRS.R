##############
#INTRODUCTION#
##############

#This data is to obtain the 2SMR and GRS variants for HOMA-IR.
#We are just going to go for: 

#What if we try use the two genome-wide significant HOMA-IR variants from Table 1 reported by Dupuis et al., based on the combined meta-analysis of discovery + replication stages:
#- GCKR rs780094
#- IGF1 rs35767

#These are from meta-analysis + replication stages.
#We might need to put the data manually.

#We have no information of SE with the replication stages.
#We are going to use the data of the original SS for these SNPs. 
#I think that is the best way.

#########
#library#
#########

library(tidyverse)
library(data.table)

###################
#Loading functions#
###################

chr_parser <- function(chr_pos){
  #A function that returns the chromosome from a chr:pos string
  
  tmp <- strsplit(chr_pos, ":")[[1]][1]
  chr_ <- strsplit(tmp, "chr")[[1]][2]
  
  return(chr_)
  
  
}

pos_parser <- function(chr_pos){
  #A function that returns the chromosome from a chr:pos string
  
  pos_ <- strsplit(chr_pos, ":")[[1]][2]

  return(pos_)
  
  
}

##############
#Loading data#
##############

homair <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/HOMA_IR/HOMA_IR_combined_Curated.txt")

############################################################################
#First we are going to take the SNPs with their original summary statistics#
############################################################################

snps_of_interest <- c("rs780094", "rs35767")

homair_gw <- homair[which(homair$snp%in%snps_of_interest),] 

########################################################################
#Now we are going to filter so that the data that we have is consistent#
########################################################################

homair_gw$N <- 37037

homair_gw$chr <- as.numeric(as.character(unlist(sapply(homair_gw$chr_pos_37, chr_parser))))
homair_gw$pos <- as.numeric(as.character(unlist(sapply(homair_gw$chr_pos_37, pos_parser))))
homair_gw$effect_allele <- as.character(unlist(sapply(homair_gw$effect_allele, toupper)))
homair_gw$other_allele <- as.character(unlist(sapply(homair_gw$other_allele, toupper)))

#Since we have very few SNPs, we can check the EAF and see if it matches with the MAF.

head(homair_gw) #for the first one, G has an EAf of 0.84~, So MAF=EAF in this case (A1=A and EAF = 0.115). For the second SNP T has an EAF of 0.40. T is our A1 and the MAF is 0.40... so all good.

homair_curated <- homair_gw %>%
  select(snp, chr, pos, effect_allele, other_allele, maf, effect, stderr, pvalue, N)

colnames(homair_curated) <- c("rsID", "Chromosome", "Position_build_37", "Effect_Allele", "Other_Allele", "EAF", "Beta", "SE", "P", "N") #we can switch MAF for EAF :)

#####################
#Let's save the data#
#####################

write.csv(homair_curated, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/FI/CLEAN_2_HOMAIR_GRS_SNPs_Manuscript.csv", quote = FALSE, row.names = FALSE)
write.table(homair_curated, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/FI/CLEAN_2_HOMAIR_SNPs_Manuscript.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)


