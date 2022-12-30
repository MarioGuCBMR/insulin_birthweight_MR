##############
#INTRODUCTION#
##############

#This is a quick code to check how much do we need to clean the CIRadjISI stratified data.

###########
#Libraries#
###########

library(data.table)
library(tidyverse)
library(TwoSampleMR)

##################
#Loading function#
##################

clean_liftover <- function(chr_pos){
  
  tmp_1 <- strsplit(chr_pos, ":")[[1]][1]
  tmp_2 <- strsplit(chr_pos, ":")[[1]][2]
  tmp_3 <- as.numeric(as.character(strsplit(tmp_2, "-")[[1]][1])) +1
  
  final <- paste(tmp_1, tmp_3, sep = ":")
  
  return(final)
  
}

parse_rsid <- function(rsid){
  
  rsid_ <- strsplit(rsid, ":")[[1]][1] #this works for those that do not present the weird formatting :)
  
  return(rsid_)
  
}

##############
#Loading data#
##############

CIRadjISI <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Glycemic_traits_raw_data/MAGIC_INSULIN_SECRETION_CIR_ISI_for_release_HMrel27.txt") 

#Let's check what columns do we have:

CIRadjISI <- CIRadjISI[order(CIRadjISI$snp),]

head(CIRadjISI) #perfect
tail(CIRadjISI) #perfect

##############################
#Filtering data and saving it#
##############################

#What we are going to do is the following:

#a) Filter for MAF/INFO if possible.
#b) remove MHC region SNPs here, if possible. If not, in the outcome.
#c) Save the data so we do not have to wait for the long query in the database.

#The CIRadjISI data has some of the SNPs without allele frequency.
#We cannot know if the imputation was done correctly...
#Hence, we are gonna remove them.

summary(CIRadjISI$maf) #perfect.

CIRadjISI_maf <- CIRadjISI[which(CIRadjISI$maf > 0.01),]

summary(CIRadjISI_maf$maf) #perfect.

#Now the INFO, do we have that data?

colnames(CIRadjISI_maf) #we do not. 

################################################
#Getting as many chr_pos possible from BMI data#
################################################

bmi <- fread("C:/Users/zlc436/Desktop/PULIT_data_reservoir/Combined/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")

#We need to be careful, some SNPs from Pulit et al do not present correct chr_pos.
#The mistake is in those without INFO data, thus, the ones in GIANT data.
#We are gonna remove those.

summary(bmi$INFO) 
summary(bmi$Freq_Tested_Allele) 
summary(bmi$CHR) 
summary(bmi$N) 

bmi <- bmi[which(is.na(bmi$INFO) == FALSE),]

summary(bmi$INFO) #all good
summary(bmi$Freq_Tested_Allele) #all good
summary(bmi$CHR) #all good
summary(bmi$N) #all good

#It seems all is OK.
#Now we need to clean the RSIDs, because the have the effect allele and the other allele.
#Which is actually fine.

bmi$RSID <- as.character(unlist(sapply(bmi$SNP, parse_rsid)))

#Now we are gonna match by SNP, 
#but we might miss some of them because they were merged.
#But we can do that later.

bmi_match <- bmi[which(bmi$RSID%in%CIRadjISI_maf$snp),] #most of them:

CIRadjISI_Pulit <- CIRadjISI_maf[which(CIRadjISI_maf$snp%in%bmi_match$RSID),]

#We have duplicates, but they are bound to have the same chr_pos, so we do not care
#about them. Let's check just in case:

dupl <- bmi_match$RSID[which(duplicated(bmi_match$RSID) == TRUE)]

bmi_dupl <- bmi_match[which(bmi_match$RSID%in%dupl),]

#Let's order them and get ready the dataframe:

bmi_dupl <- bmi_dupl[order(bmi_dupl$RSID),]

head(bmi_dupl) #same chr_pos
tail(bmi_dupl) #same chr_pos

#PERFECT. We can remove duplicates:

bmi_match <- bmi_match[which(duplicated(bmi_match$RSID) == FALSE),] #now they match.

bmi_match <- bmi_match[order(match(bmi_match$RSID, CIRadjISI_Pulit$snp)),]

#Let's check if we did it properly:

length(which(bmi_match$RSID != CIRadjISI_Pulit$snp)) #perfect.
length(which(bmi_match$RSID == CIRadjISI_Pulit$snp)) #perfect.

#We can match the chr_pos:

CIRadjISI_Pulit$chr_pos_37 <- paste("chr", bmi_match$CHR, ":", bmi_match$POS, sep = "")

#########
#PERFECT#
#########

#Let's go and get those missing:

CIRadjISI_missing <- CIRadjISI_maf[which(!(CIRadjISI_maf$snp%in%bmi_match$RSID)),]

#Let's see how many of them can be recovered from SNP Nexus:

for_nexus <- CIRadjISI_missing

for_nexus$dbsnp <- "dbsnp"

for_nexus <- for_nexus %>%
  dplyr::select(dbsnp, snp)

for_nexus <- for_nexus[order(for_nexus$snp),]

head(for_nexus)
tail(for_nexus)

write.table(for_nexus, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/SNPNexus/INPUT/CIRadjISI/CIRadjISI_combined_4_Nexus.txt", sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

####################
#Let's recover them#
####################

recovered_snps <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/SNPNexus/OUTPUT/CIRadjISI/CIRadjISI_combined_in_Nexus.txt")

recovered_snps <- recovered_snps[order(recovered_snps$`REF Allele`),]

head(recovered_snps) #nothing.
tail(recovered_snps) #nothing.

#It seems we only need to remove the deletion. 
#To do so, since we have IUPAC alleles we will just take all the alleles
#that are letters, to avoid insertions and deletions.

good_alleles <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K",
                  "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U",
                  "V", "W", "X", "Y", "Z")

recovered_snps <- recovered_snps[which(recovered_snps$`REF Allele`%in%good_alleles),] 
recovered_snps <- recovered_snps[which(recovered_snps$`ALT Allele (IUPAC)`%in%good_alleles),] 

CIRadjISI_recovered <- CIRadjISI_missing[which(CIRadjISI_missing$snp%in%recovered_snps$dbSNP),]

recovered_snps <- recovered_snps[order(match(recovered_snps$dbSNP, CIRadjISI_recovered$snp)),]

#Doing some checkity checks:

which(recovered_snps$dbSNP != CIRadjISI_recovered$MarkerName)

head(recovered_snps$dbSNP)
head(CIRadjISI_recovered$snp)

#PERFECT:

chr_ <- paste("chr", recovered_snps$Chromosome, sep = "")

recovered_snps$chr_pos <- paste(chr_, recovered_snps$Position, sep = ":")

CIRadjISI_recovered$chr_pos_37 <- recovered_snps$chr_pos

###############################
#Let's check the missing ones:#
###############################

CIRadjISI_missing_2 <- CIRadjISI_missing[which(!(CIRadjISI_missing$snp%in%recovered_snps$dbSNP)),]

check <- myvariant::queryVariants(CIRadjISI_missing_2$snp)

#Clean: 

check <- check[which(duplicated(check$query) == FALSE),]

check$chr_pos <- paste("chr", check$chrom, ":", check$hg19.end, sep = "")

check_df <- as.data.frame(check)

check_df <- check_df[order(match(check_df$query, CIRadjISI_missing_2$snp)),]

length(which(check_df$query == CIRadjISI_missing_2$snp)) #all of them

CIRadjISI_missing_2$chr_pos_37 <- check_df$chr_pos #349 still missing.

#But still, this is the best we are going to get.
#We will have to be wary with the matching, but still.

###############################
#Let's mix everything and save#
###############################

CIRadjISI_end <- rbind(CIRadjISI_Pulit, CIRadjISI_recovered, CIRadjISI_missing_2) #PERFECT number. It matches the data perfectily.

dim(CIRadjISI_end)
dim(CIRadjISI_maf)

fwrite(CIRadjISI_end, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/CIRadjISI/CIRadjISI_combined_Curated.txt")
