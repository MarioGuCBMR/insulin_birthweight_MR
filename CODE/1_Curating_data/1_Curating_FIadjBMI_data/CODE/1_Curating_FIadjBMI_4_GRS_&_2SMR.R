##############
#INTRODUCTION#
##############

#This code will serve to parse the data of FIadjBMI used in for the GRS and 2SMR.

###################
#LOADING FUNCTIONS#
###################

library(tidyverse)
library(xlsx)
library(data.table)

################
#LOAD FUNCTIONS#
################

load_independent_snps <- function(path_){
  #This function takes a path and the initial part of the name of the chromosome data (the last thing that needs to be written is chr.)
  #And loads all the data in the folder, chromosome by chromosome and adds the SNPs in a vector.
  
  #path_ <- path_fi_combined_lipid_glgc
  
  independent_snp_lax <- c()
  
  for(i in seq(1,22)){
    
    path_end <- paste(path_, i, ".clumped", sep = "")
    
    tmp_ <- fread(path_end)
    
    independent_snp_lax <- c(independent_snp_lax, tmp_$SNP)
    
  }
  
  return(independent_snp_lax)
  
}

##############
#LOADING DATA#
##############

full_raw_data <- readxl::read_xlsx("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Glycemic_traits_raw_data/MAGIC_Glycemic_Traits_Adjusted_BMI_Supplementary_Tables.xlsx", sheet = 7)

#################################################
#STEP 0: PARSING DATA SO THAT IT CAN BE READABLE#
#################################################

#Let's read the columns properly:

colnames(full_raw_data) <- as.character(full_raw_data[1,])

View(full_raw_data)

#Finally we clean the weird rows:

full_data <- full_raw_data[-1,]


#####################################
#STEP 1: GET the lead SNPs per trait#
#####################################

FI_SNPs <- full_data[which(full_data$Trait == "FI"),] #95 SNPs -> the same as we identified before.

#Now, which ones are used for GRS???

FI_SNPs_GRS <- FI_SNPs[which(FI_SNPs$`List 4: EUR GS (not heterogeneous)` == 1),] #36

FI_SNPs_GRS$chr_pos <- paste(FI_SNPs_GRS$Chr, FI_SNPs_GRS$`Pos (bp)`, sep = ":")

##################################################################
#STEP 3: retrieve the variants in the original data and save them#
##################################################################

#1. let's get the summary statistics data for each of the variants:

FIAdjBMI_original <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Romain_Mario/RAW_DATA/FIAdjBMI_2021/MAGIC1000G_FI_EUR.tsv.gz")

FIAdjBMI_original$chr_pos <- paste(FIAdjBMI_original$chromosome, FIAdjBMI_original$base_pair_location, sep = ":")

FIAdjBMI_match <- FIAdjBMI_original[which(FIAdjBMI_original$chr_pos%in%FI_SNPs_GRS$chr_pos),]

FIAdjBMI_match <- FIAdjBMI_match[order(match(FIAdjBMI_match$chr_pos, FI_SNPs_GRS$chr_pos)),]

length(which(FIAdjBMI_match$chr_pos !=FI_SNPs_GRS$chr_pos))
length(which(FIAdjBMI_match$chr_pos ==FI_SNPs_GRS$chr_pos))

FIAdjBMI_match$rsid <- FI_SNPs_GRS$`Lead/index variant (rsID)`

FIAdjBMI_match$MAF <- ifelse(as.numeric(FIAdjBMI_match$effect_allele_frequency) > 0.5, 1-as.numeric(FIAdjBMI_match$effect_allele_frequency), as.numeric(FIAdjBMI_match$effect_allele_frequency))

#2. Let's get the variables that we need:

FIAdjBMI_match_4_SM <- FIAdjBMI_match %>%
  select(rsid, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, MAF, beta, standard_error, p_value, sample_size)

colnames(FIAdjBMI_match_4_SM) <- c("rsID", "Chromosome", "Position_build_37", "Effect_Allele", "Other_Allele", "EAF", "MAF", "Beta", "SE", "P", "N")

#3. Let's change to numeric so that we do not have surprises later:

FIAdjBMI_match_4_SM$Beta <- as.numeric(FIAdjBMI_match_4_SM$Beta)
FIAdjBMI_match_4_SM$EAF <- as.numeric(FIAdjBMI_match_4_SM$EAF)
FIAdjBMI_match_4_SM$SE <- as.numeric(FIAdjBMI_match_4_SM$SE)
FIAdjBMI_match_4_SM$P <- as.numeric(FIAdjBMI_match_4_SM$P)

#This first dataframe is going to go for the analysis.

fwrite(as.data.frame(FIAdjBMI_match_4_SM), "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/FIadjBMI_independent.txt", col.names = TRUE, row.names = FALSE, quote = TRUE)

#And now we are going to do the same, but making a beautiful supplementary table.

write.xlsx(as.data.frame(FIAdjBMI_match_4_SM), file="J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/Supplementary_Tables_V1.xlsx", sheetName="Supplementary Table 1", row.names=FALSE)

###################################################
#STEP 4: CLUMP DATA WITH THE SETTINGS THAT WE NEED#
###################################################

FIAdjBMI_match_4_SM$rsid <- FIAdjBMI_match_4_SM$rsID
FIAdjBMI_match_4_SM$pval <- FIAdjBMI_match_4_SM$P

independence_test_end <- ieugwasr::ld_clump_local(FIAdjBMI_match_4_SM, bfile = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/2SMR/1000G/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/R-4.1.2/library/plinkbinr/bin/plink_Windows.exe", clump_kb = 500, clump_r2 = 0.01, clump_p = 0.05) 

#We need to check which ones are being removed.

#rs111264094: the allele frequency is really close to 0.99. You could say it is 0.99. We might want to keep it and check in LDPair. 
#rs77935490: it is triallelic. We will want to check in LDPair, for example.
#rs200678953: it is an insertion. This will be removed in the end.

#Let's check if any variant is removed due to clumping or to absences.

no_vect <- c("rs111264094", "rs77935490", "rs200678953")

independence_test_end_2 <- FIAdjBMI_match_4_SM[-which(FIAdjBMI_match_4_SM$rsID%in%no_vect),]

independence_test_end_2 <- ieugwasr::ld_clump_local(independence_test_end_2, bfile = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/2SMR/1000G/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/R-4.1.2/library/plinkbinr/bin/plink_Windows.exe", clump_kb = 500, clump_r2 = 0.01, clump_p = 0.05) 

#All the rest are present in the reference panel!!

#Let's take a close look to the others and see if the can be included.

######################
#CHECKING rs111264094#
######################

#First we are going to check the SNPs in order and see if he has any neighbours.

FIAdjBMI_match_4_SM <- FIAdjBMI_match_4_SM[order(FIAdjBMI_match_4_SM$Chromosome, FIAdjBMI_match_4_SM$Position_build_37),]

View(FIAdjBMI_match_4_SM)

21699928-48202696
66351826-48202696

#There is no issue of LD here. 
#They are all so far away there is no need to do any check.
#We will chack the EAF on the outcome every time. just in case for this one.

#####################
#CHECKING rs77935490#
#####################

#The same: it is too far away to do anything!
#We can also keep it. 

#In conclusion, we are going to add them, but doing checks for each outcome.

independence_test_end <- independence_test_end %>%
  select(-c(rsid, pval))

FIAdjBMI_match_4_SM <- FIAdjBMI_match_4_SM %>%
  select(-c(rsid, pval))

pass_filter <- c("rs77935490", "rs111264094")

included_snps <- FIAdjBMI_match_4_SM[which(FIAdjBMI_match_4_SM$rsID%in%pass_filter),]

final_set_of_indepedendent <- rbind(independence_test_end, included_snps)

fwrite(as.data.frame(final_set_of_indepedendent), "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/FIadjBMI_independent_4_GRS.csv", col.names = TRUE, row.names = FALSE, quote = TRUE)
