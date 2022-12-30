##############
#INTRODUCTION#
##############

#This code is match the chr_pos that we have from the data that has to be clumped with the bim 
#data for all the chromosome in the 1000G data Phase 3 V5. 
#The clumping is not working because the IDs are in a weird format that combines RSID, CHR_POS, and alleles.

#I cannot generate that because the reference and the alternative alleles are swapped!!
#However, this format allows us to keep the triallelic SNPs and clump them, which is
#gonna be cool.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

generate_allele_other_strand_special <- function(allele){

     list_of_bases <- unlist(strsplit(allele, ""))

     final_string <- ""

     for(base in list_of_bases){

       if(base == "A"){
          final_string =paste(final_string, "T", sep = "")
       }

       if(base == "T"){
          final_string = paste(final_string, "A", sep = "")
       }

       if(base == "C"){
          final_string =paste(final_string, "G", sep = "")
       }

       if(base == "G"){
          final_string =paste(final_string, "C", sep = "")
       }

     }

    return(final_string)

}


generate_allele_other_strand <- function(allele){

   if(nchar(allele) == 1){

     if(allele == "A"){
        return("T")
     }

     if(allele == "T"){
        return("A")
     }

     if(allele == "C"){
        return("G")
     }

     if(allele == "G"){
        return("C")
     }

   } else {

     insertion <- generate_allele_other_strand_special(allele) 
     return(insertion)

   }

}


##############
#LOADING DATA#
##############

data_4_clumping <- fread("/emc/cbmr/users/zlc436/HERMINA_PhD4_Project/CODE/1_Curating_data/1_Curating_FIadjBMI_data/data_4_clumping/independent_FIadjBMI.txt")

data_4_clumping$query_alleles_1 <- paste(data_4_clumping$effect_allele, "_", data_4_clumping$other_allele, sep = "")
data_4_clumping$query_alleles_2 <- paste(data_4_clumping$other_allele, "_", data_4_clumping$effect_allele, sep = "")

data_4_clumping$effect_allele_other_strand <- as.character(unlist(sapply(data_4_clumping$effect_allele, generate_allele_other_strand)))
data_4_clumping$other_allele_other_strand <- as.character(unlist(sapply(data_4_clumping$other_allele, generate_allele_other_strand)))

data_4_clumping$query_alleles_3 <- paste(data_4_clumping$effect_allele_other_strand, "_", data_4_clumping$other_allele_other_strand, sep = "")
data_4_clumping$query_alleles_4 <- paste(data_4_clumping$other_allele_other_strand, "_", data_4_clumping$effect_allele_other_strand, sep = "")

#And now let''s get the bim data. I could make a loop, but it is easier to just do it 22 times.

bim_1 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr1/chr1_1000G_phase_3_v5_EUR_maf.bim")
bim_2 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr2/chr2_1000G_phase_3_v5_EUR_maf.bim")
bim_3 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr3/chr3_1000G_phase_3_v5_EUR_maf.bim")
bim_4 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr4/chr4_1000G_phase_3_v5_EUR_maf.bim")
bim_5 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr5/chr5_1000G_phase_3_v5_EUR_maf.bim")
bim_6 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr6/chr6_1000G_phase_3_v5_EUR_maf.bim")
bim_7 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr7/chr7_1000G_phase_3_v5_EUR_maf.bim")
bim_8 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr8/chr8_1000G_phase_3_v5_EUR_maf.bim")
bim_9 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr9/chr9_1000G_phase_3_v5_EUR_maf.bim")
bim_10 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr10/chr10_1000G_phase_3_v5_EUR_maf.bim")
bim_11 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr11/chr11_1000G_phase_3_v5_EUR_maf.bim")
bim_12 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr12/chr12_1000G_phase_3_v5_EUR_maf.bim")
bim_13 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr13/chr13_1000G_phase_3_v5_EUR_maf.bim")
bim_14 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr14/chr14_1000G_phase_3_v5_EUR_maf.bim")
bim_15 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr15/chr15_1000G_phase_3_v5_EUR_maf.bim")
bim_16 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr16/chr16_1000G_phase_3_v5_EUR_maf.bim")
bim_17 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr17/chr17_1000G_phase_3_v5_EUR_maf.bim")
bim_18 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr18/chr18_1000G_phase_3_v5_EUR_maf.bim")
bim_19 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr19/chr19_1000G_phase_3_v5_EUR_maf.bim")
bim_20 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr20/chr20_1000G_phase_3_v5_EUR_maf.bim")
bim_21 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr21/chr21_1000G_phase_3_v5_EUR_maf.bim")
bim_22 <- fread("/emc/cbmr/users/zlc436/1000G_Phase3_v5/CURATED_DATA/PLINK_DATA/chr22/chr22_1000G_phase_3_v5_EUR_maf.bim")

#Let's merge them all:

bim_end <- rbind(bim_1, bim_2, bim_3, bim_4, bim_5, bim_6, bim_7, bim_8, bim_9,
                 bim_10, bim_11, bim_12, bim_13, bim_14, bim_15, bim_16, bim_17,
                 bim_18, bim_19, bim_20, bim_21, bim_22)

#Perfect.
#Now we need to get the chr_pos to match them:

chr_ <- paste("chr", bim_end$V1, sep = "")

bim_end$chr_pos <- paste(chr_, bim_end$V4, sep = ":")

############################################
#THIS SHOULD BE OK, let's do some checks...#
############################################

head(bim_end) #perfect.

##########################
#Let's now match the data#
##########################

bim_match <- bim_end[which(bim_end$chr_pos%in%data_4_clumping$SNP),] #42/44

##########################################
#For the time being let's do the matching#
##########################################

#Some SNPs are missing, let's check them:

data_4_clumping_missed <- data_4_clumping[which(!(data_4_clumping$SNP%in%bim_match$chr_pos)),] 

#One of them is an insertion, and the other one is a palindrome.
#Let's check in PhenoScanner: exactly the MAF < 0.01. 
#Hence why it is removed. Awesome.

data_4_clumping_ordered <- data_4_clumping[which(data_4_clumping$SNP%in%bim_match$chr_pos),]

#We have duplicates!! #but there is no need to worry. We can always match by alleles too. 
#that would make them disappear.
#Before doing that we are going to generate the other set of query alleles.

bim_match$ref_alleles <- paste(bim_match$V6, "_", bim_match$V5, sep = "") 

#Let's loop and check which can be obtained:

index_vect <- c()

for(i in seq(1,length(data_4_clumping_ordered$SNP))){

   print(i)

   index_end <- which(bim_match$chr_pos%in%data_4_clumping_ordered$SNP[i] & bim_match$ref_alleles == data_4_clumping_ordered$query_alleles_1[i] | 
                      bim_match$chr_pos%in%data_4_clumping_ordered$SNP[i] & bim_match$ref_alleles == data_4_clumping_ordered$query_alleles_2[i] |
                      bim_match$chr_pos%in%data_4_clumping_ordered$SNP[i] & bim_match$ref_alleles == data_4_clumping_ordered$query_alleles_3[i] | 
                      bim_match$chr_pos%in%data_4_clumping_ordered$SNP[i] & bim_match$ref_alleles == data_4_clumping_ordered$query_alleles_4[i] )

   if(is_empty(index_end) == FALSE){

       index_vect <- c(index_vect, index_end)

  }

}

bim_match_ordered <- bim_match[index_vect,] #this worked perfectly. There are 4 that do not match perfectly since they are insertions or deletions in the 1000G data.

data_4_clumping_ordered <- data_4_clumping[which(data_4_clumping$SNP%in%bim_match_ordered$chr_pos),]

data_4_clumping_ordered <- data_4_clumping_ordered[order(match(data_4_clumping_ordered$SNP, bim_match_ordered$chr_pos)),]

#Let's do some quick checkio:

print(length(which(bim_match_ordered$chr_pos != data_4_clumping_ordered$SNP))) #perfect match
print(length(which(bim_match_ordered$chr_pos == data_4_clumping_ordered$SNP))) #perfect match

##############################
#NOW LET'S GENERATE TWO FILES#
##############################

#We are gonna generate the file that we need to input and then a freaking dictionary to match the ID with chr_pos.
#Actually it is not extremely complicated to do the matching with some good string action, but it is easier this way.
#Also it is easier to check for mistakes if we do it this way. WIN WIN.

data_4_clumping_copy <- data_4_clumping_ordered 

data_4_clumping_copy$SNP <- bim_match_ordered$V2

write.table(data_4_clumping_copy, "/emc/cbmr/users/zlc436/HERMINA_PhD4_Project/CODE/1_Curating_data/1_Curating_FIadjBMI_data/data_4_clumping/updated_independent_FIadjBMI.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = " ")

#And now the dictionary:

dict_ <- bim_match_ordered %>%
   select(V2, chr_pos)

colnames(dict_) <- c("ID", "chr_pos")

write.table(dict_, "/emc/cbmr/users/zlc436/HERMINA_PhD4_Project/CODE/1_Curating_data/1_Curating_FIadjBMI_data/data_4_clumping/1000G_SNPID_dictionary.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = " ")



