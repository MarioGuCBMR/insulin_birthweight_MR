##############
#INTRODUCTION#
##############

#This code is to produce a dataframe with rsID, chr_pos_37 and chr_pos_38 for the HAPMAP variants with LD scores
#used by LDSC to calculate genetic correlations.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

clean_liftover_data <- function(chr_pos_liftover){
  #Takes liftover data (chr:pos1-pos2) in a build and transforms it to chr:pos data
  
  tmp_pos <- strsplit(chr_pos_liftover, "-")[[1]][2]
  chr_tmp <- strsplit(chr_pos_liftover, ":")[[1]][1]
  
  tmp_pos_clean <- as.numeric(tmp_pos)-1
  
  chr_pos <- paste(chr_tmp, ":", tmp_pos_clean, sep = "")
  
  return(chr_pos)
  
  
}

##############
#Loading data#
##############

#To properly do this we are gonna load the literal LD scores, since those files have rsID and chr_pos.
#This data is obtained by installing the LDSC github. Three essential notes:

#A) This data is from 2015, so some SNPs are old. Be wary with that.
#B) This data was used by Warrington et al 2019 to generate genetic correlations. I have replicated those genetic correlations..., so it should be fine.
#c) The data from chromosome 6 does not present variants from the HLA region due to its weird structure.

path_ <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Mario_Tools/LDSC_reference_panel/RAW_DATA/"

ldsc_df <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Mario_Tools/LDSC_reference_panel/RAW_DATA/1.l2.ldscore.gz")

#And now we loop through the rest of chromsome

for(chr in seq(2,22)){
  
  ld_path <- paste(path_, chr, ".l2.ldscore.gz", sep = "") #we set the path for the next chr
  
  tmp_df <- fread(ld_path) #we load the new df
  
  ldsc_df <- rbind(ldsc_df, tmp_df) #we update the original df
  
}

#Now we have the dataframe with position in build 37 (checked in phenoscanner.)

##############################################
#Let's generate a file to upload in liftover.#
##############################################

options(scipen = 100) #important to not make liftover crazy!

liftover_data <- paste("chr", ldsc_df$CHR, ":", as.character(as.numeric(ldsc_df$BP)-1), "-", as.character(as.numeric(ldsc_df$BP)+1), sep = "")

#Now let's save the data and load into liftover:

ldsc_df$liftover <- liftover_data

fwrite(as.data.frame(ldsc_df$liftover), "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Mario_Tools/LDSC_reference_panel/RAW_DATA/ldsc_4_liftover.txt", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

#check <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Mario_Tools/LDSC_reference_panel/RAW_DATA/ldsc_4_liftover.txt",header = FALSE)

#We had failures (like my brother) so we are gonna remove them:

failure_data <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Mario_Tools/LDSC_reference_panel/RAW_DATA/failing_conversions.txt")

ldsc_df_clean <- ldsc_df[which(!(ldsc_df$liftover%in%failure_data$`#Partially deleted in new`)),]

fwrite(as.data.frame(ldsc_df_clean$liftover), "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Mario_Tools/LDSC_reference_panel/RAW_DATA/ldsc_4_liftover_clean.txt", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

#####################################
#Let's retrieve the data in build 38#
#####################################

ldsc_38 <- as.data.frame(read.table("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Mario_Tools/LDSC_reference_panel/CURATED_DATA/LDSC_build38.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

ldsc_df_clean$liftover_38 <- ldsc_38$V1

#Now let's get the positions in the usual format read by phenoscanner, 1000G among others: chr1:1000.

ldsc_df_clean$chr_pos_38 <- as.character(unlist(sapply(ldsc_df_clean$liftover_38, clean_liftover_data)))

################################################################################
#Get positions of build 37 in the right format and let's make a final dataframe#
################################################################################

ldsc_df_clean$chr_pos_37 <- paste("chr", ldsc_df_clean$CHR, ":", ldsc_df_clean$BP, sep = "")

#I did a couple of queries in chr1 and chr22. They match perfectly!!!

#The only thing left to do is get those that could not be recovered and add them in chr_pos_38 as "-" (they could not be recovered in PS for build 38 either!, but I could in build 37. They really do not exist there.)

ldsc_df_miss <- ldsc_df[which(ldsc_df$liftover%in%failure_data$`#Partially deleted in new`),]

#Now let's get the chr_pos format:

ldsc_df_miss$chr_pos_37 <- paste("chr", ldsc_df_miss$CHR, ":", ldsc_df_miss$BP, sep = "")

#And now let's do a dummy variable of chr_pos 38.

ldsc_df_miss$chr_pos_38 <- "-"

#PERFECT! 

#######################################
#Merging the clean and the missed data#
#######################################

ldsc_df_clean_end <- ldsc_df_clean %>%
  select(SNP, CHR, BP, CM, MAF, L2, chr_pos_37, chr_pos_38)

ldsc_df_miss_end <- ldsc_df_miss %>%
  select(SNP, CHR, BP, CM, MAF, L2, chr_pos_37, chr_pos_38)

ldsc_df_end <- rbind(ldsc_df_clean_end, ldsc_df_miss_end) #same rows as the original.

#And now let's just save this baby:

fwrite(ldsc_df_end, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Mario_Tools/LDSC_reference_panel/CURATED_DATA/LDSC_HM3_SNPs_build37_&_38.txt", sep = "\t")

check <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Mario_Tools/LDSC_reference_panel/CURATED_DATA/LDSC_HM3_SNPs_build37_&_38.txt")
