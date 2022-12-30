##############
#INTRODUCTION#
##############

#This code will serve to parse the data of FI data that will be used in the GRS and 2SMR analysis.

#What we are going to do is just load the 19 FI variants and try our best.

###################
#LOADING LIBRARIES#
###################

library(tidyverse)
library(xlsx)
library(data.table)
library(jsonlite)
library(httr)

################
#LOAD FUNCTIONS#
################

ld_proxies <- function(snp){
  #Enter vector of SNPs it will output all its LD friends from 1KG European LD friends with r2 > 0.8
  #in 500kb.
  
  fake_df <- t(as.data.frame(c("d_prime", "variation2", "population_name", "r2", "variation1")))
  
  colnames(fake_df) <- fake_df[1,]
  rownames(fake_df) <- c(1)
  
  #Setting the server:
  
  server <- "http://grch37.rest.ensembl.org"
  
  for(i in snp){
    
    ext_1 <- paste("/ld/human/", i, sep = "")
    ext_2 <- paste(ext_1, "/1000GENOMES:phase_3:EUR", sep = "")
    
    r <- GET(paste(server, ext_2, sep = ""), content_type("application/json"))
    new <- fromJSON(toJSON(content(r)))
    
    fake_df <- rbind(fake_df, new)
    
  }
  
  #Now filtering for those that are in high LD:
  
  final_df <- fake_df[which(as.numeric(fake_df$r2) > 0.8),] #The NAs by coercion are the rows from the fake_df, ignore them!
  
  return(final_df)
  
}

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

##################################
#1. LOADING LAGOU ET AL 2021 DATA#
##################################

#First let's define the 19 variants that are previously reported.

lagou_snps <- c("rs2820436",                 "rs780094",                 "rs1530559",                 "rs10195252",
                "rs2943645",                 "rs17036328",                 "rs3822072",                 "rs9884482",
                "rs6822892",                 "rs4865796",                 "rs459193",                 "rs6912327",
                "rs2745353",                 "rs1167800",                 "rs983309",                 "rs7903146",
                "rs35767",                 "rs1421085",                 "rs731839")

#Now let's obtain the data for FI:

FI_lagou <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/RAW_DATA/Glycemic_traits_raw_data/FI_combined_1000G_density.txt.gz")

#We are going to remove those that are from ssIMP since they have been generated with 1000G data.

FI_lagou <- FI_lagou[which(FI_lagou$source != "SSIMP"),]

FI_lagou_df <- FI_lagou[which(FI_lagou$rsid%in%lagou_snps),] #perfect.

#Let's check that the data does not come from ssIMP which is trans-ancestry

View(FI_lagou_df) #PERFECT. All GWAS_reinserted = all from European data.

#They have everything I wanted.
#Now the only thing that we need to do is check if they are independent or not.

#We need to get the positions before we do anything:

pheno_data <- phenoscanner::phenoscanner(FI_lagou_df$rsid)$snps #one snp missing, but I think we know who that fellow is.

FI_lagou_df_1 <- FI_lagou_df[which(FI_lagou_df$rsid%in%pheno_data$snp),]
FI_lagou_df_2 <- FI_lagou_df[which(!(FI_lagou_df$rsid%in%pheno_data$snp)),]

pheno_data <- pheno_data[order(match(pheno_data$snp, FI_lagou_df_1$rsid)),]
length(which(pheno_data$snp == FI_lagou_df_1$rsid)) #all of them

FI_lagou_df_1$chr <- pheno_data$chr
FI_lagou_df_1$pos <- pheno_data$pos_hg19

#Let's put the other fella manually, who is our dear friend rs6822892.

FI_lagou_df_2$chr <- 4
FI_lagou_df_2$pos <- 157734675

FI_lagou_df <- rbind(FI_lagou_df_1, FI_lagou_df_2)

exposure_4_clumping <- FI_lagou_df %>%
  select("rsid", "p-value", "chr", "pos")

colnames(exposure_4_clumping) <- c("rsid", "pval", "chr", "pos")

exposure_ <- ieugwasr::ld_clump_local(exposure_4_clumping, bfile = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/2SMR/1000G/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/R-4.1.2/library/plinkbinr/bin/plink_Windows.exe", clump_kb = 500, clump_r2 = 0.01, clump_p = 0.05) 

#Only one variant got removed, which is the one that is not present in the 1000G
#Our old fella!
#That means that we can proceed without an issue!!!

exposure_manuscript_curated <- FI_lagou_df %>%
  select("rsid", "chr", "pos", "a2", "a1", "eaf", "beta", "se", "p-value", "n")

colnames(exposure_manuscript_curated) <- c("rsID", "Chromosome", "Position_build_37", "Effect_Allele", "Other_Allele", "EAF", "Beta", "SE", "P", "N")

#And now we save! We know there is no SNP in the MHC region here, so all is good.

write.csv(exposure_manuscript_curated, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/FI/CLEAN_19_FI_GRS_SNPs_Manuscript.csv", quote = FALSE, row.names = FALSE)
write.table(exposure_manuscript_curated, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/FI/CLEAN_19_FI_SNPs_Manuscript.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
