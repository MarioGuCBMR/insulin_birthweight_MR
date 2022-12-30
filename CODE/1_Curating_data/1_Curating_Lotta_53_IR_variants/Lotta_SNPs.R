##############
#INTRODUCTION#
##############

#This code is to obtain the data from Lotta et al:

###########
#Libraries#
###########

library(tidyverse)
library(data.table)

###################
#Loading functions#
###################

chr_parser <- function(chr_pos){
  #Take a ChrX:YYYYY format and return X.
  
  tmp <- as.character(strsplit(chr_pos, ":")[[1]])[1]
  
  tmp_2 <- as.numeric(as.character(strsplit(tmp, "Chr")[[1]]))[2]
  
  return(tmp_2)
  
}

pos_parser <- function(chr_pos){
  #Take a ChrX:YYYYY format and return X.
  
  tmp <- as.numeric(as.character(strsplit(chr_pos, ":")[[1]]))[2]
  
  return(tmp)
  
}

#########################################
#Reading SNPs and getting the data ready#
#########################################

exposure <- read.csv("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/FINAL_PIPELINE/RAW_DATA/LOTTA_SNPs_FIAdjBMI.csv", sep = ";")

colnames(exposure) <- c("rsid", "Chr_pos", "effect_allele", "OA", "beta", "pval")

exposure$rsid <- as.character(exposure$rsid)
exposure$pval <- as.numeric(exposure$pval)
library(stringr)
exposure$effect_allele <- str_replace_all(exposure$effect_allele, fixed(" "), "")
exposure$rsid <- str_replace_all(exposure$rsid, fixed(" "), "")

##############################################################
#We are gonna assume that this SNP is independent of the rest#
##############################################################

#Finding independent SNPs with stricter values:
#library(TwoSampleMR)

exposure_ <- ieugwasr::ld_clump_local(exposure, bfile = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/2SMR/1000G/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/R-4.1.2/library/plinkbinr/bin/plink_Windows.exe", clump_kb = 500, clump_r2 = 0.01, clump_p = 0.05) 

#rs6822892 because it is triallelic. 
#We have to check and see if it in 500kb distance of the others.
#If not, then there is no need to do any weird clumping.

exposure <- exposure[order(exposure$Chr_pos),]

157734675-89741269 >500000 #TRUE
157734675-3480136 >500000 #TRUE

#Then, there is no need to worry.
#We can use these fellas as much as we want.
#The main issue is that we need the SE for MR
#So we are going to use the original betas, not the ones reported this time.
#Unless Tuomas and Hermina suggest it.

###############################################
#Checking where do these bad boys come from...#
###############################################

#We are going to get Manning and Scott data and present the data as we have it there.
#Since we have both data in the same dataframe, I am going to use vector with the SNPs to properly check for them:

#Now we get the original betas, SE, and pvals:

scott <- fread("C:/Users/zlc436/Desktop/Data_4_2SMR/53_IR_Lotta/MAGIC_Metabochip_Public_data_release_25Jan/MAGIC_Scott_et_al_FI_adjBMI_Jan2013.txt")
manning <- fread("C:/Users/zlc436/Desktop/Data_4_2SMR/53_IR_Lotta/MAGIC_Manning_et_al_lnFastingInsulin_MainEffect.txt.gz")

exposure_scott <- scott[which(scott$snp%in%exposure$rsid),] #46
exposure_manning <- manning[which(manning$Snp%in%exposure$rsid),] #53

#We know that Lotta used Scott as primary and, if lacking, they used manning.
#We will see that those in manning have lower betas than Scott's:

scott_test <- exposure_scott[seq(1,5),]
manning_test <- exposure_manning[which(exposure_manning$Snp%in%scott_test$snp),]

manning_test <- manning_test[order(match(manning_test$Snp, scott_test$snp)),]

manning_test

#Snp effect_allele other_allele   maf MainEffects MainSE     MainP
#1:  rs1011685             t            c 0.133     -0.0130 0.0052 1.331e-02
#2: rs10195252             t            c 0.442      0.0130 0.0032 3.291e-05
#3:  rs1045241             t            c 0.265     -0.0040 0.0035 2.613e-01
#4: rs10995441             t            g 0.239     -0.0059 0.0039 1.312e-01
#5: rs11231693             a            g 0.036      0.0120 0.0067 6.930e-02
#BMIadjMainEffects BMIadjMainSE BMIadjMainP
#1:           -0.0093       0.0044   3.455e-02
#2:            0.0170       0.0027   1.482e-10
#3:           -0.0052       0.0030   8.726e-02
#4:           -0.0061       0.0033   6.902e-02
#5:            0.0170       0.0057   3.355e-03

scott_test

#snp effect_allele other_allele    beta     se    pvalue   maf
#1:  rs1011685             T            C -0.0110 0.0034 9.870e-04 0.133
#2: rs10195252             C            T -0.0170 0.0021 1.260e-16 0.442
#3:  rs1045241             T            C -0.0069 0.0023 2.053e-03 0.265
#4: rs10995441             G            T  0.0084 0.0025 8.500e-04 0.239
#5: rs11231693             A            G  0.0210 0.0043 7.190e-07 0.036

exposure[which(exposure$rsid%in%scott_test$snp),]

#rsid        Chr_pos effect_allele OA  beta     pval     id chr       pos    maf
#38 rs10995441 Chr10:64869239             G  T 0.014 8.50e-04 gakk7i  10  64869239 0.2485
#39 rs11231693 Chr11:63862612             A  G 0.036 7.19e-07 gakk7i  11  63862612 0.0666
#2  rs10195252 Chr2:165513091             T  C 0.029 1.26e-16 gakk7i   2 165513091 0.4423
#24  rs1045241 Chr5:118729286             C  T 0.012 2.00e-03 gakk7i   5 118729286 0.2604
#34  rs1011685  Chr8:19830769             C  T 0.019 9.80e-04 gakk7i   8  19830769 0.1302

#Exactly!!! The pvals match mostly Scott...
#but...

exposure_manning_rest <- exposure_manning[which(!(exposure_manning$Snp%in%exposure_scott$snp)),]

######################
#Checking the rest...#
######################

exposure_manning_rest

#Snp effect_allele other_allele   maf MainEffects MainSE     MainP
#1: rs11130329             a            c 0.159      0.0075 0.0048 0.1146000
#2: rs17402950             a            g 0.042     -0.0160 0.0067 0.0160800
#3:  rs7227237             t            c 0.186     -0.0120 0.0036 0.0013520
#4:  rs7323406             a            g 0.291      0.0110 0.0035 0.0010860
#5:  rs8032586             t            c 0.080     -0.0140 0.0047 0.0039020
#6:  rs8101064             t            c 0.018      0.0180 0.0083 0.0294700
#7:   rs972283             a            g 0.450     -0.0120 0.0032 0.0001564
#BMIadjMainEffects BMIadjMainSE BMIadjMainP
#1:            0.0140       0.0040   5.117e-04
#2:           -0.0160       0.0057   4.682e-03
#3:           -0.0099       0.0031   1.331e-03
#4:            0.0090       0.0030   2.674e-03
#5:           -0.0110       0.0040   4.568e-03
#6:            0.0250       0.0072   6.226e-04
#7:           -0.0130       0.0027   4.408e-06

exposure[which(exposure$rsid%in%exposure_manning_rest$Snp),]

#rsid         Chr_pos effect_allele OA  beta     pval     id chr       pos    maf
#40 rs17402950  Chr12:14571671             G  A 0.027 4.70e-03 gakk7i  12  14571671 0.0547
#43  rs7323406 Chr13:111628195             A  G 0.015 2.70e-03 gakk7i  13 111628195 0.2783
#45  rs8032586  Chr15:73081067             C  T 0.019 4.60e-03 gakk7i  15  73081067 0.1093
#47  rs7227237  Chr18:47174679             C  T 0.017 1.30e-03 gakk7i  18  47174679 0.2276
#48  rs8101064   Chr19:7293119             T  C 0.042 6.20e-04 gakk7i  19   7293119 0.0467
#18 rs11130329   Chr3:52896855             A  C 0.024 5.10e-04 gakk7i   3  52896855 0.1332
#32   rs972283  Chr7:130466854             G  A 0.022 4.41e-06 gakk7i   7 130466854 0.4394

#The p-value match perfectly!

#First we need to make sure that the alleles match in both cases:

new_a1 <- ifelse(exposure_scott$beta < 0, exposure_scott$other_allele, exposure_scott$effect_allele)
new_a2 <- ifelse(exposure_scott$beta < 0, exposure_scott$effect_allele, exposure_scott$other_allele)
new_beta <- ifelse(exposure_scott$beta < 0, exposure_scott$beta*(-1), exposure_scott$beta)

exposure_scott$A1 <- new_a1
exposure_scott$A2 <- new_a2
exposure_scott$BETA <- new_beta

#Getting ready for scott:

exposure_scott_match <- exposure[which(exposure$rsid%in%exposure_scott$snp),]
exposure_scott_match <- exposure_scott_match[order(match(exposure_scott_match$rsid, exposure_scott$snp)),]

which(exposure_scott_match$rsid != exposure_scott$snp) #perfect!
which(exposure_scott_match$effect_allele != exposure_scott$A1) #perfect too. Everything is perfect.
which(exposure_scott_match$effect_allele == exposure_scott$A1) #perfect too. Everything is perfect.

exposure_scott_match$original_beta <- exposure_scott$BETA
exposure_scott_match$original_se <- exposure_scott$se
exposure_scott_match$original_pval <- exposure_scott$pvalue
exposure_scott_match$original_maf <- exposure_scott$maf

#And now for Manning:

new_a1 <- ifelse(exposure_manning_rest$BMIadjMainEffects < 0, toupper(as.character(exposure_manning_rest$other_allele)) , toupper(as.character(exposure_manning_rest$effect_allele)))
new_a2 <- ifelse(exposure_manning_rest$BMIadjMainEffects < 0, toupper(as.character(exposure_manning_rest$effect_allele)) , toupper(as.character(exposure_manning_rest$other_allele)))
new_beta <- ifelse(exposure_manning_rest$BMIadjMainEffects < 0, exposure_manning_rest$BMIadjMainEffects*(-1), exposure_manning_rest$BMIadjMainEffects)

exposure_manning_rest$A1 <- new_a1
exposure_manning_rest$A2 <- new_a2
exposure_manning_rest$BETA <- new_beta

#Getting ready for manning:

exposure_manning_rest_match <- exposure[which(exposure$rsid%in%exposure_manning_rest$Snp),]
exposure_manning_rest_match <- exposure_manning_rest_match[order(match(exposure_manning_rest_match$rsid, exposure_manning_rest$Snp)),]

which(exposure_manning_rest_match$rsid != exposure_manning_rest$Snp) #perfect!
which(exposure_manning_rest_match$effect_allele != exposure_manning_rest$A1) #perfect too. Everything is perfect.
which(exposure_manning_rest_match$effect_allele == exposure_manning_rest$A1) #perfect too. Everything is perfect.

exposure_manning_rest_match$original_beta <- exposure_manning_rest$BETA
exposure_manning_rest_match$original_se <- exposure_manning_rest$BMIadjMainSE
exposure_manning_rest_match$original_pval <- exposure_manning_rest$BMIadjMainP
exposure_manning_rest_match$original_maf <- exposure_manning_rest$maf

##########################################################################
#Now we arrange per chromosome independently and merge with Manning first#
##########################################################################

exposure_manning_rest_match <- exposure_manning_rest_match[order(exposure_manning_rest_match$Chr_pos),]
exposure_scott_match <- exposure_scott_match[order(exposure_scott_match$Chr_pos),]

#Now we merge them all:

exposure_manning_rest_match$origin <- "Manning_et_al"
exposure_scott_match$origin <- "Scott_et_al"

exposure_manuscript_end <- rbind(exposure_manning_rest_match, exposure_scott_match)

###################################
#Let's createchr and pos variables#
###################################

exposure_manuscript_end$chr <- as.numeric(as.character(unlist(sapply(exposure_manuscript_end$Chr_pos, chr_parser))))
exposure_manuscript_end$pos <- as.numeric(as.character(unlist(sapply(exposure_manuscript_end$Chr_pos, pos_parser))))

#One of them is distorted because life is harsh, but that is no problem here.
#We query it in phenomescanner...

exposure_manuscript_end$pos[is.na(exposure_manuscript_end$pos)] <- phenoscanner::phenoscanner(exposure_manuscript_end$rsid[is.na(exposure_manuscript_end$pos)])$snps$pos_hg19

#And we check the results...

View(exposure_manuscript_end) #DAMN PERFECT.

#Finally, we need the effect allele frequency, if not we are going to have issues with the 2SMR. 

pheno_data <- phenoscanner::phenoscanner(exposure_manuscript_end$rsid)

pheno_data <- pheno_data$snps #52 SNPs because we are missing the damn missing SNP that is not in 1000G.

#In that case we should check it in the FBW data and see what allele has the closest EAF to solve it.
#In the meantime, we are going to take care of the rest:

exposure_manuscript_end_1 <- exposure_manuscript_end[which(exposure_manuscript_end$rsid%in%pheno_data$snp),]
pheno_data <- pheno_data[order(match(pheno_data$snp, exposure_manuscript_end_1$rsid)),]
length(which(exposure_manuscript_end_1$rsid == pheno_data$snp)) #all

#Let's go:

exposure_manuscript_end_1$EAF <- ifelse(exposure_manuscript_end_1$effect_allele == pheno_data$a1 & pheno_data$eur < 0.50, exposure_manuscript_end_1$original_maf, abs(1- exposure_manuscript_end_1$original_maf))

#Let's check that this has been done correctly:

head(exposure_manuscript_end_1) #worked like a freaking charm
head(pheno_data) #worked like a freaking charm.

#Only one more thing to do: check that last variant:

exposure_manuscript_end_2 <- exposure_manuscript_end[which(!(exposure_manuscript_end$rsid%in%pheno_data$snp)),]

#Let's read the FBW data and check:

fbw <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/OUTPUT/FBW/Birthweight2021_Curated.txt")

check <- fbw[which(fbw$rsID == exposure_manuscript_end_2$rsid),]

exposure_manuscript_end_2$EAF <- ifelse(exposure_manuscript_end_2$effect_allele == check$A1 & as.numeric(check$`EGG-frq`) < 50, exposure_manuscript_end_2$original_maf, abs(1- exposure_manuscript_end_2$original_maf))

#Just checked it: this work perfectly!

exposure_manuscript_end <- rbind(exposure_manuscript_end_1, exposure_manuscript_end_2)

#Finally, we are going to use the same columns as with FIadjBMI, but with origin data as well:

exposure_manuscript_curated <- exposure_manuscript_end %>%
  select("rsid", "chr", "pos", "effect_allele", "OA", "EAF", "original_maf", "beta", "original_beta", "original_se", "original_pval", "origin")

colnames(exposure_manuscript_curated) <- c("rsID", "Chromosome", "Position_build_37", "Effect_Allele", "Other_Allele", "EAF", "MAF", "Beta_SD", "Beta", "SE", "P", "origin")

#Finally, and this is really important, we do not have sample size, we need to get them from the MAGIC CONSORTIUM page.
#Damn, this is a huge pain.

#From Scott's: Fasting insulin results are for ln-transformed fasting insulin as the outcome and are adjusted for age, sex and are reported both with and without BMI adjustment. These results are from up to 108,557 individuals from 56 studies.
#From Manning's: he insulin results accounting for BMI are from an analysis of 26 studies in up to 51,750 non-diabetic participants.

exposure_manuscript_curated$N <- ifelse(exposure_manuscript_curated$origin == "Manning_et_al", 51750, 108557)

#Then we can save everything :)

write.csv(exposure_manuscript_curated, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/53_IR_Lotta/CLEAN_53_LOTTA_GRS_SNPs_Manuscript.csv", quote = FALSE, row.names = FALSE)
write.table(exposure_manuscript_curated, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/53_IR_Lotta/CLEAN_53_LOTTA_GRS_SNPs_Manuscript.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)