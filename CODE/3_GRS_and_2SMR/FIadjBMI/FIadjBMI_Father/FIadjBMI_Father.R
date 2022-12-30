##############
#INTRODUCTION#
##############

#This code will serve to parse the data of the full results of the European hits for 
#the glycemic traits found in supplementary table 2.

###################
#LOADING LIBRARIES#
###################

library(tidyverse)
library(xlsx)
library(data.table)
library(TwoSampleMR)
library(gtx)
library(meta)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(rmarkdown)

#########################
#Making space for memory#
#########################

memory.limit(size=8000000)

###################
#LOADING FUNCTIONS#
###################

match_ss_data <- function(query_chr_pos, query_alleles_1, query_alleles_2, bim_data){
  #This function takes a chr_pos and its alleles and checks whether they are present in bim data.
  #It returns the index of the bim data where the allele that you are supposed to use is located.
  
  #STEP 1: getting the data to test (it will be commented later)
  
  #query_chr_pos = "chr7:130465054"
  #query_alleles_1 = "T_G"
  #query_alleles_2 = "G_T"
  #bim_data = dictionary
  
  #STEP 2: Let's create the ref_alleles in bim_data:
  
  #bim_data$ref_alleles <- paste(bim_data$V6, "_", bim_data$V5, sep = "")
  
  #STEP 3: Now let's perform the query:
  
  bim_data_match <- bim_data[which(bim_data$chr_pos%in%query_chr_pos & bim_data$ref_alleles == query_alleles_1 | bim_data$chr_pos%in%query_chr_pos & bim_data$ref_alleles == query_alleles_2),]
  
  #Finally, let's get the index:  
  
  if(is_empty(bim_data_match) == TRUE){
    
    next()
    
  }
  
  index_end <- which(bim_data$chr_pos == bim_data_match$chr_pos & bim_data$ref_alleles == bim_data_match$ref_alleles)
  
  #And we return the index: 
  
  return(index_end)
  
}


mr_plots <- function(dat)
{
  require(TwoSampleMR)
  require(ggplot2)
  require(ggrepel)
  require(dplyr)
  
  temp <- subset(dat, outcome == outcome[1] & exposure == exposure[1])
  exposure_name <- temp$exposure[1]
  outcome_name <- temp$outcome[1]
  
  if(! "labels" %in% names(dat)) dat$labels <- NA
  
  exposure_units <- temp$units.exposure[1]
  outcome_units <- temp$units.outcome[1]
  
  mrs <- mr_singlesnp(temp, all_method=c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))
  mrl <- mr_leaveoneout(temp)
  
  mrs$index <- 1:nrow(mrs)
  mrl$index <- 1:nrow(mrl)
  
  mrs <- dplyr::arrange(merge(mrs, select(temp, SNP, labels), all.x=TRUE), index)
  mrl <- dplyr::arrange(merge(mrl, select(temp, SNP, labels), all.x=TRUE), index)
  
  mrres <- mr(temp, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))
  
  gridExtra::grid.arrange(
    mr_forest_plot(mrs)[[1]] +
      ggplot2::labs(
        title="a)",
        x=paste0("MR estimate (", outcome_units, " per ", exposure_units, ")")
      ),
    mr_scatter_plot(mrres, temp)[[1]] +
      ggplot2::labs(
        title="b)",
        x=paste0("SNP effect on ", exposure_name),
        y=paste0("SNP effect on ", outcome_name)
      ) +
      geom_label_repel(ggplot2::aes(label=labels), point.padding = unit(0.5, "lines")),
    mr_leaveoneout_plot(mrl)[[1]] +
      ggplot2::labs(
        title="c)",
        x=paste0("MR estimate (", outcome_units, " per ", exposure_units, ")"),
        y="Excluded variant"
      ),
    mr_funnel_plot(mrs)[[1]] +
      ggplot2::labs(title="d)") +
      ggplot2::theme(legend.position="none") +
      ggrepel::geom_label_repel(ggplot2::aes(label=labels), point.padding = unit(0.5, "lines")),
    ncol=2
  )
}

steiger_filtering_old <- function(dat)
{
  if(! "units.outcome" %in% names(dat))
  {
    dat$units.outcome <- NA
  }
  if(! "units.exposure" %in% names(dat))
  {
    dat$units.exposure <- NA
  }
  stopifnot(length(unique(dat$exposure)) == 1)
  stopifnot(length(unique(dat$outcome)) == 1)
  stopifnot(length(unique(dat$units.exposure)) == 1)
  stopifnot(length(unique(dat$units.outcome)) == 1)
  
  dat <- add_rsq(dat)
  
  st <- psych::r.test(
    n = dat$sample_size.exposure, 
    n2 = dat$sample_size.outcome, 
    r12 = sqrt(dat$rsq.exposure), 
    r34 = sqrt(dat$rsq.outcome)
  )
  dat$steiger_dir <- dat$rsq.exposure > dat$rsq.outcome
  dat$steiger_pval <- pnorm(-abs(st$z)) * 2
  
  return(dat)
}

###################
#STEP 1: LOAD DATA#
###################

#Get original data for exposure to obtain the samplesizes and good allele frequencies:

FIadjBMI_original <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Romain_Mario/RAW_DATA/FIAdjBMI_2021/MAGIC1000G_FI_EUR.tsv.gz")

FIadjBMI_original$chr_pos <- paste("chr", FIadjBMI_original$chromosome, ":", FIadjBMI_original$base_pair_location, sep = "")

fiadjbmi <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/FIadjBMI/FIadjBMI_independent_4_GRS.txt")
father <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/father/Birthweight_offspring_fathers2021_curated_FULL.txt")

fiadjbmi$chr_pos <- paste("chr", fiadjbmi$Chromosome, ":", fiadjbmi$Position_build_37, sep = "")

#######################
#STEP 2: MATCHING DATA#
#######################

#NOTE: the insertion/deletion are in the same format here!! (the position takes into account the amount of deletions, which is great.)

father_2_fiadjbmi <- father[which(father$rsID%in%fiadjbmi$rsID),] #30/31

#There is one missing!!

fiadjbmi[which(!(fiadjbmi$chr_pos%in%father_2_fiadjbmi$chr_pos_37)),] #we have two mismatches!

#rsID Chromosome Position_build_37 Effect_Allele Other_Allele   EAF   MAF
#1: rs77935490          2            630902             A            T 0.236 0.236
#Beta     SE         P      N     chr_pos
#1: 0.0137 0.0026 1.067e-08 151012 chr2:630902

#Let's take outliers for those missing:

#Find SNPs that are missing

missing_snps <- fiadjbmi[which(!(fiadjbmi$chr_pos%in%father_2_fiadjbmi$chr_pos_37)),] #3

#Query them in LDlink and check whether they are in 1000G

setwd("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/3_2SMR_and_GRS/FIadjBMI/FIadjBMI_Father/proxies/")

getwd()

LDlinkR::LDproxy_batch(snp = missing_snps$chr_pos, pop = "EUR", r2d = "r2", append = TRUE, token = "04cad4ca4374")

#Let's check whether we queried all SNPs.

proxies_ld_link <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/3_2SMR_and_GRS/FIadjBMI/FIadjBMI_Father/proxies/combined_query_snp_list.txt")

#Let's see which variants did we cover:

query_snps <- proxies_ld_link[which(proxies_ld_link$query_snp == proxies_ld_link$Coord & proxies_ld_link$Distance == 0 & proxies_ld_link$R2 == 1),]

#We covered all of the variants except that I/D.
#Which can be ignorend.

###########################################
#In any case let's clean the proxies first#
###########################################

#Hence let's first clean the proxies that we find with LD_link.
#Afterwards, we are going to check if we can get the proxies for this fella.

proxies_df_clean <- proxies_ld_link[which(!(proxies_ld_link$query_snp == proxies_ld_link$Coord & proxies_ld_link$Distance == 0 & proxies_ld_link$R2 == 1)),]

#Now we are only interested in those with a r2 > 0.8

proxies_df_clean_curated <- proxies_df_clean[which(proxies_df_clean$R2 > 0.8),]

#Just in case, let's check the MAF:

summary(proxies_df_clean_curated$MAF) #they are all perfect.

#Let's get the proxies in FIAdjBMI_

proxies_in_FIadjBMI <- FIadjBMI_original[which(FIadjBMI_original$chr_pos%in%proxies_df_clean_curated$Coord),] #230!! Plus some duplicates, it seems.

#We have duplicates, but we know those do not exist in LDLink so we are going to remove them:

#The duplicates rise from insertions and deletions.
#We are going to match first the proxies in father and then work with the mismatches.

proxies_in_father <- father[which(father$chr_pos_37%in%proxies_df_clean_curated$Coord),] #199

#Let's match those found in father with FIadjBMI:

proxies_in_FIadjBMI <- proxies_in_FIadjBMI[which(proxies_in_FIadjBMI$chr_pos%in%proxies_in_father$chr_pos_37),] #225/227!! We lose some!

proxies_in_father <- proxies_in_father[which(proxies_in_father$chr_pos_37%in%proxies_in_FIadjBMI$chr_pos),] #not changing.

#We have a severe case of duplicates.

#In any case, that makes my life easier:

missing_snp <- proxies_in_father[which(!(proxies_in_father$chr_pos_37%in%proxies_in_FIadjBMI$chr_pos)),] #none. It is a duplicate.
dupl_snp <- proxies_in_father$chr_pos_37[which(duplicated(proxies_in_father$chr_pos_37) == TRUE)]

#Let's check the correct allele.

proxies_in_father[which(proxies_in_father$chr_pos_37%in%dupl_snp),] #it is actually a tri-allelic SNP, but the G/T is not counted because G MAF = 0 in European. Hence why it is included in LDlink.

#Many of them are insertions and deletions. Or variants with exclamation marks,
#which I am never too sure what they mean. 

#Let's check them in FIadjBMI:

proxies_in_FIadjBMI[which(proxies_in_FIadjBMI$chr_pos%in%dupl_snp),]

#None of them match...
#I think the best is doing some cleaning.

##################################################################################
#CHECKING THAT THE ALLELES ARE OKAY BETWEEN REFERENCE PANEL, EXPOSURE AND OUTCOME#
##################################################################################

#To do so we need to do the following. 

#1. Get the final proxies' original data:

proxies_df_check <- proxies_df_clean[which(proxies_df_clean$Coord%in%proxies_in_father$chr_pos_37),]

check_fiadjbmi <- proxies_in_FIadjBMI
check_father <- proxies_in_father

#2. Let's check the alleles for FIAdjBMI and see if the mismatches are due to formating and not due to wrong SNP.

check_fiadjbmi$alleles_1 <- paste("(", check_fiadjbmi$effect_allele, "/", check_fiadjbmi$other_allele, ")", sep = "")
check_fiadjbmi$alleles_2 <- paste("(", check_fiadjbmi$other_allele, "/", check_fiadjbmi$effect_allele, ")", sep = "")

proxies_df_check$chr_pos <- proxies_df_check$Coord
proxies_df_check$ref_alleles <- proxies_df_check$Alleles

#And now check everything!

index_vect <- c()

for(i in seq(1,length(check_fiadjbmi$chr_pos))){
  
  print(i)
  
  index_tmp <- as.numeric(as.character(unlist(sapply(check_fiadjbmi$chr_pos[i], match_ss_data, query_alleles_1 = check_fiadjbmi$alleles_1[i], query_alleles_2 = check_fiadjbmi$alleles_2[i], bim_data = proxies_df_check))))
  
  index_vect <- c(index_vect, index_tmp)
  
}

length(index_vect) #199! 

proxies_df_check_clean <- proxies_df_check[index_vect,]
proxies_df_check_wrong <- proxies_df_check[-index_vect,]

proxies_in_FIadjBMI_clean <- proxies_in_FIadjBMI[which(proxies_in_FIadjBMI$chr_pos%in%proxies_df_check_clean$Coord),] #perfect
proxies_in_FIadjBMI_wrong <- proxies_in_FIadjBMI[which(!(proxies_in_FIadjBMI$chr_pos%in%proxies_df_check_clean$Coord)),] #we removed the weird ones.

#There is a duplicate that we need to solve.
#But still let's check which one is right and the wrong one:

proxies_in_FIadjBMI_clean[which(duplicated(proxies_in_FIadjBMI_clean$chr_pos) == TRUE),]

#Let's check the right allele.
#This SNP should be in both: the clean and the wrong sets:

proxies_df_check_clean[which(proxies_df_check_clean$chr_pos == "chr2:630903"),] #A/G #It is the one we see as duplicate!

proxies_in_FIadjBMI_clean[which(proxies_in_FIadjBMI_clean$chr_pos == "chr2:630903"),]

#The insertion and deletion!

proxies_in_FIadjBMI <- proxies_in_FIadjBMI_clean[which(proxies_in_FIadjBMI_clean$effect_allele != "D"),] #199! Done.

#Let's match the proxies_father and proxies_df_check with this new set!

proxies_df_check <- proxies_df_check[which(proxies_df_check$Coord%in%proxies_in_FIadjBMI$chr_pos),] #204
proxies_in_father <- proxies_in_father[which(proxies_in_father$chr_pos_37%in%proxies_in_FIadjBMI$chr_pos),] #208!

#Everything is perfect!!

#3. Let's check the alleles for father and see if the mismatches are due to formating and not due to wrong SNP.

check_father$alleles_1 <- paste("(", check_father$A1, "/", check_father$A0, ")", sep = "")
check_father$alleles_2 <- paste("(", check_father$A0, "/", check_father$A1, ")", sep = "")

#And now check everything!

index_vect <- c()

for(i in seq(1,length(check_father$chr_pos_37))){
  
  print(i)
  
  index_tmp <- as.numeric(as.character(unlist(sapply(check_father$chr_pos_37[i], match_ss_data, query_alleles_1 = check_father$alleles_1[i], query_alleles_2 = check_father$alleles_2[i], bim_data = proxies_df_check))))
  
  index_vect <- c(index_vect, index_tmp)
  
}

length(index_vect) #199/199! Perfect.

proxies_df_check_clean <- proxies_df_check[index_vect,]
proxies_df_check_wrong <- proxies_df_check[-index_vect,]

proxies_in_father_clean <- proxies_in_father[which(proxies_in_father$chr_pos_37%in%proxies_df_check_clean$Coord),] #199/199!! PERFECT FIT.
proxies_in_father_wrong <- proxies_in_father[which(!(proxies_in_father$chr_pos_37%in%proxies_df_check_clean$Coord)),] #only one insertion.

#Let's check the wrong ones and see what happens:

proxies_in_father_wrong #none because they were erased
proxies_df_check_wrong #none because they were erased

#In any case, we are don people!

proxies_df_check <- proxies_df_check_clean
proxies_in_father <- proxies_in_father_clean
proxies_in_FIadjBMI <- proxies_in_FIadjBMI[which(proxies_in_FIadjBMI$chr_pos%in%proxies_in_father$chr_pos_37),] #199

#PERFECT! ALLELES CHECKED. ALL GOOD.
#Now we can proceed to perform the clumping on FIadjBMI!

########################################################
#CLUMPING TO OBTAIN THE PROXIES OF THE MISSING VARIANTS#
########################################################

fiadjbmi_in_father <- fiadjbmi[which(fiadjbmi$chr_pos%in%father_2_fiadjbmi$chr_pos_37),]

fiadjbmi_in_father_4_clumping <- fiadjbmi_in_father %>%
  select(rsID, P)

colnames(fiadjbmi_in_father_4_clumping) <- c("rsid", "pval")

#We need RSID so...

proxies_df_check <- proxies_df_check[order(match(proxies_df_check$Coord, proxies_in_FIadjBMI$chr_pos)),]

length(which(proxies_df_check$chr_pos == proxies_in_FIadjBMI$chr_pos)) #all of them.

proxies_in_FIadjBMI$rsid <- proxies_df_check$RS_Number

proxies_4_clumping <- proxies_in_FIadjBMI %>%
  select(rsid, p_value)

colnames(proxies_4_clumping) <- c("rsid", "pval")

variants_4_clumping <- rbind(proxies_4_clumping, fiadjbmi_in_father_4_clumping)

#Let's proceed to the clumping being really wary that at least two SNPs are going to be removed.
#Due to not being in the reference panel.

independence_test_end <- ieugwasr::ld_clump_local(variants_4_clumping, bfile = "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Sedentary_BMI/CODE_4_replication/2SMR/1000G/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/R-4.1.2/library/plinkbinr/bin/plink_Windows.exe", clump_kb = 500, clump_r2 = 0.01, clump_p = 0.05) 

#PERFECT. Only two SNPs are missing from the dataset, but not the one I was expecting.
#The reason behind this is that rs77935490 is a missing SNP in the data.

#Let's check the SNPs that we have obtained:

independence_test_end[which(!(independence_test_end$rsid%in%fiadjbmi_in_father_4_clumping$rsid)),]

#We have SNP rs5017304, a SNP that was added. If it is in chromosome 2, that is the substitute to our fella.
#This is perfect.
#We know that one SNP was ignored: rs74832795

#There is one SNP that is missing, though, and we do not know why: rs74832795. Let's check why.

proxies_df_check[which(proxies_df_check$RS_Number == "rs74832795"),] #another of the candidates.

#Exactly, it is a proxy of our missing SNP.
#Now we need to check whether any of the proxies for this SNP have been recovered.
#It seems two SNPs were recovered. We need to remember that one SNP was lost.
#Let's check if the missing SNP should have been the lead SNP.

#To check we need to get the proxy for rs77935490 and compare:

#The proxy we have is rs5017304 (chr2:2:630903)
#The one we lost is rs74832795 (chr2:632592)

#Let's see which one is the lead SNP in FIadjBMI

check_lead <- c("chr2:630903", "chr2:632592")

FIadjBMI_original[which(FIadjBMI_original$chr_pos%in%check_lead),]

#chromosome base_pair_location effect_allele other_allele
#1:          2             630903             D            I
#2:          2             630903             A            G
#3:          2             632592             A            G
#effect_allele_frequency    beta standard_error   p_value sample_size
#1:                   0.240  0.0138         0.0029 2.475e-07      119003
#2:                   0.235  0.0138         0.0026 1.344e-08      148101
#3:                   0.770 -0.0119         0.0025 2.896e-07      150846
#chr_pos
#1: chr2:630903
#2: chr2:630903
#3: chr2:632592

#This implies that the lead SNP is chosen correctly despite one variant being chosen over the other.
#We have recovered the loci for the analysis!
#Jesus Christ, we are almost there!!!

#################################################
#Let's get the FIadjBMI for all clumped variants#
#################################################

proxies_recovered <- c("chr2:630903") #chr_pos obtain from checking in FBW which has both. WHAT A PAIN.

proxies_in_FIadjBMI_selected <- proxies_in_FIadjBMI[which(proxies_in_FIadjBMI$chr_pos%in%proxies_recovered),]

#With this we can merge the data, though is wise to clean the separately and merge them later..., so we will do that!

##########################
#STEP 3: CLEANING COLUMNS#
##########################

#First for those in FIadjBMI that were in FBW:

exposure_fiadjbmi_1 <- fiadjbmi_in_father %>%
  select("rsID", "Effect_Allele", "Other_Allele", "EAF",
         "Beta", "SE", "P", "N", "chr_pos")

colnames(exposure_fiadjbmi_1) <- c("SNP", "effect_allele.exposure", "other_allele.exposure",
                                   "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "chr_pos")

#Now for those recovered:

exposure_fiadjbmi_2 <- proxies_in_FIadjBMI_selected %>%
  select(rsid, effect_allele, other_allele, effect_allele_frequency,
         beta, standard_error,  p_value, sample_size, chr_pos)

colnames(exposure_fiadjbmi_2)  <- c("SNP", "effect_allele.exposure", "other_allele.exposure",
                                    "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "chr_pos")

exposure_fiadjbmi <- rbind(exposure_fiadjbmi_1, exposure_fiadjbmi_2)

#And now the same for the father data:

outcome_father <- father[which(father$rsID%in%exposure_fiadjbmi$SNP),] #31

outcome_father <- outcome_father %>%
  select(rsID, A1,A0, `IS-frq`, `Beta-A1`, SE, P, `IS-info`)

colnames(outcome_father) <- c("SNP", "effect_allele.outcome", "other_allele.outcome",
                           "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "info.outcome")

#############################
#STEP 4: QC on EXPOSURE DATA#
#############################

#1. Check MAF:

summary(exposure_fiadjbmi$eaf.exposure) #perfect!

#2. Let's check for MHC:

fiadjbmi[which(str_detect(exposure_fiadjbmi$chr_pos, "chr6:") == TRUE),] #none in MHC! No problemo.

############################
#STEP 5: QC on OUTCOME DATA#
############################

#1. Let's check the MAF:

summary(as.numeric(outcome_father$eaf.outcome)) #there is one NA, but since the MAF in exposure is OK we are gonna keept it.

outcome_father$eaf.outcome <- as.numeric(outcome_father$eaf.outcome)/100

summary(as.numeric(outcome_father$eaf.outcome)) #all in order, but wildly different max of frequencies. Maybe due to different alleles???

#2. Let's check the INFO:

summary(outcome_father$info.outcome) #all really high values. PERFECT.

################################################
#STEP 6: ALIGN DATA TO THE INCREASE IR VARIANTS#
################################################

new_a1 <- ifelse(as.numeric(exposure_fiadjbmi$beta.exposure) < 0, exposure_fiadjbmi$other_allele.exposure, exposure_fiadjbmi$effect_allele.exposure)
new_a2 <- ifelse(as.numeric(exposure_fiadjbmi$beta.exposure) < 0, exposure_fiadjbmi$effect_allele.exposure, exposure_fiadjbmi$other_allele.exposure)
new_beta <- ifelse(as.numeric(exposure_fiadjbmi$beta.exposure) < 0, (-1)*(exposure_fiadjbmi$beta.exposure), exposure_fiadjbmi$beta.exposure)
new_eaf <- ifelse(as.numeric(exposure_fiadjbmi$beta.exposure) < 0, 1-exposure_fiadjbmi$eaf.exposure, exposure_fiadjbmi$eaf.exposure)

exposure_fiadjbmi_pos <- exposure_fiadjbmi

exposure_fiadjbmi_pos$effect_allele.exposure <- new_a1
exposure_fiadjbmi_pos$other_allele.exposure <- new_a2
exposure_fiadjbmi_pos$beta.exposure <- new_beta
exposure_fiadjbmi_pos$eaf.exposure <- new_eaf

#Now let's harmonise the data taking that into account:

exposure_fiadjbmi_pos$id.exposure <- "FIadjBMI"
exposure_fiadjbmi_pos$exposure <- "FIadjBMI"

outcome_father$id.outcome <- "Father effect on offspring's birth weight"
outcome_father$outcome <- "Father effect on offspring's birth weight"

dat_1_pre <- harmonise_data(exposure_fiadjbmi_pos, outcome_father, action = 3) #15: perfect.

#Two are palindromic..., but what kind?

dat_1_filt <- dat_1_pre[which(dat_1_pre$remove == FALSE),] #31/31
dat_1_filt <- dat_1_filt[which(dat_1_filt$palindromic == FALSE & dat_1_filt$ambiguous == FALSE | dat_1_filt$palindromic == TRUE & dat_1_filt$ambiguous == FALSE),] #31/31

#Let's check whether the effects are IR-increasing:

summary(dat_1_filt$beta.exposure) #YES. All aligneds to the increasing allele.

#Finally, let's remove the variants without SE:

dat_1_filt_clean <- dat_1_filt[which(!(is.na(dat_1_filt$se.outcome) == TRUE)),]

#Let's see which one was removed:

dat_1_filt[which(!(dat_1_filt$SNP%in%dat_1_filt_clean$SNP)),]

saveRDS(dat_1_filt_clean, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FIadjBMI/Father/FIadjBMI_father_merged_df")

#################
#STEP 7: RUN GRS#
#################

#We need an approximation of the samplesizes: 

dat_1_filt_clean$samplesize.outcome <- 60134

dat_1_filt_clean$unweighted_beta <- 1

grs_weighted <- grs.summary(dat_1_filt_clean$beta.exposure, dat_1_filt_clean$beta.outcome, dat_1_filt_clean$se.outcome, dat_1_filt_clean$samplesize.outcome)
grs_unweighted <- grs.summary(dat_1_filt_clean$unweighted_beta, dat_1_filt_clean$beta.outcome, dat_1_filt_clean$se.outcome, dat_1_filt_clean$samplesize.outcome)

saveRDS(grs_weighted, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FIadjBMI/Father/FIadjBMI_father_weighted")
saveRDS(grs_unweighted, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FIadjBMI/Father/FIadjBMI_father_unweighted")

########################################
#STEP 8: Check validity of IVS for 2SMR#
########################################

#We are gonna test the validity by:

#1. Calculating Steiger to check that there is no bidirectionality:

dat_1_filt_clean$units.exposure <- "SD"
dat_1_filt_clean$units.outcome <- "SD"

dat_1_filt_steiger <- steiger_filtering(dat_1_filt_clean) #some NaNs were produced and a variant has been removed. I think one SE was NA, hence the issue.
dat_1_filt_steiger <- dat_1_filt_steiger[which(dat_1_filt_steiger$steiger_dir == TRUE),] #30/31. Perfect!!

dat_1_filt_steiger$mr_keep <- TRUE #this allows us to keep the unambiguous.

saveRDS(dat_1_filt_steiger, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/Father/FIadjBMI_father_merged_df")

dat_1_filt_steiger <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/Father/FIadjBMI_father_merged_df")

#We have 21/31

#Let's check which ones have been removed:

dat_1_filt_clean$SNP[which(!(dat_1_filt_clean$SNP%in%dat_1_filt_steiger$SNP))]

#2. Let's calculate the F-statistic:

F_ = ((dat_1_filt_steiger$beta.exposure)^2)/((dat_1_filt_steiger$se.exposure)^2)

summary(F_) #all of them above 10!

mF  = mean(F_)

print(mF)
#58.06345

Isq(dat_1_filt_steiger$beta.exposure, dat_1_filt_steiger$se.exposure)
#0.8698437

##################
#STEP 9: RUN 2SMR#
##################

mr_results <- mr(dat_1_filt_steiger)

saveRDS(mr_results, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/Father/FIadjBMI_father_original")

###############################
#STEP 10: sensitivity analysis#
###############################

#1: Rucker test:

rucker <- mr_rucker(dat_1_filt_steiger) #Rucker sees a bit of difference, but still chooses IVW random effects over Egger.

saveRDS(rucker, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/Father/FIadjBMI_father_rucker_test_original")

#2: I2 and Cochran's Q with meta package

dat_1_filt_steiger$mr <- dat_1_filt_steiger$beta.outcome/dat_1_filt_steiger$beta.exposure
dat_1_filt_steiger$mr_se <- ((dat_1_filt_steiger$mr*((dat_1_filt_steiger$se.exposure/dat_1_filt_steiger$beta.exposure)^2+(dat_1_filt_steiger$se.outcome/dat_1_filt_steiger$beta.outcome)^2)^0.5)^2)^0.5
het_Q_Isq <- metagen(dat_1_filt_steiger$mr, dat_1_filt_steiger$mr_se)

saveRDS(het_Q_Isq, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/Father/FIadjBMI_father_Q_Isq_original")

#3: Sensitivity plots:

dat_1_filt_steiger$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/sensitivity_plots/IR/FIadjBMI_father_original.tiff ", units="in", width=10, height=10, res=300)
mr_plots(dat_1_filt_steiger)
dev.off()

#############################
#STEP 10: DETECTING OUTLIERS#
#############################

library(RadialMR)

radial_input <- format_radial(dat_1_filt_steiger$beta.exposure, dat_1_filt_steiger$beta.outcome, dat_1_filt_steiger$se.exposure, dat_1_filt_steiger$se.outcome, dat_1_filt_steiger$SNP)

radial_output <- ivw_radial(radial_input, alpha = 0.05, weights = 3) 

outliers <- radial_output$outliers$SNP 

#"rs116141873"

dat_1_post <- dat_1_filt_steiger[-which(dat_1_filt_steiger$SNP%in%outliers),]

saveRDS(dat_1_post, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/Father/FIadjBMI_Father_merged_df_post_outliers")

#########################################
#STEP 11: Check validity of IVS for 2SMR#
#########################################

#We have 20/36

#2. Let's calculate the F-statistic:

F_ = ((dat_1_post$beta.exposure)^2)/((dat_1_post$se.exposure)^2)

summary(F_) #all of them above 10!

mF  = mean(F_)

print(mF)
#58.33543

Isq(dat_1_post$beta.exposure, dat_1_post$se.exposure)
#0.8587708

#Not the best Isq I have seen..., that is for sure.
#Check thresholds just in case...

##################
#STEP 9: RUN 2SMR#
##################

mr_results <- mr(dat_1_post)

saveRDS(mr_results, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/Father/FIadjBMI_Father_post")

###############################
#STEP 12: sensitivity analysis#
###############################

#1: Rucker test:

rucker <- mr_rucker(dat_1_post)

saveRDS(rucker, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/Father/FIadjBMI_Father_rucker_test_post")

#2: I2 and Cochran's Q with meta package

dat_1_post$mr <- dat_1_post$beta.outcome/dat_1_post$beta.exposure
dat_1_post$mr_se <- ((dat_1_post$mr*((dat_1_post$se.exposure/dat_1_post$beta.exposure)^2+(dat_1_post$se.outcome/dat_1_post$beta.outcome)^2)^0.5)^2)^0.5
het_Q_Isq <- metagen(dat_1_post$mr, dat_1_post$mr_se)

saveRDS(het_Q_Isq, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/Father/FIadjBMI_Father_Q_Isq_post")

#3: Sensitivity plots:

dat_1_post$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/sensitivity_plots/IR/FIadjBMI_Father_post.tiff ", units="in", width=10, height=10, res=300)
mr_plots(dat_1_post)
dev.off()


###################
#STEP 13: MR-CLUST#
###################

library(mrclust)

bx = dat_1_filt_steiger$beta.exposure
bxse = dat_1_filt_steiger$se.exposure
by = dat_1_filt_steiger$beta.outcome
byse = dat_1_filt_steiger$se.outcome
ratio_est = by/bx
ratio_est_se = byse/abs(bx)

snp_names = dat_1_filt_steiger$SNP

res_em = mr_clust_em(theta = ratio_est, theta_se = ratio_est_se, bx = bx, by = by, bxse = bxse, byse = byse, obs_names = snp_names)

mr_clust_res <-res_em$results$best

saveRDS(mr_clust_res, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/Father/FIadjBMI_father_MR_CLUST")
