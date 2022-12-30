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
library(jsonlite)
library(httr)

#########################
#Making space for memory#
#########################

memory.limit(size=8000000)

###################
#LOADING FUNCTIONS#
###################

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

fiadjbmi <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/FIadjBMI_independent.txt")
fbw <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/FBWAdj/FBWadj_curated.txt")

fiadjbmi$chr_pos <- paste("chr", fiadjbmi$Chromosome, ":", fiadjbmi$Position_build_37, sep = "")

#######################
#STEP 2: MATCHING DATA#
#######################

#NOTE: the insertion/deletion are in the same format here!! (the position takes into account the amount of deletions, which is great.)

fbw_2_fiadjbmi <- fbw[which(fbw$chr_pos%in%fiadjbmi$chr_pos),] #32/36

#Let's take outliers for those missing:

missing_snps <- fiadjbmi[which(!(fiadjbmi$chr_pos%in%fbw_2_fiadjbmi$chr_pos)),]

proxies_df <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/3_2SMR_and_GRS/FIadjBMI/proxies/combined_query_snp_list.txt")

query_snps <- proxies_df[which(proxies_df$Distance == 0),] #only 2 queries

proxies_df_clean <- proxies_df[which(!(proxies_df$Distance == 0)),]

proxies_df_clean <- proxies_df_clean[which(proxies_df_clean$R2 > 0.8),]

#######################################
#STEP 3: get the proxies found in both#
#######################################

FIadjBMI_original <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Romain_Mario/RAW_DATA/FIAdjBMI_2021/MAGIC1000G_FI_EUR.tsv.gz")
FIadjBMI_original$chr_pos <- paste("chr",FIadjBMI_original$chromosome, ":",FIadjBMI_original$base_pair_location, sep = "")

FIadjBMI_proxies <- FIadjBMI_original[which(FIadjBMI_original$chr_pos%in%proxies_df_clean$Coord),] #304
fbw_2_proxies <- fbw[which(fbw$chr_pos%in%proxies_df_clean$Coord),] #247

#Let's clean the MAF for both sets:

FIadjBMI_proxies <- FIadjBMI_proxies[which(FIadjBMI_proxies$effect_allele_frequency > 0.01 & FIadjBMI_proxies$effect_allele_frequency < 0.99),]
fbw_2_proxies <- fbw_2_proxies[which(fbw_2_proxies$eaf > 0.01 & fbw_2_proxies$eaf < 0.99),] #247

FIadjBMI_proxies_match <- FIadjBMI_proxies[which(FIadjBMI_proxies$chr_pos%in%fbw_2_proxies$chr_pos),] #perfect
fbw_2_proxies_match <- fbw_2_proxies[which(fbw_2_proxies$chr_pos%in%FIadjBMI_proxies_match$chr_pos),]

#Let's get the rsIDs first:

FIadjBMI_proxies_match <- FIadjBMI_proxies_match[order(match(FIadjBMI_proxies_match$chr_pos, fbw_2_proxies_match$chr_pos)),]

length(which(FIadjBMI_proxies_match$chr_pos == fbw_2_proxies_match$chr_pos)) #all match

FIadjBMI_proxies_match$rsid <- fbw_2_proxies_match$RSID

#The all match perfectly!!
#Let's clump the data. To do so we need the pre-proxies (the 32 matching)
#and the proxies for the missing ones.
#In the best case scenario we should get 35 variants.

FIadjBMI_pre_proxies <- fiadjbmi[which(!(fiadjbmi$chr_pos%in%missing_snps$chr_pos)),]

#AWESOME, let's get the data necessary to clump these fellas.

FIadjBMI_pre_proxies_4_clump <- FIadjBMI_pre_proxies %>%
  select(rsID, P)

FIadjBMI_proxies_match_4_clump <- FIadjBMI_proxies_match %>%
  select(rsid, p_value)

colnames(FIadjBMI_pre_proxies_4_clump) <- c("rsid", "pval")
colnames(FIadjBMI_proxies_match_4_clump) <- c("rsid", "pval")

FIadjBMI_4_clumping <- rbind(FIadjBMI_pre_proxies_4_clump, FIadjBMI_proxies_match_4_clump)

FIadjBMI_independent <- ieugwasr::ld_clump_local(FIadjBMI_4_clumping, bfile = "C:/Users/zlc436/Downloads/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/win-library/4.0/plinkbinr/bin/plink_Windows.exe", clump_kb = 500, clump_r2 = 0.01, clump_p = 0.99) 

#####################################
#Let's get the data set for these 29#
#####################################

outcome_fbw <- fbw[which(fbw$RSID%in%FIadjBMI_independent$rsid),] #29/29
exposure_fiadjbmi <- FIadjBMI_original[which(FIadjBMI_original$chr_pos%in%outcome_fbw$chr_pos),]

#PERFECT: 29 matches.

##########################
#STEP 3: CLEANING COLUMNS#
##########################

#We need to parse the data for 2SMR functions..., 
#so let's get started:

outcome_fbw <- outcome_fbw %>%
  select(chr_pos, ea,nea, eaf, beta, se, p, n_ownBW)

colnames(outcome_fbw) <- c("SNP", "effect_allele.outcome", "other_allele.outcome",
                           "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome")

#Now we do the same for fiadjbmi:

exposure_fiadjbmi <- exposure_fiadjbmi %>%
  select(chr_pos, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size)

colnames(exposure_fiadjbmi) <- c("SNP", "effect_allele.exposure", "other_allele.exposure",
                                 "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "sample_size.exposure")


#############################
#STEP 4: QC on EXPOSURE DATA#
#############################

#1. Check MAF:

summary(exposure_fiadjbmi$eaf.exposure) #perfect!

#2. Let's check for MHC:

exposure_fiadjbmi[which(str_detect(exposure_fiadjbmi$SNP, "chr6:") == TRUE),] #none in MHC! No problemo.

############################
#STEP 5: QC on OUTCOME DATA#
############################

#1. Let's check the MAF:

summary(as.numeric(outcome_fbw$eaf.outcome)) 

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

outcome_fbw$id.outcome <- "Fetal effect on Birth Weight adjusted for Maternal Effect"
outcome_fbw$outcome <- "Fetal effect on Birth Weight adjusted for Maternal Effect"

dat_1_pre <- harmonise_data(exposure_fiadjbmi_pos, outcome_fbw, action = 3) #29: perfect.

#Two are palindromic..., but what kind?

dat_1_filt <- dat_1_pre[which(dat_1_pre$remove == FALSE),] #29/29
dat_1_filt <- dat_1_filt[which(dat_1_filt$palindromic == FALSE & dat_1_filt$ambiguous == FALSE | dat_1_filt$palindromic == TRUE & dat_1_filt$ambiguous == FALSE),] #15/36

#Let's check whether the effects are IR-increasing:

summary(dat_1_filt$beta.exposure) #YES. All aligneds to the increasing allele.

#Finally, let's remove the variants without SE:

dat_1_filt_clean <- dat_1_filt[which(!(is.na(dat_1_filt$se.outcome) == TRUE)),] #29/29

saveRDS(dat_1_filt_clean, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_merged_df")

#################
#STEP 7: RUN GRS#
#################

#We need an approximation of the samplesizes: 

dat_1_filt_clean$unweighted_beta <- 1

grs_weighted <- grs.summary(dat_1_filt_clean$beta.exposure, dat_1_filt_clean$beta.outcome, dat_1_filt_clean$se.outcome, dat_1_filt_clean$samplesize.outcome)
grs_unweighted <- grs.summary(dat_1_filt_clean$unweighted_beta, dat_1_filt_clean$beta.outcome, dat_1_filt_clean$se.outcome, dat_1_filt_clean$samplesize.outcome)

saveRDS(grs_weighted, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_weighted")
saveRDS(grs_unweighted, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_unweighted")

########################################
#STEP 8: Check validity of IVS for 2SMR#
########################################

#We are gonna test the validity by:

#1. Calculating Steiger to check that there is no bidirectionality:

dat_1_filt_clean$units.exposure <- "SD"
dat_1_filt_clean$units.outcome <- "SD"

dat_1_filt_steiger <- steiger_filtering_old(dat_1_filt_clean) 
dat_1_filt_steiger <- dat_1_filt_steiger[which(dat_1_filt_steiger$steiger_dir == TRUE),] #27/36. Not bad!!

dat_1_filt_steiger$mr_keep <- TRUE #this allows us to keep the unambiguous.

saveRDS(dat_1_filt_steiger, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_merged_df")

dat_1_filt_steiger <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_merged_df")

#We have 25/29

#2. Let's calculate the F-statistic:

F_ = ((dat_1_filt_steiger$beta.exposure)^2)/((dat_1_filt_steiger$se.exposure)^2)

summary(F_) #all of them above 10!

mF  = mean(F_)

print(mF)
#56.15594

Isq(dat_1_filt_steiger$beta.exposure, dat_1_filt_steiger$se.exposure)
#0.8266854

##################
#STEP 9: RUN 2SMR#
##################

mr_results <- mr(dat_1_filt_steiger)

saveRDS(mr_results, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_original")

###############################
#STEP 10: sensitivity analysis#
###############################

#1: Rucker test:

rucker <- mr_rucker(dat_1_filt_steiger) #Rucker sees a bit of difference, but still chooses IVW random effects over Egger.

saveRDS(rucker, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_rucker_test_original")

#2: I2 and Cochran's Q with meta package

dat_1_filt_steiger$mr <- dat_1_filt_steiger$beta.outcome/dat_1_filt_steiger$beta.exposure
dat_1_filt_steiger$mr_se <- ((dat_1_filt_steiger$mr*((dat_1_filt_steiger$se.exposure/dat_1_filt_steiger$beta.exposure)^2+(dat_1_filt_steiger$se.outcome/dat_1_filt_steiger$beta.outcome)^2)^0.5)^2)^0.5
het_Q_Isq <- metagen(dat_1_filt_steiger$mr, dat_1_filt_steiger$mr_se)

saveRDS(het_Q_Isq, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_Q_Isq_original")

#3: Sensitivity plots:

dat_1_filt_steiger$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_original.tiff ", units="in", width=10, height=10, res=300)
mr_plots(dat_1_filt_steiger)
dev.off()

#When we remove two variants they assocation becomes significant.
#They are bound to be outliers.
#funnel plot also shows asymmetry towards the negative side.
#There are bonafide outliers that are not good for MR.

#############################
#STEP 10: DETECTING OUTLIERS#
#############################

library(RadialMR)

radial_input <- format_radial(dat_1_filt_steiger$beta.exposure, dat_1_filt_steiger$beta.outcome, dat_1_filt_steiger$se.exposure, dat_1_filt_steiger$se.outcome, dat_1_filt_steiger$SNP)

radial_output <- ivw_radial(radial_input, alpha = 0.05, weights = 3) 

outliers <- radial_output$outliers$SNP 

#"chr19:33899065" 
#"chr2:227099534" 
#"chr6:34222201" 

#The outliers are a mixed bag. Some of them are good hits, but many are
#actually not outliers, but huge hits.

dat_1_post <- dat_1_filt_steiger[-which(dat_1_filt_steiger$SNP%in%outliers),]

saveRDS(dat_1_post, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_merged_df_post_outliers")

#########################################
#STEP 11: Check validity of IVS for 2SMR#
#########################################

#We have 25/36

#2. Let's calculate the F-statistic:

F_ = ((dat_1_post$beta.exposure)^2)/((dat_1_post$se.exposure)^2)

summary(F_) #all of them above 10!

mF  = mean(F_)

print(mF)
#51.90189

Isq(dat_1_post$beta.exposure, dat_1_post$se.exposure)
#0.7684725

#Not the best Isq I have seen..., that is for sure.
#Check thresholds just in case...

##################
#STEP 9: RUN 2SMR#
##################

mr_results <- mr(dat_1_post)

saveRDS(mr_results, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_post")

###############################
#STEP 12: sensitivity analysis#
###############################

#1: Rucker test:

rucker <- mr_rucker(dat_1_post)

saveRDS(rucker, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_rucker_test_post")

#2: I2 and Cochran's Q with meta package

dat_1_post$mr <- dat_1_post$beta.outcome/dat_1_post$beta.exposure
dat_1_post$mr_se <- ((dat_1_post$mr*((dat_1_post$se.exposure/dat_1_post$beta.exposure)^2+(dat_1_post$se.outcome/dat_1_post$beta.outcome)^2)^0.5)^2)^0.5
het_Q_Isq <- metagen(dat_1_post$mr, dat_1_post$mr_se)

saveRDS(het_Q_Isq, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_Q_Isq_post")

#3: Sensitivity plots:

dat_1_post$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_post.tiff ", units="in", width=10, height=10, res=300)
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

#According the MR-clust we have just one cluster!!!
#The others are actually NULLs!!!

#Let's plot these fellas:

res_em$plots

#Which are the real ones??

cluster_1 <- res_em$results$best[which(res_em$results$best$cluster_class == 1),] #7 variants associated

summary(cluster_1$probability) #not all of them have a superhigh probablity, but I guess it is fine.

dat_1_cluster_1 <- dat_1_filt_steiger[which(dat_1_filt_steiger$SNP%in%cluster_1$observation),]

#########################################
#STEP 14: Check validity of IVS for 2SMR#
#########################################

#2. Let's calculate the F-statistic:

F_ = ((dat_1_cluster_1$beta.exposure)^2)/((dat_1_cluster_1$se.exposure)^2)

summary(F_) #all of them above 10!

mF  = mean(F_)

print(mF)
#70.28662

Isq(dat_1_cluster_1$beta.exposure, dat_1_cluster_1$se.exposure)
#0.8512054

#Not the best Isq I have seen..., that is for sure.
#Check thresholds just in case...

###################
#STEP 15: RUN 2SMR#
###################

mr_results <- mr(dat_1_cluster_1) #results are really good.

saveRDS(mr_results, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_cluster_1")

###############################
#STEP 16: sensitivity analysis#
###############################

#1: Rucker test:

rucker <- mr_rucker(dat_1_cluster_1)

saveRDS(rucker, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_rucker_test_clust_1")

#2: I2 and Cochran's Q with meta package

dat_1_cluster_1$mr <- dat_1_cluster_1$beta.outcome/dat_1_cluster_1$beta.exposure
dat_1_cluster_1$mr_se <- ((dat_1_cluster_1$mr*((dat_1_cluster_1$se.exposure/dat_1_cluster_1$beta.exposure)^2+(dat_1_cluster_1$se.outcome/dat_1_cluster_1$beta.outcome)^2)^0.5)^2)^0.5
het_Q_Isq <- metagen(dat_1_cluster_1$mr, dat_1_cluster_1$mr_se)

saveRDS(het_Q_Isq, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_Q_Isq_clust_1")

#3: Sensitivity plots:

dat_1_cluster_1$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBWadj/FIadjBMI_FBWadj_clust_1.tiff ", units="in", width=10, height=10, res=300)
mr_plots(dat_1_cluster_1)
dev.off()

