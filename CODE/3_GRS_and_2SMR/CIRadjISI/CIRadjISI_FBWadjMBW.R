##############
#INTRODUCTION#
##############

#This code will serve to parse the data of the full results of the European hits for 
#the glycemic traits found in supplementary table 2.

###################
#LOADING LIBRARIES#
###################

library(tidyverse)
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

CIRadjISI <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/CIRadjISI/CLEAN_8_CIRadjISI_GRS_SNPs_Manuscript.csv")
fbwadjmbw <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_GW_analysis/FBWAdj/FBWAdj_curated.txt")

#######################
#STEP 2: MATCHING DATA#
#######################

#NOTE: the insertion/deletion are in the same format here!! (the position takes into account the amount of deletions, which is great.)

fbwadjmbw_2_CIRadjISI <- fbwadjmbw[which(fbwadjmbw$RSID%in%CIRadjISI$snp),] #8/8

#Since we know that they are independent from each other, there is no need to do anything else.

##########################
#STEP 3: CLEANING COLUMNS#
##########################

#We need to parse the data for 2SMR functions..., 
#so let's get started:

outcome_fbwadjmbw <- fbwadjmbw_2_CIRadjISI %>%
  select(RSID, ea,nea, eaf, beta, se, p, n_ownBW)

colnames(outcome_fbwadjmbw) <- c("SNP", "effect_allele.outcome", "other_allele.outcome",
                           "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome")

#Now we do the same for CIRadjISI:

exposure_CIRadjISI <- CIRadjISI %>%
  select(snp, Final_A1, Final_A2, Final_EAF, Final_BETA, stderr, pvalue)

exposure_CIRadjISI$N <- 4789

colnames(exposure_CIRadjISI) <- c("SNP", "effect_allele.exposure", "other_allele.exposure",
                                 "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "sample_size.exposure")

#############################
#STEP 4: QC on EXPOSURE DATA#
#############################

#1. Check MAF:

summary(exposure_CIRadjISI$eaf.exposure) #perfect!

#2. We know none of them are in the MHC!

############################
#STEP 5: QC on OUTCOME DATA#
############################

#1. Let's check the MAF:

summary(as.numeric(outcome_fbwadjmbw$eaf.outcome)) #there is one NA, but since the MAF in exposure is OK we are gonna keept it.

################################################
#STEP 6: ALIGN DATA TO THE INCREASE IR VARIANTS#
################################################

new_a1 <- ifelse(as.numeric(exposure_CIRadjISI$beta.exposure) < 0, exposure_CIRadjISI$other_allele.exposure, exposure_CIRadjISI$effect_allele.exposure)
new_a2 <- ifelse(as.numeric(exposure_CIRadjISI$beta.exposure) < 0, exposure_CIRadjISI$effect_allele.exposure, exposure_CIRadjISI$other_allele.exposure)
new_beta <- ifelse(as.numeric(exposure_CIRadjISI$beta.exposure) < 0, (-1)*(exposure_CIRadjISI$beta.exposure), exposure_CIRadjISI$beta.exposure)
new_eaf <- ifelse(as.numeric(exposure_CIRadjISI$beta.exposure) < 0, 1-exposure_CIRadjISI$eaf.exposure, exposure_CIRadjISI$eaf.exposure)

exposure_CIRadjISI_pos <- exposure_CIRadjISI

exposure_CIRadjISI_pos$effect_allele.exposure <- new_a1
exposure_CIRadjISI_pos$other_allele.exposure <- new_a2
exposure_CIRadjISI_pos$beta.exposure <- new_beta
exposure_CIRadjISI_pos$eaf.exposure <- new_eaf

#Now let's harmonise the data taking that into account:

exposure_CIRadjISI_pos$id.exposure <- "CIRadjISI"
exposure_CIRadjISI_pos$exposure <- "CIRadjISI"

outcome_fbwadjmbw$id.outcome <- "Fetal Own Birth Weight adjusted for maternal effect"
outcome_fbwadjmbw$outcome <- "Fetal Own Birth Weight adjusted for maternal effect"

dat_1_pre <- harmonise_data(exposure_CIRadjISI_pos, outcome_fbwadjmbw, action = 3) #15: perfect.

#Two are palindromic..., but what kind?

dat_1_filt <- dat_1_pre[which(dat_1_pre$remove == FALSE),] #8/8
dat_1_filt <- dat_1_filt[which(dat_1_filt$palindromic == FALSE & dat_1_filt$ambiguous == FALSE | dat_1_filt$palindromic == TRUE & dat_1_filt$ambiguous == FALSE),] #8/8

#Let's check whether the effects are IR-increasing:

summary(dat_1_filt$beta.exposure) #YES. All aligneds to the increasing allele.

#Finally, let's remove the variants without SE:

dat_1_filt_clean <- dat_1_filt[which(!(is.na(dat_1_filt$se.outcome) == TRUE)),]

#Which one is removed?

dat_1_filt[which(is.na(dat_1_filt$se.outcome) == TRUE),]

saveRDS(dat_1_filt_clean, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/CIRadjISI/FBWadjMBW/CIRadjISI_FBWadjMBW_merged_df")

#################
#STEP 7: RUN GRS#
#################

dat_1_filt_clean$unweighted_beta <- 1

grs_weighted <- grs.summary(dat_1_filt_clean$beta.exposure, dat_1_filt_clean$beta.outcome, dat_1_filt_clean$se.outcome, dat_1_filt_clean$samplesize.outcome)
grs_unweighted <- grs.summary(dat_1_filt_clean$unweighted_beta, dat_1_filt_clean$beta.outcome, dat_1_filt_clean$se.outcome, dat_1_filt_clean$samplesize.outcome)

saveRDS(grs_weighted, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/CIRadjISI/FBWadjMBW/CIRadjISI_FBWadjMBW_weighted")
saveRDS(grs_unweighted, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/CIRadjISI/FBWadjMBW/CIRadjISI_FBWadjMBW_unweighted")

########################################
#STEP 8: Check validity of IVS for 2SMR#
########################################

#We are gonna test the validity by:

#1. Calculating Steiger to check that there is no bidirectionality:

dat_1_filt_clean$units.exposure <- "SD"
dat_1_filt_clean$units.outcome <- "SD"

dat_1_filt_steiger <- steiger_filtering_old(dat_1_filt_clean) 
dat_1_filt_steiger <- dat_1_filt_steiger[which(dat_1_filt_steiger$steiger_dir == TRUE),] #7/7. Not bad!!

dat_1_filt_steiger$mr_keep <- TRUE #this allows us to keep the unambiguous.

saveRDS(dat_1_filt_steiger, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/CIRadjISI/FBWadjMBW/CIRadjISI_fbwadjmbw_merged_df")

dat_1_filt_steiger <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/CIRadjISI/FBWadjMBW/CIRadjISI_fbwadjmbw_merged_df")

#We have 8/8

#2. Let's calculate the F-statistic:

F_ = ((dat_1_filt_steiger$beta.exposure)^2)/((dat_1_filt_steiger$se.exposure)^2)

summary(F_) #one below 10!

#Which one?

dat_1_filt_steiger[which(F_ < 10),]

dat_1_filt_steiger <- dat_1_filt_steiger[which(F_ > 10),]

F_ = ((dat_1_filt_steiger$beta.exposure)^2)/((dat_1_filt_steiger$se.exposure)^2)

summary(F_) #one below 10!

mF  = mean(F_)

print(mF)
#25.1928

Isq(dat_1_filt_steiger$beta.exposure, dat_1_filt_steiger$se.exposure)
#0.5646962

##################
#STEP 9: RUN 2SMR#
##################

mr_results <- mr(dat_1_filt_steiger)

saveRDS(mr_results, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/CIRadjISI/FBWadjMBW/CIRadjISI_fbwadjmbw_original")

#There are mixed-results.

###############################
#STEP 10: sensitivity analysis#
###############################

#1: Rucker test:

rucker <- mr_rucker(dat_1_filt_steiger) #Rucker sees a bit of difference, but still chooses IVW random effects over Egger.

#Huge pleiotropy, but it seems that is homogeneous across the variants. IVW with random effects is better.

saveRDS(rucker, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/CIRadjISI/FBWadjMBW/CIRadjISI_fbwadjmbw_rucker_test_original")

#2: I2 and Cochran's Q with meta package

dat_1_filt_steiger$mr <- dat_1_filt_steiger$beta.outcome/dat_1_filt_steiger$beta.exposure
dat_1_filt_steiger$mr_se <- ((dat_1_filt_steiger$mr*((dat_1_filt_steiger$se.exposure/dat_1_filt_steiger$beta.exposure)^2+(dat_1_filt_steiger$se.outcome/dat_1_filt_steiger$beta.outcome)^2)^0.5)^2)^0.5
het_Q_Isq <- metagen(dat_1_filt_steiger$mr, dat_1_filt_steiger$mr_se)

saveRDS(het_Q_Isq, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/CIRadjISI/FBWadjMBW/CIRadjISI_fbwadjmbw_Q_Isq_original")

#3: Sensitivity plots:

dat_1_filt_steiger$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/sensitivity_plots/IS/CIRadjISI_fbwadjmbw_original.tiff ", units="in", width=10, height=10, res=300)
mr_plots(dat_1_filt_steiger)
dev.off()

#############################
#STEP 10: DETECTING OUTLIERS#
#############################

library(RadialMR)

radial_input <- format_radial(dat_1_filt_steiger$beta.exposure, dat_1_filt_steiger$beta.outcome, dat_1_filt_steiger$se.exposure, dat_1_filt_steiger$se.outcome, dat_1_filt_steiger$SNP)

radial_output <- ivw_radial(radial_input, alpha = 0.05, weights = 3) 

outliers <- radial_output$outliers$SNP 

#"rs10830963"
#"rs12549902"
#"rs4502156"
#"rs7756992"
#"rs7923866" 
#"rs933360"

#The outliers are a mixed bag. Some of them are good hits, but many are
#actually not outliers, but huge hits.

dat_1_post <- dat_1_filt_steiger[-which(dat_1_filt_steiger$SNP%in%outliers),]

saveRDS(dat_1_post, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/CIRadjISI/FBWadjMBW/CIRadjISI_fbwadjmbw_merged_df_post_outliers")

#########################################
#STEP 11: Check validity of IVS for 2SMR#
#########################################

#We have 14/36

#2. Let's calculate the F-statistic:

F_ = ((dat_1_post$beta.exposure)^2)/((dat_1_post$se.exposure)^2)

summary(F_) #all of them above 10!

mF  = mean(F_)

print(mF)
#15.34028

Isq(dat_1_post$beta.exposure, dat_1_post$se.exposure)
#1

#Not the best Isq I have seen..., that is for sure.
#Check thresholds just in case...

##################
#STEP 9: RUN 2SMR#
##################

mr_results <- mr(dat_1_post)

saveRDS(mr_results, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/CIRadjISI/FBWadjMBW/CIRadjISI_fbwadjmbw_post")

###############################
#STEP 12: sensitivity analysis#
###############################

#1: Rucker test:

#rucker <- mr_rucker(dat_1_post)
#saveRDS(rucker, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/CIRadjISI/father/CIRadjISI_father_rucker_test_post")

#2: I2 and Cochran's Q with meta package

dat_1_post$mr <- dat_1_post$beta.outcome/dat_1_post$beta.exposure
dat_1_post$mr_se <- ((dat_1_post$mr*((dat_1_post$se.exposure/dat_1_post$beta.exposure)^2+(dat_1_post$se.outcome/dat_1_post$beta.outcome)^2)^0.5)^2)^0.5
het_Q_Isq <- metagen(dat_1_post$mr, dat_1_post$mr_se) #there is no heterogeneity here my friends :) -> just one SNP

#saveRDS(het_Q_Isq, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/CIRadjISI/father/CIRadjISI_father_Q_Isq_post")

#3: Sensitivity plots:

dat_1_post$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/sensitivity_plots/CIRadjISI_fbwadjmbw_post.tiff ", units="in", width=10, height=10, res=300)
mr_plots(dat_1_post)
dev.off()

#CANNOT BE RUN!

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

saveRDS(mr_clust_res, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/CIRadjISI/FBWadjMBW/CIRadjISI_FBWadjMBW_MR_CLUST")
