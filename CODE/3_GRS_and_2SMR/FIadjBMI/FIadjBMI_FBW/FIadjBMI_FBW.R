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

fiadjbmi <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/FIadjBMI/FIadjBMI_independent_4_GRS.csv")
fbw <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Liftover/OUTPUT/FBW/Birthweight2021_Curated.txt")

fiadjbmi$chr_pos <- paste("chr", fiadjbmi$Chromosome, ":", fiadjbmi$Position_build_37, sep = "")

#######################
#STEP 2: MATCHING DATA#
#######################

#NOTE: the insertion/deletion are in the same format here!! (the position takes into account the amount of deletions, which is great.)

fbw_2_fiadjbmi <- fbw[which(fbw$chr_pos_37%in%fiadjbmi$chr_pos),] #31/31

#We need to change some rsIDs since they are merged and won't be found in the 1000G reference panel.

fiadjbmi$rsID[which(fiadjbmi$rsID == "rs77935490")] <- "rs5017305"

#Thera are two other mismatches: one being an insertion/deletion: rs200678953
#And the other because of the MAF: rs111264094

##########################
#STEP 3: CLEANING COLUMNS#
##########################

#We need to parse the data for 2SMR functions..., 
#so let's get started:

outcome_fbw <- fbw_2_fiadjbmi %>%
  select(rsID, A1,A0, `EGG-frq`, `Beta-A1`, SE, P, `IS-info`)

colnames(outcome_fbw) <- c("SNP", "effect_allele.outcome", "other_allele.outcome",
                           "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "info.outcome")

#Now we do the same for fiadjbmi:

exposure_fiadjbmi <- fiadjbmi %>%
  select(rsID, Effect_Allele, Other_Allele, EAF, Beta, SE, P, N)

colnames(exposure_fiadjbmi) <- c("SNP", "effect_allele.exposure", "other_allele.exposure",
                                 "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "sample_size.exposure")

#############################
#STEP 4: QC on EXPOSURE DATA#
#############################

#1. Check MAF:

summary(exposure_fiadjbmi$eaf.exposure) #perfect!

#2. Let's check for MHC:

fiadjbmi[which(str_detect(exposure_fiadjbmi$SNP, "chr6:") == TRUE),] #none in MHC! No problemo.

############################
#STEP 5: QC on OUTCOME DATA#
############################

#1. Let's check the MAF:

summary(as.numeric(outcome_fbw$eaf.outcome)) #there is one NA, but since the MAF in exposure is OK we are gonna keept it.

outcome_fbw$eaf.outcome <- as.numeric(outcome_fbw$eaf.outcome)/100

summary(as.numeric(outcome_fbw$eaf.outcome)) #all in order, but wildly different max of frequencies. Maybe due to different alleles???

#2. Let's check the INFO:

summary(outcome_fbw$info.outcome) #all really high values. PERFECT.

#3. Also the heterogeneity, though that value is in the original data:

summary(as.numeric(fbw_2_fiadjbmi$I2)) #we have some of them with huge heterogeneity.
summary(as.numeric(fbw_2_fiadjbmi$`P-het`)) #But only very few with significant heterogenety values.

#4. Some of the variants have beta = 0. We will remove them after merging them.

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

outcome_fbw$id.outcome <- "Fetal effect on own Birth Weight"
outcome_fbw$outcome <- "Fetal effect on own Birth Weight"

dat_1_pre <- harmonise_data(exposure_fiadjbmi_pos, outcome_fbw, action = 3) 

#Two are palindromic..., but what kind?

dat_1_filt <- dat_1_pre[which(dat_1_pre$remove == FALSE),] #31/31
dat_1_filt <- dat_1_filt[which(dat_1_filt$palindromic == FALSE & dat_1_filt$ambiguous == FALSE | dat_1_filt$palindromic == TRUE & dat_1_filt$ambiguous == FALSE),] #31/31

#Let's check whether the effects are IR-increasing:

summary(dat_1_filt$beta.exposure) #YES. All aligneds to the increasing allele.

#Finally, let's remove the variants without SE:

dat_1_filt_clean <- dat_1_filt[which(!(is.na(dat_1_filt$se.outcome) == TRUE)),]

#Which one is being removed?

dat_1_filt$SNP[which(!(dat_1_filt$SNP%in%dat_1_filt_clean$SNP))]

saveRDS(dat_1_filt_clean, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FIadjBMI/FBW/FIadjBMI_FBW_merged_df")

#################
#STEP 7: RUN GRS#
#################

#We need an approximation of the samplesizes: 

dat_1_filt_clean$samplesize.outcome <- 423683

dat_1_filt_clean$unweighted_beta <- 1

grs_weighted <- grs.summary(dat_1_filt_clean$beta.exposure, dat_1_filt_clean$beta.outcome, dat_1_filt_clean$se.outcome, dat_1_filt_clean$samplesize.outcome)
grs_unweighted <- grs.summary(dat_1_filt_clean$unweighted_beta, dat_1_filt_clean$beta.outcome, dat_1_filt_clean$se.outcome, dat_1_filt_clean$samplesize.outcome)

saveRDS(grs_weighted, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FIadjBMI/FBW/FIadjBMI_FBW_weighted")
saveRDS(grs_unweighted, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FIadjBMI/FBW/FIadjBMI_FBW_unweighted")

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

saveRDS(dat_1_filt_steiger, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBW/FIadjBMI_FBW_merged_df")

dat_1_filt_steiger <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBW/FIadjBMI_FBW_merged_df")

#Let's check the ones that have been removed:

dat_1_filt_clean$SNP[which(!(dat_1_filt_clean$SNP%in%dat_1_filt_steiger$SNP))]

#We have 24/31

#2. Let's calculate the F-statistic:

F_ = ((dat_1_filt_steiger$beta.exposure)^2)/((dat_1_filt_steiger$se.exposure)^2)

summary(F_) #all of them above 10!

mF  = mean(F_)

print(mF)
#57.4801

Isq(dat_1_filt_steiger$beta.exposure, dat_1_filt_steiger$se.exposure)
#0.8332531

##################
#STEP 9: RUN 2SMR#
##################

mr_results <- mr(dat_1_filt_steiger)

saveRDS(mr_results, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBW/FIadjBMI_FBW_original")

###############################
#STEP 10: sensitivity analysis#
###############################

#1: Rucker test:

rucker <- mr_rucker(dat_1_filt_steiger) #Rucker sees a bit of difference, but still chooses IVW random effects over Egger.

saveRDS(rucker, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBW/FIadjBMI_FBW_rucker_test_original")

#2: I2 and Cochran's Q with meta package

dat_1_filt_steiger$mr <- dat_1_filt_steiger$beta.outcome/dat_1_filt_steiger$beta.exposure
dat_1_filt_steiger$mr_se <- ((dat_1_filt_steiger$mr*((dat_1_filt_steiger$se.exposure/dat_1_filt_steiger$beta.exposure)^2+(dat_1_filt_steiger$se.outcome/dat_1_filt_steiger$beta.outcome)^2)^0.5)^2)^0.5
het_Q_Isq <- metagen(dat_1_filt_steiger$mr, dat_1_filt_steiger$mr_se)

saveRDS(het_Q_Isq, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBW/FIadjBMI_FBW_Q_Isq_original")

#3: Sensitivity plots:

dat_1_filt_steiger$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/sensitivity_plots/IR/FIadjBMI_FBW_original.tiff ", units="in", width=10, height=10, res=300)
mr_plots(dat_1_filt_steiger)
dev.off()

#According to the sensitivity plots we have very strong variants
#associated with increase birthweight 
#and a whole variants that carries the IVW: "rs1260326".

#############################
#STEP 10: DETECTING OUTLIERS#
#############################

library(RadialMR)

radial_input <- format_radial(dat_1_filt_steiger$beta.exposure, dat_1_filt_steiger$beta.outcome, dat_1_filt_steiger$se.exposure, dat_1_filt_steiger$se.outcome, dat_1_filt_steiger$SNP)

radial_output <- ivw_radial(radial_input, alpha = 0.05, weights = 3) 

outliers <- radial_output$outliers$SNP 

#[1] "rs11727676"  "rs118164457" "rs1260326"  
#[4] "rs2943646"   "rs459193"    "rs6487237"  
#[7] "rs7133378"   "rs731839"    "rs860598"  

#The outliers are a mixed bag. Some of them are good hits, but many are
#actually not outliers, but huge hits.

dat_1_post <- dat_1_filt_steiger[-which(dat_1_filt_steiger$SNP%in%outliers),]

saveRDS(dat_1_post, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBW/FIadjBMI_FBW_merged_df_post_outliers")

#########################################
#STEP 11: Check validity of IVS for 2SMR#
#########################################

#We have 15/36

#2. Let's calculate the F-statistic:

F_ = ((dat_1_post$beta.exposure)^2)/((dat_1_post$se.exposure)^2)

summary(F_) #all of them above 10!

mF  = mean(F_)

print(mF)
#50.29824

Isq(dat_1_post$beta.exposure, dat_1_post$se.exposure)
#0.8314253

#Not the best Isq I have seen..., that is for sure.
#Check thresholds just in case...

##################
#STEP 9: RUN 2SMR#
##################

mr_results <- mr(dat_1_post)

saveRDS(mr_results, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBW/FIadjBMI_FBW_post")

###############################
#STEP 12: sensitivity analysis#
###############################

#1: Rucker test:

rucker <- mr_rucker(dat_1_post)

saveRDS(rucker, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBW/FIadjBMI_FBW_rucker_test_post")

#2: I2 and Cochran's Q with meta package

dat_1_post$mr <- dat_1_post$beta.outcome/dat_1_post$beta.exposure
dat_1_post$mr_se <- ((dat_1_post$mr*((dat_1_post$se.exposure/dat_1_post$beta.exposure)^2+(dat_1_post$se.outcome/dat_1_post$beta.outcome)^2)^0.5)^2)^0.5
het_Q_Isq <- metagen(dat_1_post$mr, dat_1_post$mr_se)

saveRDS(het_Q_Isq, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBW/FIadjBMI_FBW_Q_Isq_post")

#3: Sensitivity plots:

dat_1_post$labels <- NA

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/sensitivity_plots/IR/FIadjBMI_FBW_post.tiff ", units="in", width=10, height=10, res=300)
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

saveRDS(mr_clust_res, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/FIadjBMI/FBW/FIadjBMI_FBW_MR_CLUST")
