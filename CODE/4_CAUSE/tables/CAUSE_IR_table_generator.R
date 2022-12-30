##############
#INTRODUCTION#
##############

#This code is to generate in a comprehensive manner the results for IR GC

###########
#Libraries#
###########

library(data.table)
library(tidyverse)
library(readRDS)

###################
#Loading functions#
###################

elpd_parser <- function(outcome_elpd, outcome_name){
  #We are going to parse the ELPD data for each outcome (combining the three exposures.)
  
  #We need a special dataframe to add things:
  
  seed_elpd <- outcome_elpd[-seq(1,9),]
  
  seed_elpd <- rbind(seed_elpd, as.data.frame(t(colnames(seed_elpd))))
  
  #Important to make sure that we make sure to show which outcome we are using:
  
  colnames(seed_elpd) <- rep(outcome_name,5)
  
  colnames(outcome_elpd) <-  rep(outcome_name,5)
  
  #And we are going to merge it
  
  outcome_elpd <- rbind(seed_elpd, outcome_elpd)
  
  return(outcome_elpd)
  
}

cause_res_parser <- function(exp_1, exp_2, exp_3, outcome_name){
  #This function takes the results for three exposures with the same outcome
  #And merges and cleans the dataframes.
  
  #STEP 1: we need a seed dataframe so that we can save things accordingly.
  
  library(cause)
  
  seed_cause <- summary(exp_1)$tab[-seq(1,2),]
  
  seed_cause <- rbind(seed_cause, as.data.frame(t(colnames(seed_cause))))
  
  seed_cause$V5 <- "p"
  
  colnames(seed_cause) <- rep(outcome_name, 5)
  
  #Now we are ready to merge the data and clean it.
  
  #First we take the results for the models
  
  res_df_1 <- summary(exp_1)$tab
  res_df_2 <- summary(exp_2)$tab
  res_df_3 <- summary(exp_3)$tab
  
  #Then the shared Vs causal model p-value
  
  p_df_1 <- as.data.frame(t(as.data.frame(summary(exp_1)$p)))
  p_df_2 <- as.data.frame(t(as.data.frame(summary(exp_2)$p)))
  p_df_3 <- as.data.frame(t(as.data.frame(summary(exp_3)$p)))
  
  #We do some formating in the pval so that we can merge the data...
  
  p_df_1 <- rbind(p_df_1, p_df_1)
  p_df_2 <- rbind(p_df_2, p_df_2)
  p_df_3 <- rbind(p_df_3, p_df_3)
  
  colnames(p_df_1) <- ("p")
  colnames(p_df_2) <- ("p")
  colnames(p_df_3) <- ("p")
  
  #Now we combine the same-type dataframes, to make things easier
  
  res_df_end <- rbind(res_df_1, res_df_2, res_df_3)
  p_df_end <- rbind(p_df_1, p_df_2, p_df_3)
  
  #And now we combine it all:
  
  final_df <- cbind(res_df_end, p_df_end)
  
  #Finally, we change colnames so that we can merge with the seed:
  
  colnames(final_df) <- rep(outcome_name, 5)
  
  final_df_end <- rbind(seed_cause, final_df)
  
  #And we are ready to send this off:
  
  return(final_df_end)
  
}

###############
#Loading files#
###############

#Now for IR:

FIadjBMI_fbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/5_CAUSE/IR/res_fiadjbmi_FBW.rds")
FIadjBMI_mbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/5_CAUSE/IR/res_fiadjbmi_MBW.rds")
FIadjBMI_fbwadjmbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/5_CAUSE/IR/res_fiadjbmi_FBWadjMBW.rds")
FIadjBMI_mbwadjfbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/5_CAUSE/IR/res_fiadjbmi_MBWadjFBW.rds")
FIadjBMI_father <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/5_CAUSE/IR/res_fiadjbmi_Father.rds")

FI_fbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/5_CAUSE/IR/res_fi_FBW.rds")
FI_mbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/5_CAUSE/IR/res_fi_MBW.rds")
FI_fbwadjmbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/5_CAUSE/IR/res_fi_FBWadjMBW.rds")
FI_mbwadjfbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/5_CAUSE/IR/res_fi_MBWadjFBW.rds")
FI_father <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/5_CAUSE/IR/res_fi_Father.rds")

homair_fbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/5_CAUSE/IR/res_homair_FBW.rds")
homair_mbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/5_CAUSE/IR/res_homair_MBW.rds")
homair_fbwadjmbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/5_CAUSE/IR/res_homair_FBWadjMBW.rds")
homair_mbwadjfbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/5_CAUSE/IR/res_homair_MBWadjFBW.rds")
homair_father <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/5_CAUSE/IR/res_homair_Father.rds")

##############################################
#Let's parse the ELDPs and make them in table#
##############################################

#This one is actually gonna be relatively simple.
#First let's make dataframes for each outcome:

fbw_elpd <- rbind(FIadjBMI_fbw$elpd, FI_fbw$elpd, homair_fbw$elpd)
fbwadjmbw_elpd <- rbind(FIadjBMI_fbwadjmbw$elpd, FI_fbwadjmbw$elpd, homair_fbwadjmbw$elpd)
mbw_elpd <- rbind(FIadjBMI_mbw$elpd, FI_mbw$elpd, homair_mbw$elpd)
mbwadjfbw_elpd <- rbind(FIadjBMI_mbwadjfbw$elpd, FI_mbwadjfbw$elpd, homair_mbwadjfbw$elpd)
father_elpd <- rbind(FIadjBMI_father$elpd, FI_father$elpd, homair_father$elpd)

#Now let's parse it beautifully:

fbw_elpd_end <- elpd_parser(fbw_elpd, "FBW")
fbwadjmbw_elpd_end <- elpd_parser(fbwadjmbw_elpd, "FBWadjMBW")
mbw_elpd_end <- elpd_parser(mbw_elpd, "MBW")
mbwadjfbw_elpd_end <- elpd_parser(mbwadjfbw_elpd, "MBWadjFBW")
father_elpd_end <- elpd_parser(father_elpd, "Father")

#We need a special dataframe to add things:

exposure_column <- c("value", "FIadjBMI", "FIadjBMI", "FIadjBMI", "FI", "FI", "FI", "HOMA-IR", "HOMA-IR", "HOMA-IR")

final_elpd_end <- cbind(exposure_column, fbw_elpd_end, fbwadjmbw_elpd_end, mbw_elpd_end, mbwadjfbw_elpd_end, father_elpd_end)

#And this is done!!

fwrite(final_elpd_end, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/CAUSE/raw/ELPD_IR.csv", sep = ",", col.names = TRUE)

################################################
#Let's parse the results and make them in table#
################################################

library(cause)

fbw_cause <- cause_res_parser(FIadjBMI_fbw, FI_fbw, homair_fbw, "FBW")
fbwadjmbw_cause <- cause_res_parser(FIadjBMI_fbwadjmbw, FI_fbwadjmbw, homair_fbwadjmbw, "FBWadjMBW")
mbw_cause <- cause_res_parser(FIadjBMI_mbw, FI_mbw, homair_mbw, "MBW")
mbwadjfbw_cause <- cause_res_parser(FIadjBMI_mbwadjfbw, FI_mbwadjfbw, homair_mbwadjfbw, "MBWadjFBW")
father_cause <- cause_res_parser(FIadjBMI_father, FI_father, homair_father, "Father")

#Now let's make the columns for the exposure and make the final dataframe:

exposure_column <- c("value", "FIadjBMI", "FIadjBMI", "FI", "FI", "HOMA-IR", "HOMA-IR")

final_cause_end <- cbind(exposure_column, fbw_cause, fbwadjmbw_cause, mbw_cause, mbwadjfbw_cause, father_cause)

fwrite(final_cause_end, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/CAUSE/raw/CAUSE_IR.csv", sep = ",", col.names = TRUE, row.names = FALSE, na = "-")

