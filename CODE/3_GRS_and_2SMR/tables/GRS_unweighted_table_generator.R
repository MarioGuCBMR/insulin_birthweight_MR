##############
#INTRODUCTION#
##############

#This code is to generate in a comprehensive manner the results for unweighted GRS.

###########
#Libraries#
###########

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

lower_ci <- function(beta_, se_){
  #A function that takes the beta and se from summary statistics
  #and generates the lower CI.
  
  lower_ci <- beta_ - qnorm(0.975)*se_
  
  return(lower_ci)
  
}

upper_ci <- function(beta_, se_){
  #A function that takes the beta and se from summary statistics
  #and generates the upper CI.
  
  upper_ci <- beta_ + qnorm(0.975)*se_
  
  return(upper_ci)
  
}

results_parser <- function(results_list){
  #This function takes the result for GRS and extracts, beta, SE and PVAL.
  #Then, it calculate the confindence intervals.
  #and finally it outputs a vector of two values:
  #1) a string like: "beta (CI_l, CI_u)"
  #2) a pval.
  
  #STEP 1: unlist the object:
  
  results_df <- results_list[[1]]
  
  #STEP 2: define the confidence intervals:
  
  results_df$lower_ci <- lower_ci(results_df$ahat, results_df$aSE)
  results_df$upper_ci <- upper_ci(results_df$ahat, results_df$aSE)
  
  #STEP 3: generate beta (CI_l, CI_u) string and clean p-value:
  
  #First we need to check if any value (beta, lower_ci or upper_ci need some formatting tweak)
  
  beta_ci_raw <- c(results_df$ahat, results_df$lower_ci, results_df$upper_ci)
  
  index_format_e <- which(abs(beta_ci_raw) < 0.01) #position of the one that needs especial formatting.
  
  beta_ci_vect <- c() #we open a vector to store the results.
  
  for(value_index in seq(1,3)){
    
    if(value_index%in%index_format_e){
      
      final_value <- formatC(beta_ci_raw[value_index], format = "e", digits = 4)
      
      beta_ci_vect <- c(beta_ci_vect, final_value)
      
    } else {
      
      final_value <- round(as.numeric(beta_ci_raw[value_index]), digits = 4)
      
      beta_ci_vect <- c(beta_ci_vect, final_value)
      
    }
    
  }
  
  #Finally, let's save the data into beta_ci vect
  
  beta_ci <- paste(beta_ci_vect[1], " (", beta_ci_vect[2], ":", beta_ci_vect[3], ")", sep = "")
  
  #Finally, we do something similar with the p-value, but it is way easier:
  
  if(results_df$pval < 0.01){
    
    pval_ <- formatC(results_df$pval, format = "e", digits = 2)
    
  } else {
    
    pval_ <- formatC(round(results_df$pval, digits = 2), digits = 2)
    
  }
  
  #STEP 4: return a vector:
  
  final_vect <- c(beta_ci, pval_)
  
  return(final_vect)
  
}


###############
#Loading files#
###############

#First for IS traits:

CIR_fbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/CIR/FBW/CIR_FBW_unweighted")
CIR_mbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/CIR/MBW/CIR_MBW_unweighted")
CIR_fbwadjmbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/CIR/FBWadjMBW/CIR_FBWadjMBW_unweighted")
CIR_mbwadjfbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/CIR/MBWadjFBW/CIR_mbwadjfbw_unweighted")
CIR_father <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/CIR/Father/CIR_father_unweighted")

CIRadjISI_fbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/CIRadjISI/FBW/CIRadjISI_FBW_unweighted")
CIRadjISI_mbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/CIRadjISI/MBW/CIRadjISI_MBW_unweighted")
CIRadjISI_fbwadjmbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/CIRadjISI/FBWadjMBW/CIRadjISI_FBWadjMBW_unweighted")
CIRadjISI_mbwadjfbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/CIRadjISI/MBWadjFBW/CIRadjISI_mbwadjfbw_unweighted")
CIRadjISI_father <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/CIRadjISI/Father/CIRadjISI_father_unweighted")

DI_fbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/DI/FBW/DI_FBW_unweighted")
DI_mbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/DI/MBW/DI_MBW_unweighted")
DI_fbwadjmbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/DI/FBWadjMBW/DI_FBWadjMBW_unweighted")
DI_mbwadjfbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/DI/MBWadjFBW/DI_mbwadjfbw_unweighted")
DI_father <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/DI/Father/DI_father_unweighted")

#Now for IR:

FIadjBMI_fbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FIadjBMI/FBW/FIadjBMI_FBW_unweighted")
FIadjBMI_mbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FIadjBMI/MBW/FIadjBMI_MBW_unweighted")
FIadjBMI_fbwadjmbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FIadjBMI/FBWadjMBW/FIadjBMI_FBWadjMBW_unweighted")
FIadjBMI_mbwadjfbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FIadjBMI/MBWadjFBW/FIadjBMI_mbwadjfbw_unweighted")
FIadjBMI_father <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FIadjBMI/Father/FIadjBMI_father_unweighted")

FI_fbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FI/FBW/FI_FBW_unweighted")
FI_mbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FI/MBW/FI_MBW_unweighted")
FI_fbwadjmbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FI/FBWadjMBW/FI_FBWadjMBW_unweighted")
FI_mbwadjfbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FI/MBWadjFBW/FI_mbwadjfbw_unweighted")
FI_father <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/FI/Father/FI_father_unweighted")

homair_fbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/HOMAIR/FBW/HOMAIR_FBW_unweighted")
homair_mbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/HOMAIR/MBW/HOMAIR_MBW_unweighted")
homair_fbwadjmbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/HOMAIR/FBWadjMBW/HOMAIR_FBWadjMBW_unweighted")
homair_mbwadjfbw <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/HOMAIR/MBWadjFBW/HOMAIR_mbwadjfbw_unweighted")
homair_father <- readRDS("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/HOMAIR/Father/HOMAIR_father_unweighted")

##################################################
#Making a dataframe for each exposure and outcome#
##################################################

#Let's put all data into a list to be able to loop this:

fbw_list <- list(FI_fbw, FIadjBMI_fbw, homair_fbw, CIR_fbw, CIRadjISI_fbw, DI_fbw)
fbwadjmbw_list <- list(FI_fbwadjmbw, FIadjBMI_fbwadjmbw, homair_fbwadjmbw, CIR_fbwadjmbw, CIRadjISI_fbwadjmbw, DI_fbwadjmbw)
mbw_list <- list(FI_mbw, FIadjBMI_mbw, homair_mbw, CIR_mbw, CIRadjISI_mbw, DI_mbw)
mbwadjfbw_list <- list(FI_mbwadjfbw, FIadjBMI_mbwadjfbw, homair_mbwadjfbw, CIR_mbwadjfbw, CIRadjISI_mbwadjfbw, DI_mbwadjfbw)
father_list <- list(FI_father, FIadjBMI_father, homair_father, CIR_father, CIRadjISI_father, DI_father)

######################################################################
#Looping to make a dataframe of six rows and two columns: beta and CI#
######################################################################

#We are going to do this with a function that curates our data.
#To speed up the process we are going to do the following

fbw_results_list <- as.data.frame(t(c("GRS effect size (95% CI)", "P-value")))
fbwadjmbw_results_list <- as.data.frame(t(c("GRS effect size (95% CI)", "P-value")))
mbw_results_list <- as.data.frame(t(c("GRS effect size (95% CI)", "P-value")))
mbwadjfbw_results_list <- as.data.frame(t(c("GRS effect size (95% CI)", "P-value")))
father_results_list <- as.data.frame(t(c("GRS effect size (95% CI)", "P-value")))

colnames(fbw_results_list) <- c("FBW", "FBW")
colnames(fbwadjmbw_results_list) <- c("FBWadjMBW", "FBWadjMBW")
colnames(mbw_results_list) <- c("MBW", "MBW")
colnames(mbwadjfbw_results_list) <- c("MBWadjFBW", "MBWadjFBW")
colnames(father_results_list) <- c("Father", "Father")

for(res_index in seq(1,6)){
  
  fbw_results_list_tmp <- results_parser(fbw_list[res_index])
  fbwadjmbw_results_list_tmp <- results_parser(fbwadjmbw_list[res_index])
  mbw_results_list_tmp <- results_parser(mbw_list[res_index])
  mbwadjfbw_results_list_tmp <- results_parser(mbwadjfbw_list[res_index])
  father_results_list_tmp <- results_parser(father_list[res_index])
  
  #And now we update the original df:
  
  fbw_results_list <- rbind(fbw_results_list, fbw_results_list_tmp)
  fbwadjmbw_results_list <- rbind(fbwadjmbw_results_list, fbwadjmbw_results_list_tmp)
  mbw_results_list <- rbind(mbw_results_list, mbw_results_list_tmp)
  mbwadjfbw_results_list <- rbind(mbwadjfbw_results_list, mbwadjfbw_results_list_tmp)
  father_results_list <- rbind(father_results_list, father_results_list_tmp)
  
}

#################################
#Let's clean all of these fellas#
#################################

rownames(fbw_results_list) <- c("value", "FI", "FIadjBMI", "HOMA-IR", "CIR", "CIRadjISI", "DI")
rownames(mbw_results_list) <- c("value", "FI", "FIadjBMI", "HOMA-IR",  "CIR", "CIRadjISI", "DI")
rownames(fbwadjmbw_results_list) <- c("value",  "FI", "FIadjBMI", "HOMA-IR", "CIR", "CIRadjISI", "DI")
rownames(mbwadjfbw_results_list) <- c("value",  "FI", "FIadjBMI", "HOMA-IR", "CIR", "CIRadjISI", "DI")
rownames(father_results_list) <- c("value",  "FI", "FIadjBMI", "HOMA-IR", "CIR", "CIRadjISI", "DI")

####################################################################
#PERFECT now we just need to merge the data and we will have it all#
####################################################################

final_grs_results <- cbind(fbw_results_list, fbwadjmbw_results_list, mbw_results_list, mbwadjfbw_results_list, father_results_list)

View(final_grs_results)

fwrite(final_grs_results, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/GRS_&_2SMR/3_results_4_GRS/results_GRS_unweighted_w_homair.csv", row.names = TRUE)
