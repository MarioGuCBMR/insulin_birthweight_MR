##############
#INTRODUCTION#
##############

#This code is to generate in a comprehensive manner the results for IS GC

###########
#Libraries#
###########

library(data.table)
library(tidyverse)
library(readtext)

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

gc_parser <- function(text_){
  #This function searches in a text the rsults for the Genetic Covariance and returns the beta and the confidence intervals.
  
  #FISst we separete by lines:
  
  tmp <- strsplit(text_, "\n")[[1]]
  
  #We take the GC and the PVAL:
  
  index_gc <- which(str_detect(tmp, "Genetic Correlation:") == TRUE)
  index_p <- which(str_detect(tmp, "P:") == TRUE)
  
  gc <- tmp[index_gc]
  pval <- tmp[index_p]
  
  #We clean the effect size and SE so that we can make the CIs
  
  effect_size <- as.numeric(as.character(strsplit(gc, " ")[[1]][3]))
  se_ <- as.character(strsplit(gc, " ")[[1]][4])
  se_ <- as.character(strsplit(se_, "[(]")[[1]][2])
  se_ <- as.numeric(as.character(strsplit(se_, "[)]")[[1]][1]))
  
  pval_ <- as.character(strsplit(pval, " ")[[1]][2])
  
  lower_ci_ <- lower_ci(effect_size, se_)
  upper_ci_ <- upper_ci(effect_size, se_)
  
  final_vect <- c(effect_size, lower_ci_, upper_ci_, pval_)
  
  return(final_vect)
  
}

###############
#Loading files#
###############

#Now for IS:

CIR_fbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/LDSC/CIR_FBW_LDSC.txt.log")
CIR_mbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/LDSC/CIR_MBW_LDSC.txt.log")
CIR_fbwadjmbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/LDSC/CIR_FBWadj_LDSC.txt.log")
CIR_mbwadjfbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/LDSC/CIR_MBWadj_LDSC.txt.log")
CIR_father <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/LDSC/CIR_Father_LDSC.txt.log")

CIRadjISI_fbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/LDSC/CIRadjISI_FBW_LDSC.txt.log")
CIRadjISI_mbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/LDSC/CIRadjISI_MBW_LDSC.txt.log")
CIRadjISI_fbwadjmbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/LDSC/CIRadjISI_FBWadj_LDSC.txt.log")
CIRadjISI_mbwadjfbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/LDSC/CIRadjISI_MBWadj_LDSC.txt.log")
CIRadjISI_father <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/LDSC/CIRadjISI_Father_LDSC.txt.log")

DI_fbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/LDSC/DI_FBW_LDSC.txt.log")
DI_mbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/LDSC/DI_MBW_LDSC.txt.log")
DI_fbwadjmbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/LDSC/DI_FBWadj_LDSC.txt.log")
DI_mbwadjfbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/LDSC/DI_MBWadj_LDSC.txt.log")
DI_father <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/LDSC/DI_Father_LDSC.txt.log")

#####################################
#Let's parse the genetic correlation#
#####################################

CIR_fbw_gc <- gc_parser(CIR_fbw$text)
CIR_mbw_gc <- gc_parser(CIR_mbw$text)
CIR_fbwadjmbw_gc <- gc_parser(CIR_fbwadjmbw$text)
CIR_mbwadjfbw_gc <- gc_parser(CIR_mbwadjfbw$text)
CIR_father_gc <- gc_parser(CIR_father$text)

CIRadjISI_fbw_gc <- gc_parser(CIRadjISI_fbw$text)
CIRadjISI_mbw_gc <- gc_parser(CIRadjISI_mbw$text)
CIRadjISI_fbwadjmbw_gc <- gc_parser(CIRadjISI_fbwadjmbw$text)
CIRadjISI_mbwadjfbw_gc <- gc_parser(CIRadjISI_mbwadjfbw$text)
CIRadjISI_father_gc <- gc_parser(CIRadjISI_father$text)

DI_fbw_gc <- gc_parser(DI_fbw$text)
DI_mbw_gc <- gc_parser(DI_mbw$text)
DI_fbwadjmbw_gc <- gc_parser(DI_fbwadjmbw$text)
DI_mbwadjfbw_gc <- gc_parser(DI_mbwadjfbw$text)
DI_father_gc <- gc_parser(DI_father$text)

##################################################
#Making a dataframe for each exposure and outcome#
##################################################

#Let's put all data into a list to be able to loop this:

fbw_list <- list(CIR_fbw_gc, CIRadjISI_fbw_gc, DI_fbw_gc)
fbwadjmbw_list <- list(CIR_fbwadjmbw_gc, CIRadjISI_fbwadjmbw_gc, DI_fbwadjmbw_gc)
mbw_list <- list(CIR_mbw_gc, CIRadjISI_mbw_gc, DI_mbw_gc)
mbwadjfbw_list <- list(CIR_mbwadjfbw_gc, CIRadjISI_mbwadjfbw_gc, DI_mbwadjfbw_gc)
father_list <- list(CIR_father_gc, CIRadjISI_father_gc, DI_father_gc)

######################################################################
#Looping to make a dataframe of six rows and two columns: beta and CI#
######################################################################

#We are going to do this with a function that curates our data.
#To speed up the process we are going to do the following

fbw_results_list <- as.data.frame(t(c("Genetic Correlation", "95% LCI", "95% UCI", "P-value")))
fbwadjmbw_results_list <- as.data.frame(t(c("Genetic Correlation", "95% LCI", "95% UCI", "P-value")))
mbw_results_list <- as.data.frame(t(c("Genetic Correlation", "95% LCI", "95% UCI", "P-value")))
mbwadjfbw_results_list <- as.data.frame(t(c("Genetic Correlation", "95% LCI", "95% UCI", "P-value")))
father_results_list <- as.data.frame(t(c("Genetic Correlation", "95% LCI", "95% UCI", "P-value")))

colnames(fbw_results_list) <- c("FBW", "FBW", "FBW", "FBW")
colnames(fbwadjmbw_results_list) <- c("FBWadjMBW", "FBWadjMBW", "FBWadjMBW", "FBWadjMBW")
colnames(mbw_results_list) <- c("MBW", "MBW", "MBW", "MBW")
colnames(mbwadjfbw_results_list) <- c("MBWadjFBW", "MBWadjFBW", "MBWadjFBW", "MBWadjFBW")
colnames(father_results_list) <- c("Father", "Father", "Father", "Father")

for(i in seq(1,3)){
  
  tmp_fbw <- as.data.frame(t(as.data.frame(unlist(fbw_list[i]))))
  colnames(tmp_fbw) <- c("FBW", "FBW", "FBW", "FBW")
  
  tmp_fbwadjmbw <- as.data.frame(t(as.data.frame(unlist(fbwadjmbw_list[i]))))
  colnames(tmp_fbwadjmbw) <- c("FBWadjMBW", "FBWadjMBW", "FBWadjMBW", "FBWadjMBW")
  
  tmp_mbw <- as.data.frame(t(as.data.frame(unlist(mbw_list[i]))))
  colnames(tmp_mbw) <- c("MBW", "MBW", "MBW", "MBW")
  
  tmp_mbwadjfbw <- as.data.frame(t(as.data.frame(unlist(mbwadjfbw_list[i]))))
  colnames(tmp_mbwadjfbw) <- c("MBWadjFBW", "MBWadjFBW", "MBWadjFBW", "MBWadjFBW")
  
  tmp_father <- as.data.frame(t(as.data.frame(unlist(father_list[i]))))
  colnames(tmp_father) <- c("Father", "Father", "Father", "Father")
  
  fbw_results_list <- rbind(fbw_results_list, tmp_fbw)
  fbwadjmbw_results_list <- rbind(fbwadjmbw_results_list, tmp_fbwadjmbw)
  mbw_results_list <- rbind(mbw_results_list, tmp_mbw)
  mbwadjfbw_results_list <- rbind(mbwadjfbw_results_list, tmp_mbwadjfbw)
  father_results_list <- rbind(father_results_list, tmp_father)
  
}


#################################
#Let's clean all of these fellas#
#################################

rownames(fbw_results_list) <- c("value", "CIR", "CIRadjISI", "DI")
rownames(mbw_results_list) <- c("value", "CIR", "CIRadjISI", "DI")
rownames(fbwadjmbw_results_list) <- c("value", "CIR", "CIRadjISI", "DI")
rownames(mbwadjfbw_results_list) <- c("value", "CIR", "CIRadjISI", "DI")
rownames(father_results_list) <- c("value", "CIR", "CIRadjISI", "DI")

####################################################################
#PERFECT now we just need to merge the data and we will have it all#
####################################################################

final_gc_IS_results <- cbind(fbw_results_list, fbwadjmbw_results_list, mbw_results_list, mbwadjfbw_results_list, father_results_list)

View(final_gc_IS_results)

fwrite(final_gc_IS_results, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/GC/results_IS_GC.csv", sep = ",", row.names = TRUE, col.names = TRUE)
