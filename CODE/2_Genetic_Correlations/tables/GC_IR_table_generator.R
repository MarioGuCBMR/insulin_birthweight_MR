##############
#INTRODUCTION#
##############

#This code is to generate in a comprehensive manner the results for IR GC

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
  
  #First we separete by lines:
  
  tmp <- strsplit(text_, "\n")[[1]]
  
  #We take the GC and the PVAL:
  
  index_gc <- which(str_detect(tmp, "Genetic Correlation:") == TRUE)
  index_p <- which(str_detect(tmp, "P:") == TRUE)
  
  gc <- tmp[index_gc]
  pval <- tmp[index_p]
  
  #We clean the effect size and SE so that we can make the CIs
  
  effect_size <- as.numeric(as.character(strsplit(gc, " ")[[1]][4]))
  se_ <- as.character(strsplit(gc, " ")[[1]][5])
  se_ <- as.character(strsplit(se_, "[(]")[[1]][2])
  se_ <- as.numeric(as.character(strsplit(se_, "[)]")[[1]][1]))
  
  pval_ <- as.character(strsplit(pval, " ")[[1]][3])
  
  lower_ci_ <- lower_ci(effect_size, se_)
  upper_ci_ <- upper_ci(effect_size, se_)
  
  final_vect <- c(effect_size, lower_ci_, upper_ci_, pval_)
  
  return(final_vect)
  
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
      
      final_value <- formatC(beta_ci_raw[value_index], format = "e", digits = 2)
      
      beta_ci_vect <- c(beta_ci_vect, final_value)
      
    } else {
      
      final_value <- round(as.numeric(beta_ci_raw[value_index]), digits = 2)
      
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

#Now for IR:

FIadjBMI_fbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/fiadjbmi_FBW.txt")
FIadjBMI_mbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/fiadjbmi_MBW.txt")
FIadjBMI_fbwadjmbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/fiadjbmi_FBWadj.txt")
FIadjBMI_mbwadjfbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/fiadjbmi_MBWadj.txt")
FIadjBMI_father <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/fiadjbmi_Father.txt")

FI_fbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/FI_FBW.txt")
FI_mbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/FI_MBW.txt")
FI_fbwadjmbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/FI_FBWadj.txt")
FI_mbwadjfbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/FI_MBWadj.txt")
FI_father <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/FI_Father.txt")

homair_fbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/homair_FBW.txt")
homair_mbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/homair_MBW.txt")
homair_fbwadjmbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/homair_FBWadj.txt")
homair_mbwadjfbw <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/homair_MBWadj.txt")
homair_father <- readtext("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CODE/2_Genetic_Correlations/output/homair_Father.txt")

#####################################
#Let's parse the genetic correlation#
#####################################

fiadjbmi_fbw_gc <- gc_parser(FIadjBMI_fbw$text)
fiadjbmi_mbw_gc <- gc_parser(FIadjBMI_mbw$text)
fiadjbmi_fbwadjmbw_gc <- gc_parser(FIadjBMI_fbwadjmbw$text)
fiadjbmi_mbwadjfbw_gc <- gc_parser(FIadjBMI_mbwadjfbw$text)
fiadjbmi_father_gc <- gc_parser(FIadjBMI_father$text)

fi_fbw_gc <- gc_parser(FI_fbw$text)
fi_mbw_gc <- gc_parser(FI_mbw$text)
fi_fbwadjmbw_gc <- gc_parser(FI_fbwadjmbw$text)
fi_mbwadjfbw_gc <- gc_parser(FI_mbwadjfbw$text)
fi_father_gc <- gc_parser(FI_father$text)

homair_fbw_gc <- gc_parser(homair_fbw$text)
homair_mbw_gc <- gc_parser(homair_mbw$text)
homair_fbwadjmbw_gc <- gc_parser(homair_fbwadjmbw$text)
homair_mbwadjfbw_gc <- gc_parser(homair_mbwadjfbw$text)
homair_father_gc <- gc_parser(homair_father$text)

##################################################
#Making a dataframe for each exposure and outcome#
##################################################

#Let's put all data into a list to be able to loop this:

fbw_list <- list(fiadjbmi_fbw_gc, fi_fbw_gc, homair_fbw_gc)
fbwadjmbw_list <- list(fiadjbmi_fbwadjmbw_gc, fi_fbwadjmbw_gc, homair_fbwadjmbw_gc)
mbw_list <- list(fiadjbmi_mbw_gc, fi_mbw_gc, homair_mbw_gc)
mbwadjfbw_list <- list(fiadjbmi_mbwadjfbw_gc, fi_mbwadjfbw_gc, homair_mbwadjfbw_gc)
father_list <- list(fiadjbmi_father_gc, fi_father_gc, homair_father_gc)

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

rownames(fbw_results_list) <- c("value", "FIadjBMI", "FI", "HOMA-IR")
rownames(mbw_results_list) <- c("value", "FIadjBMI", "FI", "HOMA-IR")
rownames(fbwadjmbw_results_list) <- c("value", "FIadjBMI", "FI", "HOMA-IR")
rownames(mbwadjfbw_results_list) <- c("value", "FIadjBMI", "FI", "HOMA-IR")
rownames(father_results_list) <- c("value", "FIadjBMI", "FI", "HOMA-IR")

####################################################################
#PERFECT now we just need to merge the data and we will have it all#
####################################################################

final_gc_ir_results <- cbind(fbw_results_list, fbwadjmbw_results_list, mbw_results_list, mbwadjfbw_results_list, father_results_list)

View(final_gc_ir_results)

fwrite(final_gc_ir_results, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/GC/results_IR_GC.csv", sep = ",", row.names = TRUE, col.names = TRUE)
