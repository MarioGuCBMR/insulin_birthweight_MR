##############
#INTRODUCTION#
##############

#This code takes as input:

#The heterogeneity results for IR and reports them in a supercool table.

###########
#LIBRARIES#
###########

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

loading_rds <- function(path_, file_){
  
  final_path <- paste(path_, file_, sep = "")
  
  final_df <- try(readRDS(final_path), silent = TRUE)
  
  return(final_df)
  
}

round_mgu_ci <- function(digit){
  
  if(is.na(digit)){
    
    return(digit)
    
  }
  
  digit_ <- round(digit, digits = 2)
  
  if(digit_ < 0){
    
    check <- abs(digit_)
    
  } else {
    
    check <- digit_
    
  }
  
  if(nchar(as.character(check)) < 4 & digit_ != 0 & digit_ != 1){
    
    #options(digits=9)
    digit_new <- as.character(paste(as.character(digit_), as.character(0), sep = ""))
    
    return(digit_new)
    
  } else {
    
    return(digit_)
    
  }
  
}



round_mgu_pval <- function(digit){
  
  if(is.na(digit)){
    
    return(digit)
    
  }
  
  if(digit < 0.01){
    
    digit_ <- formatC(digit, format = "e", digits = 2)
    
    return(digit_)
    
  } else {
    
    digit_ <- round(digit, digits = 2)
    
    return(digit_)
    
  }
  
}


organizing_het_data <- function(het_data, rucker_data, het_post = list(0), rucker_post = list(0)){
  
  if(het_post[[1]][1] == 0 & unlist(rucker_post[[1]][1])[1] ==  0){ #this won't ever be 0 in our data, so it works just fine.
    
    cochran_original <- paste(as.character(round_mgu_pval(het_data$Q)), " (P = ", as.character(round_mgu_pval(het_data$pval.Q)), ")", sep = "")
    cochran_post <- "NA"
    isq_original <- paste(as.character(round_mgu_pval(het_data$I2)), " (", as.character(round_mgu_pval(het_data$lower.I2)), ",",  as.character(round_mgu_pval(het_data$upper.I2)), ")", sep = "")
    isq_post <- "NA"
    egger_original <- paste(as.character(round_mgu_pval(rucker_data[[1]]$intercept$Estimate[1])),  " (P = ", as.character(round_mgu_pval(rucker_data[[1]]$intercept$P[1])), ")", sep = "")
    egger_post <- "NA"
    rucker_original_ <- paste(as.character(round_mgu_pval(rucker_data[[1]]$Q$Q[3])),  " (P = ", as.character(round_mgu_pval(rucker_data[[1]]$Q$P[3])), ")", sep = "")
    rucker_post_ <- "NA"
    
  } 
  
  if(het_post[[1]][1] != 0 & unlist(rucker_post[[1]][1])[1] ==  0){ #this won't ever be 0 in our data, so it works just fine.
    
    cochran_original <- paste(as.character(round_mgu_pval(het_data$Q)), " (P = ", as.character(round_mgu_pval(het_data$pval.Q)), ")", sep = "")
    cochran_post <- paste(as.character(round_mgu_pval(het_post$Q)), " (P = ", as.character(round_mgu_pval(het_post$pval.Q)), ")", sep = "")
    isq_original <- paste(as.character(round_mgu_ci(het_data$I2)), " (", as.character(round_mgu_ci(het_data$lower.I2)), ",",  as.character(round_mgu_ci(het_data$upper.I2)), ")", sep = "")
    isq_post <-  paste(as.character(round_mgu_ci(het_post$I2)), " (", as.character(round_mgu_ci(het_post$lower.I2)), ",",  as.character(round_mgu_ci(het_post$upper.I2)), ")", sep = "")
    egger_original <- paste(as.character(round_mgu_pval(rucker_data[[1]]$intercept$Estimate[1])),  " (P = ", as.character(round_mgu_pval(rucker_data[[1]]$intercept$P[1])), ")", sep = "")
    egger_post <- "NA"
    rucker_original_ <- paste(as.character(round_mgu_pval(rucker_data[[1]]$Q$Q[3])),  " (P = ", as.character(round_mgu_pval(rucker_data[[1]]$Q$P[3])), ")", sep = "")
    rucker_post_ <- "NA"
    
  } 
  
  
  if(het_post[[1]][1] != 0 & unlist(rucker_post[[1]][1])[1] !=  0){
    
    cochran_original <- paste(as.character(round_mgu_pval(het_data$Q)), " (P = ", as.character(round_mgu_pval(het_data$pval.Q)), ")", sep = "")
    cochran_post <- paste(as.character(round_mgu_pval(het_post$Q)), " (P = ", as.character(round_mgu_pval(het_post$pval.Q)), ")", sep = "")
    isq_original <- paste(as.character(round_mgu_ci(het_data$I2)), " (", as.character(round_mgu_ci(het_data$lower.I2)), ",",  as.character(round_mgu_ci(het_data$upper.I2)), ")", sep = "")
    isq_post <-  paste(as.character(round_mgu_ci(het_post$I2)), " (", as.character(round_mgu_ci(het_post$lower.I2)), ",",  as.character(round_mgu_ci(het_post$upper.I2)), ")", sep = "")
    egger_original <- paste(as.character(round_mgu_pval(rucker_data[[1]]$intercept$Estimate[1])),  " (P = ", as.character(round_mgu_pval(rucker_data[[1]]$intercept$P[1])), ")", sep = "")
    egger_post <- paste(as.character(round_mgu_pval(rucker_post[[1]]$intercept$Estimate[1])),  " (P = ", as.character(round_mgu_pval(rucker_post[[1]]$intercept$P[1])), ")", sep = "")
    rucker_original_ <- paste(as.character(round_mgu_pval(rucker_data[[1]]$Q$Q[3])),  " (P = ", as.character(round_mgu_pval(rucker_data[[1]]$Q$P[3])), ")", sep = "")
    rucker_post_ <- paste(as.character(round_mgu_pval(rucker_post[[1]]$Q$Q[3])),  " (P = ", as.character(round_mgu_pval(rucker_post[[1]]$Q$P[3])), ")", sep = "")
    
  }
  
  #Finally we return a six-item vector:
  
  cochran_final <- paste(cochran_original, " / ", cochran_post, sep = "") 
  isq_final <- paste(isq_original, " / ", isq_post, sep = "") 
  egger_final <- paste(egger_original, " / ", egger_post, sep = "") 
  rucker_final <- paste(rucker_original_, " / ", rucker_post_, sep = "") 
  
  final_vect <- c(cochran_final, isq_final, egger_final, rucker_final)
  
  return(final_vect)
  
}

heterogeneity_data_parser <- function(exposure, outcome){
  
  original_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/"
  
  exposure_path <- paste(original_path, exposure, "/", outcome, "/", sep = "")
  
  files_in_path <- list.files(exposure_path)
  
  het_input_file_original <- files_in_path[which(str_detect(files_in_path, "Q_Isq_original") == TRUE)]
  het_input_file_post <- files_in_path[which(str_detect(files_in_path, "Q_Isq_post") == TRUE)]
  rucker_input_file_original <- files_in_path[which(str_detect(files_in_path, "rucker_test_original") == TRUE)]
  rucker_input_file_post <- files_in_path[which(str_detect(files_in_path, "rucker_test_post") == TRUE)]
  
  #There will always be "original", but we do not always have post.
  
  if(is_empty(het_input_file_post) == FALSE){
    
    if(is_empty(rucker_input_file_post) == FALSE){
    
    het_file_post <- loading_rds(exposure_path, het_input_file_post)
    rucker_file_post <- loading_rds(exposure_path, rucker_input_file_post)
    het_file_original <- loading_rds(exposure_path, het_input_file_original)
    rucker_file_original <- loading_rds(exposure_path, rucker_input_file_original)
    
    final_vect <- organizing_het_data(het_file_original, rucker_file_original, het_file_post, rucker_file_post)
    
    } else { #this is a case in which Rucker could not be produced because we had less than SNPs
    
      het_file_post <- loading_rds(exposure_path, het_input_file_post)
      het_file_original <- loading_rds(exposure_path, het_input_file_original)
      rucker_file_original <- loading_rds(exposure_path, rucker_input_file_original)
      
      final_vect <- organizing_het_data(het_file_original, rucker_file_original, het_file_post)
      
  }
    
  } else {
    
    het_file_original <- loading_rds(exposure_path, het_input_file_original)
    rucker_file_original <- loading_rds(exposure_path, rucker_input_file_original)
    
    final_vect <- organizing_het_data(het_file_original, rucker_file_original)
    
    
  }
  
  return(final_vect)
  
}

######
#MAIN#
######

#Let's get the data for all our fellas:
#We are going to do it for each outcome:

outcome_ <- c("FBW", "FBWadjMBW", "MBW", "MBWadjFBW", "Father")
exposure_ <- c("CIR", "CIRadjISI", "DI")

#Let's get the data for each outcome.
#To do so we are going to do the 53_IR_Lotta_SNPs first:

IS_FBW <- heterogeneity_data_parser(exposure_[1], outcome_[1])
IS_FBWadjMBW <- heterogeneity_data_parser(exposure_[1], outcome_[2])

IS_MBW <- heterogeneity_data_parser(exposure_[1], outcome_[3])
IS_MBWadjFBW <- heterogeneity_data_parser(exposure_[1], outcome_[4])
IS_Father <- heterogeneity_data_parser(exposure_[1], outcome_[5])

#We change it to dataframes...

IS_FBW_df <- as.data.frame(t(IS_FBW))
IS_FBWadjMBW_df <- as.data.frame(t(IS_FBWadjMBW))
IS_MBW_df <- as.data.frame(t(IS_MBW))
IS_MBWadjFBW_df <- as.data.frame(t(IS_MBWadjFBW))
IS_Father_df <- as.data.frame(t(IS_Father))

#Now we generate a template to have the table perfectly done ;)

template_df <- IS_FBW_df
template_df <- rbind(template_df, template_df)
template_df$V1 <- c("Cochran's Q (P-value)", "(before outlier removal / after outlier removal)")
template_df$V2 <- c("Isq (CI)",  "(before outlier removal / after outlier removal)")
template_df$V3 <- c("Egger's intercept test (P-value)", "(before outlier removal / after outlier removal)")
template_df$V4 <- c("Rucker's test (P-value)", "(before outlier removal / after outlier removal)")

#And now we get these fellas:

IS_FBW_df <- rbind(template_df, IS_FBW_df)
IS_FBWadjMBW_df <- rbind(template_df, IS_FBWadjMBW_df)
IS_MBW_df <- rbind(template_df, IS_MBW_df)
IS_MBWadjFBW_df <- rbind(template_df, IS_MBWadjFBW_df)
IS_Father_df <- rbind(template_df, IS_Father_df)

#Cool! Now we have to do the loop.

for(i in seq(2, 3)){
  
  print(i)
  
  print(1)
  tmp_FBW <- heterogeneity_data_parser(exposure_[i], outcome_[1])
  print(2)
  tmp_FBWadjMBW <- heterogeneity_data_parser(exposure_[i], outcome_[2])
  print(3)
  tmp_MBW <- heterogeneity_data_parser(exposure_[i], outcome_[3])
  print(4)
  tmp_MBWadjFBW <- heterogeneity_data_parser(exposure_[i], outcome_[4])
  print(5)
  tmp_Father <- heterogeneity_data_parser(exposure_[i], outcome_[5])
  
  #We add these fellas...
  
  tmp_FBW_df <- as.data.frame(t(tmp_FBW))
  tmp_FBWadjMBW_df <- as.data.frame(t(tmp_FBWadjMBW))
  tmp_MBW_df <- as.data.frame(t(tmp_MBW))
  tmp_MBWadjFBW_df <- as.data.frame(t(tmp_MBWadjFBW))
  tmp_Father_df <- as.data.frame(t(tmp_Father))
  
  #And now we merge...
  
  IS_FBW_df <- rbind(IS_FBW_df, tmp_FBW_df)
  IS_FBWadjMBW_df <- rbind(IS_FBWadjMBW_df, tmp_FBWadjMBW_df)
  IS_MBW_df <- rbind(IS_MBW_df, tmp_MBW_df)
  IS_MBWadjFBW_df <- rbind(IS_MBWadjFBW_df, tmp_MBWadjFBW_df)
  IS_Father_df <- rbind(IS_Father_df, tmp_Father_df)
  
}

#Let's change the names...

colnames(IS_FBW_df) <- rep("FBW", 4)
colnames(IS_FBWadjMBW_df) <- rep("FBWadjMBW", 4)
colnames(IS_MBW_df) <- rep("MBW", 4)
colnames(IS_MBWadjFBW_df) <- rep("MBWadjFBW", 4)
colnames(IS_Father_df) <- rep("Father", 4)

#And finally we merge the data:

IS_final_df <- cbind(IS_FBW_df, IS_FBWadjMBW_df, IS_MBW_df, IS_MBWadjFBW_df, IS_Father_df)

#We got it :)

fwrite(IS_final_df, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/GRS_&_2SMR/5_heterogeneity_tests/raw/het_2SMR_IS.csv", sep = ",", col.names = TRUE, row.names = FALSE, na = "-")
