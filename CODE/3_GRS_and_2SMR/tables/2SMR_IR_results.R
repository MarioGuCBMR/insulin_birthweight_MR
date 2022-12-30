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

#We want the raw data for the tables, we can deal with the numbers ourselves later.

#round_mgu_ci <- function(digit){
#  
#  if(is.na(digit)){
#    
#    return(digit)
#    
#
#  }
#  
#  digit_ <- round(digit, digits = 2)
#  
#  if(digit_ < 0){
#    
#    check <- abs(digit_)
#    
#  } else {
#    
#    check <- digit_
#    
#  }
#  
#  if(nchar(as.character(check)) < 4 & digit_ != 0 & digit_ != 1){
#    
#    #options(digits=9)
#    digit_new <- as.character(paste(as.character(digit_), as.character(0), sep = ""))
#    
#    return(digit_new)
#    
#  } else {
#    
#    return(digit_)
#    
#  }
#  
#}

#round_mgu_pval <- function(digit){
#  
#  if(is.na(digit)){
#    
#    return(digit)
#    
#  }
#  
#  if(digit < 0.01){
#    
#    digit_ <- formatC(digit, format = "e", digits = 2)
#    
#    return(digit_)
#    
#  } else {
#    
#    digit_ <- round(digit, digits = 2)
#    
#    return(digit_)
#    
#  }
#  
#}

mr_file_finder <- function(list_of_files, path_, exposure, outcome){
  
  #First, let's check if we have post files:
  
  check <- length(which(str_detect(list_of_files, "post") == TRUE))
  
  #If we do we are gonna work with post
  
  if(check > 0){
    
    file_pattern <- "post"
    
  } else {
    
    file_pattern <- "original"
    
  }
  
  #Let's take the files that we are interested in
  
  files_we_want <- list_of_files[which(str_detect(list_of_files, file_pattern) == TRUE)]
  
  #And now... let's remove the annoying ones:
  
  files_we_want <- files_we_want[which(str_detect(files_we_want, "tiff") == FALSE)] #out with the figures
  files_we_want <- files_we_want[which(str_detect(files_we_want, "merged") == FALSE)] #remove the data
  files_we_want <- files_we_want[which(str_detect(files_we_want, "rucker") == FALSE)] #remove the het 1/2
  files_we_want <- files_we_want[which(str_detect(files_we_want, "Isq") == FALSE)] #remove the het 2/2
  
  setwd(path_)
  
  final_df <- readRDS(files_we_want)
  
  final_df$lower_ci <- lower_ci(final_df$b, final_df$se)
  final_df$upper_ci <- upper_ci(final_df$b, final_df$se)
  
  #Indexing so that IVW comes first:
  
  if(length(final_df$b) > 1){
    
    index_method <- c(3, 1, 2, 5)
    
    final_df <- final_df[index_method,]
    
  }
  
  #Before we do anything, we are going to order them in IVW, 
  
  #Now let's clean this data a little bit
  
  if(outcome == "FBW"){
  
  final_df_clean <- final_df %>%
    select(exposure, method, nsnp, b, lower_ci, upper_ci, se, pval)
  
  } else {
    
    final_df_clean <- final_df %>%
      select(method, nsnp, b, lower_ci, upper_ci, se, pval)
    
  }
  
  #Finally, we are going to do some tricks to make the data beautiful
  #So that we can save it properly, but only for FI, which is the first one:
  
  if(exposure == "FI"){
  
     final_df_clean_copy <- final_df_clean
     final_df_clean <- final_df_clean[-seq(1,length(final_df_clean$method)),]
     final_df_clean[1,] <- colnames(final_df_clean)
     final_df_clean <- rbind(final_df_clean, final_df_clean_copy)
  
  }
  
  #Let's remove "Simple Mode":
  
  if("Simple mode"%in%final_df_clean$method){
    
    final_df_clean <- final_df_clean[which(final_df_clean$method != "Simple mode"),]
    
  }
  
  return(final_df_clean)
  
}

mr_results_parser <- function(outcome){
  
  original_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/"
  
  #Now let's loop this baby through the outcomes:
  
  exposures <- c("FI", "FIadjBMI", "HOMAIR")
  
  for(exposure in exposures){
    
    if(exposure == "FI"){
  
      exposure_path <- paste(original_path, exposure, "/", outcome, "/", sep = "")
  
      files_in_path <- list.files(exposure_path)
  
      results_end <- mr_file_finder(files_in_path, exposure_path, exposure, outcome)
  
    } else {
      
      exposure_path <- paste(original_path, exposure, "/", outcome, "/", sep = "")
      
      files_in_path <- list.files(exposure_path)
      
      results <- mr_file_finder(files_in_path, exposure_path, exposure, outcome)
      
      results_end <- rbind(results_end, results)
      
    }
  
  }
  
  #And now we change the names of the colums, by the name of the outcome, so that we can clean things
  #properly later :)
  
  colnames(results_end) <- rep(outcome, length(colnames(results_end)))
  
  return(results_end)

}


######
#MAIN#
######

IR_FBW <- mr_results_parser("FBW")
IR_FBWadjMBW <- mr_results_parser("FBWadjMBW")
IR_MBW <- mr_results_parser("MBW")
IR_MBWadjFBW <- mr_results_parser("MBWadjFBW")
IR_Father <- mr_results_parser("Father")

IR_df <- cbind(IR_FBW, IR_FBWadjMBW, IR_MBW, IR_MBWadjFBW, IR_Father)

#We save these fellas:

fwrite(IR_df, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/GRS_&_2SMR/4_2SMR_results/2SMR_IR_results_w_homair.csv", sep = ",", col.names = TRUE, row.names = FALSE, na = "-")

















