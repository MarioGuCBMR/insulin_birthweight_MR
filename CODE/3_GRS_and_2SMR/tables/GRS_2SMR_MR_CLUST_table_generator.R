##############
#INTRODUCTION#
##############

#This code takes as input:

#an input file for the GRS/2SMR analysis.
#the input file for GRS.
#the input file for 2SMR before outlier extraction
#the input file for 2SMR after outlier extraction
#MR-CLUST results.

#The objective is to have a table that reads the data and makes it as clear as possible
#what happened with that variant.

#The output will have the same format as the input.
#In the example of DI: 8 rows + all the info column.
#Then we will add:

#FBW: GRS+MR
#FBWadjMBW: GRS+MR
#MBW: GRS+MR
#MBWadjFBW: GRS+MR
#Father: GRS+MR

###########
#LIBRARIES#
###########

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

fiadjbmi_input_cleaner <- function(FIadjBMI_input){
  #This function takes the FIadjBMI data and cleans it to have the same nomenclature and format
  #as the rest.
  #It also works for FI!! 
  
  #STEP 1: aligning the data to the FIadjBMI-increasing allele.
  
  new_a1 <- ifelse(as.numeric(FIadjBMI_input$Beta) < 0, FIadjBMI_input$Other_Allele, FIadjBMI_input$Effect_Allele)
  new_a2 <- ifelse(as.numeric(FIadjBMI_input$Beta) < 0, FIadjBMI_input$Effect_Allele, FIadjBMI_input$Other_Allele)
  new_beta <- ifelse(as.numeric(FIadjBMI_input$Beta) < 0, (-1)*(FIadjBMI_input$Beta), FIadjBMI_input$Beta)
  new_eaf <- ifelse(as.numeric(FIadjBMI_input$Beta) < 0, 1-FIadjBMI_input$EAF, FIadjBMI_input$EAF)
  
  FIadjBMI_input_pos <- FIadjBMI_input
  
  FIadjBMI_input_pos$final_effect_allele <- new_a1
  FIadjBMI_input_pos$final_other_allele <- new_a2
  FIadjBMI_input_pos$final_beta <- new_beta
  FIadjBMI_input_pos$final_eaf <- new_eaf

  #STEP 2: select the needed columns:
  
  FIadjBMI_end <- FIadjBMI_input_pos %>%
    select(rsID, Chromosome, Position_build_37, final_effect_allele, final_other_allele, final_eaf, final_beta, SE, P, N)
  
  #STEP 3: change the names:
  
  colnames(FIadjBMI_end) <- c("snp", "CHR", "POS", "Final_A1", "Final_A2", "Final_EAF", "Final_BETA", "stderr", "pvalue",    "N")
  
  #STEP 4: return the data:
  
  return(FIadjBMI_end)
  
}


ir_input_cleaner <- function(IR_input){
  #This function takes the FIadjBMI data and cleans it to have the same nomenclature and format
  #as the rest.
  #It also works for FI!! 
  
  #STEP 1: select the needed columns:
  
  IR_input_end <- IR_input %>%
    select(rsID, Chromosome, Position_build_37, Effect_Allele, Other_Allele, EAF, Beta, SE, P, N, origin)
  
  #STEP 3: change the names:
  
  colnames(IR_input_end) <- c("snp", "CHR", "POS", "Final_A1", "Final_A2", "Final_EAF", "Final_BETA", "stderr", "pvalue",    "N", "origin")
  
  #STEP 4: return the data:
  
  return(IR_input_end)
  
}



grs_data_searcher <- function(exposure, outcome){
  
  original_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/3_GRS_output/"
  
  exposure_path <- paste(original_path, exposure, "/", outcome, "/", sep = "")
  
  files_in_path <- list.files(exposure_path)
  
  grs_input_file <- files_in_path[which(str_detect(files_in_path, "merged") == TRUE)]
  
  final_path <- paste(exposure_path, grs_input_file, sep = "")
  
  final_df <- readRDS(final_path)
  
  return(final_df)
  
}

mrclust_data_searcher <- function(exposure, outcome){
  
  original_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/OUTPUT_DATA/4_2SMR_output/"
  
  exposure_path <- paste(original_path, exposure, "/", outcome, "/", sep = "")
  
  files_in_path <- list.files(exposure_path)
  
  mrclust_input_file <- files_in_path[which(str_detect(files_in_path, "CLUST") == TRUE)]
  
  final_path <- paste(exposure_path, mrclust_input_file, sep = "")
  
  final_df <- readRDS(final_path)
  
  return(final_df)
  
}

adding_cluster_info<- function(file, cluster_data){
  
  final_cluster_vect <- c()
  
  for(snp in file$snp){
    
    if(snp%in%cluster_data$observation){
      
      cluster_data_tmp <- cluster_data[which(cluster_data$observation == snp),]
      
      cluster_data_retrieved <- paste(cluster_data_tmp$cluster_class, "/", cluster_data_tmp$probability, "/", cluster_data_tmp$cluster_mean, sep = "")
      
      final_cluster_vect <- c(final_cluster_vect, cluster_data_retrieved)
      
    } else {
      
      final_cluster_vect <- c(final_cluster_vect, "-")
      
    }
    
  }
  
  return(final_cluster_vect)
  
}
  
updating_grs_2smr <- function(file, exposure){
  #This function cleans the data and gets it ready for the paper.
  #It changes the format so that later in the paper it is easier to parse.
  #Additionally, it adds info on which analysis it has been used.
  #Finally, it adds information on which cluster has been reported.
  
  ###################################################################################
  #STEP 0: adding internal variables to run the function: to be commented afterwards#
  ###################################################################################
  
  #file = FIadjBMI_input 
  #exposure = "FIadjBMI"
  
  ###############################################
  #STEP 1: LOAD INFO FROM GRS/2SMR/MR-CLUST INFO#
  ###############################################
  
  #1.1: first we are going to obtain the data for the GRS:
  
  fbw_grs <- grs_data_searcher(exposure, "FBW")
  fbwadjmbw_grs <- grs_data_searcher(exposure, "FBWadjMBW")
  mbw_grs <- grs_data_searcher(exposure, "MBW")
  mbwadjfbw_grs <- grs_data_searcher(exposure, "MBWadjFBW")
  father_grs <- grs_data_searcher(exposure, "Father")
  
  #1.2: Now we are going to obtain the data for the MR-CLUST:
  
  fbw_cluster <- mrclust_data_searcher(exposure, "FBW")
  fbwadjmbw_cluster <- mrclust_data_searcher(exposure, "FBWadjMBW")
  mbw_cluster <- mrclust_data_searcher(exposure, "MBW")
  mbwadjfbw_cluster <- mrclust_data_searcher(exposure, "MBWadjFBW")
  father_cluster <- mrclust_data_searcher(exposure, "Father")
  
  ################################################
  #STEP 2: put info for GRS and 2SMR and MR-CLUST#
  ################################################
  
  #We need to do it for each outcome.
  #We could make another function, but I think it is OK to do it here.
  
  #FBW:
  
  file$GRS_2SMR_FBW <- ifelse(file$snp%in%fbw_grs$SNP, "GRS", "-")
  file$GRS_2SMR_FBW <- ifelse(file$snp%in%fbw_cluster$observation, "GRS,2SMR,MR-CLUST", "GRS")
  file$cluster_and_effect_FBW <- adding_cluster_info(file, fbw_cluster)
    
  #FBWadjMBW.
    
  file$GRS_2SMR_FBWadjMBW <- ifelse(file$snp%in%fbwadjmbw_grs$SNP, "GRS", "-")
  file$GRS_2SMR_FBWadjMBW <- ifelse(file$snp%in%fbwadjmbw_cluster$observation, "GRS,2SMR, MR-CLUST", "GRS")
  file$cluster_and_effect_FBWadjMBW <- adding_cluster_info(file, fbwadjmbw_cluster)
  
  #MBW
  
  file$GRS_2SMR_MBW <- ifelse(file$snp%in%mbw_grs$SNP, "GRS", "-")
  file$GRS_2SMR_MBW <- ifelse(file$snp%in%mbw_cluster$observation, "GRS,2SMR,MR-CLUST", "GRS")
  file$cluster_and_effect_MBW <- adding_cluster_info(file, mbw_cluster)
  
  #MBWadjFBW
  
  file$GRS_2SMR_MBWadjFBW <- ifelse(file$snp%in%mbwadjfbw_grs$SNP, "GRS", "-")
  file$GRS_2SMR_MBWadjFBW <- ifelse(file$snp%in%mbwadjfbw_cluster$observation, "GRS,2SMR,MR-CLUST", "GRS")
  file$cluster_and_effect_MBWadjFBW <- adding_cluster_info(file, mbwadjfbw_cluster)
  
  #Father
  
  file$GRS_2SMR_Father <- ifelse(file$snp%in%father_grs$SNP, "GRS", "-")
  file$GRS_2SMR_Father <- ifelse(file$snp%in%father_cluster$observation, "GRS,2SMR, MR-CLUST", "GRS")
  file$cluster_and_effect_Father <- adding_cluster_info(file, father_cluster)

  ###################################################################################
  #STEP 3: CHANGE THE FORMAT OF THE FILE SO THAT WE CAN WORK WITH IT BETTER IN EXCEL#
  ###################################################################################
  
  #Now we can split the data that we have:
  
  index_cols_info <- which(str_detect(colnames(file), "GRS") == FALSE & str_detect(colnames(file), "cluster") == FALSE)
  index_cols_results <- which(str_detect(colnames(file), "GRS") == TRUE | str_detect(colnames(file), "cluster") == TRUE)
  
  file_snp_info <- file[, ..index_cols_info]
  file_results <- file[, ..index_cols_results]
  
  #First let's change our SNPs info:
  
  colnames_ <- colnames(file_snp_info)
  colnames_ <- as.data.frame(t(colnames_))
  colnames(colnames_) <- colnames_
  
  file_formated_1 <- rbind(colnames_, file_snp_info)
  
  colnames(file_formated_1) <- rep("SNP_INFO", length(colnames(file_formated_1))) #this will help us distribute the info better in the excel sheet.
  
  #Now let's change the data:
  
  colnames_ <- rep(c("Analysis", "Cluster/Probability/cluster's mean causal effect"), 5)
  colnames(file_results) <- unlist(colnames_)
  colnames_ <- as.data.frame(t(colnames_))
  colnames(colnames_) <- colnames_
  
  file_formated_2 <- rbind(colnames_, file_results)
  
  colnames(file_formated_2) <- c("FBW", "FBW", "FBWadjMBW", "FBWadjMBW", "MBW", "MBW", "MBWadjFBW", "MBWadjFBW", "Father", "Father")
  
  ###################################
  #Let's mix this data: this is over#
  ###################################
  
  file_formated_end <- cbind(file_formated_1, file_formated_2)
  
  return(file_formated_end)
  
}

######
#MAIN#
######

#################################
#Let's do it for IS traits first#
#################################

#Let's first do DI:

DI_input <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/DI/CLEAN_7_DI_GRS_SNPs_Manuscript.csv")
DI_input$N <- 5130
DI_updated <- updating_grs_2smr(DI_input, "DI")
fwrite(DI_updated, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/GRS_&_2SMR/2_updated_data_4_GRS_&_2MSR_and_MRCLUST_results/UPDATED_CLEAN_7_DI_GRS_SNPs_Manuscript.csv")

#Now let's do CIR:

#We will change the column names for CIR to those with DI because they do not match...

CIR_input <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/CIR/CLEAN_8_CIR_GRS_SNPs_Manuscript.csv")
CIR_input$N <- 5318
colnames(CIR_input) <- colnames(DI_input)
CIR_updated <- updating_grs_2smr(CIR_input, "CIR")
fwrite(CIR_updated, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/GRS_&_2SMR/2_updated_data_4_GRS_&_2MSR_and_MRCLUST_results/UPDATED_CLEAN_8_CIR_GRS_SNPs_Manuscript.csv")

#Finally, let's do CIRadjISI:

CIRadjISI_input <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/CIRadjISI/CLEAN_8_CIRadjISI_GRS_SNPs_Manuscript.csv")
CIRadjISI_input$N <- 4789
CIRadjISI_updated <- updating_grs_2smr(CIRadjISI_input, "CIRadjISI")
fwrite(CIRadjISI_updated, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/GRS_&_2SMR/2_updated_data_4_GRS_&_2MSR_and_MRCLUST_results/UPDATED_CLEAN_8_CIRadjISI_GRS_SNPs_Manuscript.csv")

#################################
#Let's do it for IR traits first#
#################################

#We change to FIadjBMI

FIadjBMI_input <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/FIadjBMI/FIadjBMI_independent_4_GRS.csv")
FIadjBMI_input <- fiadjbmi_input_cleaner(FIadjBMI_input)
FIadjBMI_input <- updating_grs_2smr(FIadjBMI_input, "FIadjBMI")
fwrite(FIadjBMI_input, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/GRS_&_2SMR/2_updated_data_4_GRS_&_2MSR_and_MRCLUST_results/UPDATED_CLEAN_31_FIadjBMI_GRS_SNPs_Manuscript.csv")

#We change to FI:

FI_input <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/FI/CLEAN_19_FI_GRS_SNPs_Manuscript.csv")
FI_input <- fiadjbmi_input_cleaner(FI_input)
FI_updated <- updating_grs_2smr(FI_input, "FI")
fwrite(FI_updated, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/GRS_&_2SMR/2_updated_data_4_GRS_&_2MSR_and_MRCLUST_results/UPDATED_CLEAN_19_FI_GRS_SNPs_Manuscript.csv")

#Finally, let's do IR:

IR_input <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/CURATED_DATA/Data_4_2SMR/53_IR_Lotta/CLEAN_53_LOTTA_GRS_SNPs_Manuscript.csv")
IR_input <- ir_input_cleaner(IR_input)
IR_updated <- updating_grs_2smr(IR_input, "53_IR_Lotta_SNPs")
fwrite(IR_updated, "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/GRS_&_2SMR/2_updated_data_4_GRS_&_2MSR_and_MRCLUST_results/UPDATED_CLEAN_53_IR_GRS_SNPs_Manuscript.csv")
