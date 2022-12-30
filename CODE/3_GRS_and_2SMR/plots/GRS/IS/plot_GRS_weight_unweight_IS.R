##############
#INTRODUCTION#
##############

#This code will take all the GRS results and make forest plots automatically.
#This is going to take the weighted GRS only.
#We will do the unweighted later.

###########
#Libraries#
###########

library(tidyverse)
library(ggplot2)
library(data.table)
library(gridExtra)

###################
#Loading functions#
###################

beta_parser <- function(effect_size){
  
  tmp <- as.character(strsplit(effect_size, " ")[[1]][1])
  
  return(tmp)
  
}

lower_ci_parser <- function(effect_size){
  
  tmp <- as.character(strsplit(effect_size, ":")[[1]][1])
  tmp <- as.character(strsplit(tmp, "[(]")[[1]][2])
  
  return(tmp)
  
}

upper_ci_parser <- function(effect_size){
  
  tmp <- as.character(strsplit(effect_size, ":")[[1]][2])
  tmp <- as.character(strsplit(tmp, "[)]")[[1]][1])
  
  return(tmp)
  
}

data_parsers <- function(results){
  #This function is going to return a whole dataframe with the following vectors:
  #exposure
  #outcome
  #beta
  #lower_ci
  #upper_ci
  #id (exposure-outcome)
  #pvals
  
  #First we make empty vectors to save all:
  
  outcomes <- c()
  exposures <- c()
  betas <- c()
  lower_cis <- c()
  upper_cis <- c()
  ids <- c()
  pvals <- c()
  effects <- c()
  
  #First we take the exposures!!
  
  exposure <- results$V1
  exposure <- exposure[-1]
  exposures <- rep(exposure, 5)
  
  #and we get the final results without the exposures:
  
  results <- results[,-1]
  
  for(i in seq(1,10, by = 2)){
    
    j=i+1
    
    tmp=results[,seq(i,j)]
    
    #Now we are going to obtain the outcome, the beta, the lower_ci, upper_ci and pval!
    
    #First we get the outcomes, which are in the colnames
    
    outcome <- colnames(tmp)[1]
    outcome <- rep(outcome, 3)
    outcomes <- c(outcomes, outcome)
    
    #Now we remove these colnames! And get the data ready to extract the numbers:
    
    col_name_tmp <- tmp[1,]
    colnames(tmp) <- col_name_tmp
    tmp <- tmp[-1,] #removing the colnames now that we have them
    
    beta <- as.numeric(as.character(unlist(sapply(tmp$`GRS effect size (95% CI)`, beta_parser))))
    lower_ci <- as.numeric(as.character(unlist(sapply(tmp$`GRS effect size (95% CI)`, lower_ci_parser))))
    upper_ci <- as.numeric(as.character(unlist(sapply(tmp$`GRS effect size (95% CI)`, upper_ci_parser))))
    
    betas <- c(betas, beta)
    lower_cis <- c(lower_cis, lower_ci)
    upper_cis <- c(upper_cis, upper_ci)
    pvals <- c(pvals, tmp$`P-value`)
    effects <- c(effects, tmp$`GRS effect size (95% CI)`)
    
    
    #And here we end the loop!
    
  }
  
  #And now we loop to get the ids
  
  for(index in seq(1, length(betas))){
    
    id <- paste(exposures[index], "-", outcomes[index], sep = "")
    
    ids <- c(ids, id)
    
  }
  
  ids <- as.data.frame(ids)
  exposures <- as.data.frame(exposures)
  outcomes <- as.data.frame(outcomes)
  betas <- as.data.frame(betas)
  lower_cis <- as.data.frame(lower_cis)
  upper_cis <- as.data.frame(upper_cis)
  pvals <- as.data.frame(pvals)
  
  #Before that we are going to clean the effects:
  effects <- gsub(" - ", ",", effects)
  
  effects <- as.data.frame(effects)
  
  
  final_df <- cbind(ids, exposures, outcomes, betas, lower_cis, upper_cis, pvals, effects)
  
  return(final_df)
  
}

round_mgu <- function(digit){
  
  digit_ <- round(digit, digits = 2)
  
  if(digit_ < 0){
    
    check <- abs(digit_)
    
  } else {
    
    check <- digit_
    
  }
  
  if(nchar(as.character(check)) < 4){
    
    #options(digits=9)
    digit_new <- as.character(paste(as.character(digit_), as.character(0), sep = ""))
    
    return(digit_new)
    
  } else {
    
    return(digit_)
    
  }
  
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

##############
#Loading data#
##############

grs_results_weighted <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/GRS_&_2SMR/3_results_4_GRS/results_GRS_weighted.csv", header = TRUE)
grs_results_unweighted <- fread("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/supplementary_tables_per_method/GRS_&_2SMR/3_results_4_GRS/results_GRS_unweighted.csv", header = TRUE)

grs_results_w <- as.data.frame(grs_results_weighted)
grs_results_unw <- as.data.frame(grs_results_unweighted)

grs_results_w_is <- grs_results_w[c(1,5,6,7),]
grs_results_unw_is <- grs_results_unw[c(1,5,6,7),]

grs_results_w_is <- data_parsers(grs_results_w_is)
grs_results_unw_is <- data_parsers(grs_results_unw_is)

#########################################################################################################
#Due to some issues with the data we are going to clean the effects that we are going to put in the plot#
#########################################################################################################

#First weighted betas:

weighted_betas <- as.character(unlist(sapply(grs_results_w_is$betas, round_mgu_ci)))
weighted_lower_ci <- as.character(unlist(sapply(grs_results_w_is$lower_cis, round_mgu_ci)))
weighted_upper_ci <- as.character(unlist(sapply(grs_results_w_is$upper_cis, round_mgu_ci)))

grs_results_w_is$effects <- paste(weighted_betas, " (",  weighted_lower_ci, ";", weighted_upper_ci, ")", sep = "")

#Then unweighted betas:

unweighted_betas <- as.character(unlist(sapply(grs_results_unw_is$betas, round_mgu_pval)))
unweighted_lower_ci <- as.character(unlist(sapply(grs_results_unw_is$lower_cis, round_mgu_pval)))
unweighted_upper_ci <- as.character(unlist(sapply(grs_results_unw_is$upper_cis, round_mgu_pval)))

grs_results_unw_is$effects <- paste(unweighted_betas, " (",  unweighted_lower_ci, ";", unweighted_upper_ci, ")", sep = "")

#################################
#Make plots for weighted section#
#################################

grs_results_og <- grs_results_w_is

grs_results_og$alphabetical_id <- c("a1", "a2", "a3", "a4", "a5",
                                         "b1", "b2", "b3", "b4", "b5",
                                         "c1", "c2", "c3", "c4", "c5")

plot1 <- ggplot(grs_results_og, aes(y = fct_rev(alphabetical_id), x = as.numeric(betas), shape = exposures)) +
  geom_point(size = 2.7, position = "identity") +  
  geom_errorbarh(aes(xmin = as.numeric(lower_cis), xmax = as.numeric(upper_cis)), height = 0.25) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("GRS effect size (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_blank(),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.ticks.y = element_blank(),
        legend.position = "none")



grs_results_og$final_ids <- c("", "FBW", "",
                                   "", "FBWadjMBW", "",
                                   "", "MBW", "",
                                   "", "MBWadjFBW", "",
                                   "", "Father", "")

table_base <- ggplot(grs_results_og, aes(y=alphabetical_id)) +
  ylab(NULL) + xlab("  ") + 
  theme(plot.title = element_text(hjust = 0.5, size=12), 
        axis.text.x = element_text(color="white", hjust = -3, size = 25), ## This is used to help with alignment
        axis.line = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.y = element_blank(), 
        legend.position = "none",
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.background = element_blank())

tab0 <- table_base + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = final_ids), size = 3, lineheight = .2,
            inherit.aes = TRUE, vjust = 1)


## 95% CI table
tab1 <- table_base +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = effects), size = 2.7,lineheight = 1,
            inherit.aes = TRUE,vjust = 1)

tab2 <- table_base +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = pvals), size = 2.7,lineheight = 1,
            inherit.aes = TRUE, vjust = 1)

lay <-  matrix(c(1,1,2,2,2,2,2,2,3,3,4,4,4,4,4,4,4,4,4,4,4,4), nrow = 1)

grid_1 <- grid.arrange(tab0, tab1, tab2,plot1, layout_matrix = lay)

##############################
#Now let's add the unweighted#
##############################

grs_results_og_2 <- grs_results_unw_is

grs_results_og_2$alphabetical_id <- c("a1", "a2", "a3", "a4", "a5",
                                    "b1", "b2", "b3", "b4", "b5",
                                    "c1", "c2", "c3", "c4", "c5")

plot2 <- ggplot(grs_results_og_2, aes(y = fct_rev(alphabetical_id), x = as.numeric(betas), shape = exposures)) +
  geom_point(size = 2.7, position = "identity") +  
  geom_errorbarh(aes(xmin = as.numeric(lower_cis), xmax = as.numeric(upper_cis)), height = 0.25) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  xlab("GRS effect size (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_blank(),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.ticks.y = element_blank(),
        legend.position = "right")



grs_results_og_2$final_ids <- c("", "FBW", "",
                              "", "FBWadjMBW", "",
                              "", "MBW", "",
                              "", "MBWadjFBW", "",
                              "", "Father", "")

table_base_2 <- ggplot(grs_results_og_2, aes(y=alphabetical_id)) +
  ylab(NULL) + xlab("  ") + 
  theme(plot.title = element_text(hjust = 0.5, size=12), 
        axis.text.x = element_text(color="white", hjust = -3, size = 25), ## This is used to help with alignment
        axis.line = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.y = element_blank(), 
        legend.position = "none",
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.background = element_blank())

tab0_2 <- table_base_2 + 
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = final_ids), size = 3, lineheight = .2,
            inherit.aes = TRUE, vjust = 1)


## 95% CI table
tab1_2 <- table_base_2 +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = effects), size = 2.7,lineheight = 1,
            inherit.aes = TRUE,vjust = 1)

tab2_2 <- table_base_2 +
  geom_text(aes(y = rev(alphabetical_id), x = 1, label = pvals), size = 2.7,lineheight = 1,
            inherit.aes = TRUE, vjust = 1)

lay_2 <-  matrix(c(1,1,2,2,2,2,2,2,3,3,4,4,4,4,4,4,4,4,4,4,4,4), nrow = 1)

grid_2 <- grid.arrange(tab0_2, tab1_2, tab2_2, plot2, layout_matrix = lay_2)

###################################
#And now let's merge the two grids#
###################################

tiff("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina_and_Mario/IR_BW_AIM4/Manuscript/figures/GRS/IS_GRS_plot.tiff", width = 6500, height=3200, res = 300)
grid.arrange(grid_1, grid_2, layout_matrix = matrix(c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2), nrow=1))
dev.off()
